"""Validated sample interfaces for the portable assembly workflow."""
import copy
import csv
import hashlib
import json
import re
from pathlib import Path

import jsonschema
import yaml
from assembly.scripts.prepared_input_utils import describe_prepared_input


class ConfigError(ValueError):
    pass


DEFAULTS = yaml.safe_load(Path(__file__).with_name('defaults.yaml').read_text())


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + '\n'


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def write_stable(path, value):
    path = Path(path)
    content = canonical(value)
    if path.exists() and path.read_text() == content:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(content)
    tmp.replace(path)


def merge(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def identifier(value, label):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value or ''):
        raise ConfigError(f'{label}: use a nonempty ASCII identifier (letters, digits, _, ., -)')
    return value


def table(path, required, allowed):
    try:
        with path.open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)) or not required <= set(fields):
                raise ConfigError(f'{path}: required unique columns: {sorted(required)}')
            if set(fields) - allowed:
                raise ConfigError(f'{path}: unknown columns: {sorted(set(fields) - allowed)}')
            result = []
            for line, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise ConfigError(f'{path}:{line}: malformed TSV row')
                if not any(v.strip() for v in row.values()):
                    continue
                result.append({k: v.strip() for k, v in row.items()})
            return result
    except OSError as exc:
        raise ConfigError(str(exc)) from exc


def source_signature(path):
    stat = Path(path).stat()
    # st_dev identifies a local mount and changes between HPC clients of the
    # same shared filesystem. It must not invalidate a cross-node resume.
    return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
            'ctime_ns': stat.st_ctime_ns, 'inode': stat.st_ino}


def resolve_config(config_path, repository, run_dir, selected, cores, mem_mb):
    """Metadata-only preflight; full checksums/FASTQ parsing belong in compute rules."""
    config_path, repository, run_dir = map(lambda p: Path(p).resolve(),
                                          (config_path, repository, run_dir))
    if cores < 1 or mem_mb < 1024:
        raise ConfigError('Execution budget requires cores >= 1 and mem_mb >= 1024')
    if run_dir == repository or run_dir.is_relative_to(repository):
        raise ConfigError('Run directory must be outside the workflow code repository')
    if (run_dir / 'layout_migration/status.json').exists():
        raise ConfigError('This completed run was migrated and is preserved as a historical result. '
                          'Read its numbered outputs; use a new --run-dir for new computation. '
                          'The original assembly must not be rerun just to update layout metadata.')
    previous = run_dir / 'config/plan.json'
    if previous.is_file():
        old_plan = json.loads(previous.read_text())
        if old_plan.get('layout_version', 1) < 2 and (run_dir / 'samples').is_dir():
            raise ConfigError('Legacy output layout detected. Migrate completed runs with '
                              'assembly/migrate_layout.py, or choose a new --run-dir. '
                              'No legacy outputs or Snakemake metadata were overwritten.')
    try:
        raw = yaml.safe_load(config_path.read_text())
        schema = yaml.safe_load((repository / 'assembly/schema.yaml').read_text())
        jsonschema.Draft202012Validator(schema).validate(raw)
    except (OSError, yaml.YAMLError, jsonschema.ValidationError) as exc:
        raise ConfigError(f'Configuration validation failed: {exc}') from exc
    sample_path = (config_path.parent / raw['samples_table']).resolve()
    readset_path = (config_path.parent / raw['readsets_table']).resolve()
    rows = table(sample_path, {'sample_id', 'species', 'profile'},
                 {'sample_id', 'species', 'profile'})
    samples = {}
    for row in rows:
        sid = identifier(row['sample_id'], 'sample_id')
        if sid in samples:
            raise ConfigError(f'Duplicate sample_id: {sid}')
        samples[sid] = row
    if not samples:
        raise ConfigError('samples.tsv is empty')
    unknown = set(raw.get('samples', {})) - set(samples)
    if unknown:
        raise ConfigError(f'Sample overrides refer to unknown samples: {sorted(unknown)}')
    selected = selected or list(samples)
    if len(selected) != len(set(selected)) or set(selected) - set(samples):
        raise ConfigError(f'Duplicate or unknown selected samples: {selected}')
    reads = table(readset_path,
                  {'sample_id', 'readset_id', 'library_id', 'technology', 'format', 'file1'},
                  {'sample_id', 'readset_id', 'library_id', 'technology', 'format',
                   'file1', 'file2', 'md5_1', 'md5_2', 'prepared_manifest'})
    by_sample, keys = {}, set()
    for row in reads:
        sid = row['sample_id']
        if sid not in samples:
            raise ConfigError(f'Readset refers to unknown sample: {sid}')
        rid = identifier(row['readset_id'], 'readset_id')
        identifier(row['library_id'], 'library_id')
        if (sid, rid) in keys:
            raise ConfigError(f'Duplicate readset_id for {sid}: {rid}')
        keys.add((sid, rid))
        if row['technology'] not in {'hifi', 'ont', 'hic', 'wgs'}:
            raise ConfigError(f'{sid}/{rid}: unknown technology')
        by_sample.setdefault(sid, []).append(row)
    result = {'schema_version': 1, 'layout_version': 2, 'config_path': str(config_path),
              'repository': str(repository), 'run_dir': str(run_dir),
              'budget': {'cores': cores, 'mem_mb': mem_mb}, 'sample_order': selected,
              'samples': {}}
    used_paths = set()
    for sid in selected:
        row = samples[sid]
        profile = row['profile']
        profiles = raw.get('analysis_profiles', {})
        if profile and profile not in profiles:
            raise ConfigError(f'{sid}: unknown analysis profile {profile}')
        settings = merge(DEFAULTS, raw.get('defaults', {}))
        settings = merge(settings, profiles.get(profile, {}))
        settings = merge(settings, raw.get('samples', {}).get(sid, {}))
        asm = settings['assembly']
        polyploid = (settings['ploidy'] == 4 and settings['scaffolding']['backend'] == 'gphase'
                     and asm['mode'] == 'hifi_ul_hic' and asm['representation'] == 'unitigs')
        if settings['ploidy'] is not None and settings['ploidy'] > 2 and not polyploid:
            raise ConfigError(f'{sid}: polyploid assembly strategy requires its validated backend; '
                              'use the explicit tetraploid GPhase unitig route; '
                              'the diploid backend defaults must not be applied silently')
        if settings['scaffolding']['backend'] == 'gphase' and not polyploid:
            raise ConfigError(f'{sid}: GPhase requires ploidy=4, hifi_ul_hic and representation=unitigs')
        joint = asm['mode'] in {'hifi_hic', 'hifi_ul_hic'}
        if not ((asm['mode'] == 'hifi_only' and asm['representation'] == 'primary_alternate')
                or (joint and asm['representation'] == 'phased') or polyploid):
            raise ConfigError(f'{sid}: supported routes: hifi_only/primary_alternate or '
                              'hifi_hic,hifi_ul_hic/phased (no automatic fallback)')
        scaffold = settings['scaffolding']
        if scaffold['hic_reference'] == 'per_haplotype' and not (
                settings['ploidy'] == 2 and joint and asm['representation'] == 'phased'
                and scaffold['backend'] == 'haphic'):
            raise ConfigError(f'{sid}: per_haplotype Hi-C mapping requires the diploid phased HapHiC route')
        if polyploid and asm['purge_level'] not in (None,0):
            raise ConfigError(f'{sid}: raw unitig route requires purge_level=0 or null')
        if polyploid and scaffold['gap_size'] != 100:
            raise ConfigError(f'{sid}: GPhase native gap policy currently requires gap_size=100')
        if polyploid and (scaffold['review_agp'] or scaffold['review_manifest']):
            raise ConfigError(f'{sid}: reviewed polyploid AGP import is not implemented')
        if polyploid and scaffold['correct_nrounds']:
            raise ConfigError(f'{sid}: GPhase uses its native YaHS correction; set correct_nrounds=0')
        if joint and ((settings['ploidy'] != 2 and not polyploid) or not settings['chromosomes_per_haplotype']):
            raise ConfigError(f'{sid}: joint phased scaffolding requires explicit ploidy=2 '
                              'and chromosomes_per_haplotype')
        if joint != scaffold['enabled']:
            raise ConfigError(f'{sid}: joint phased route requires scaffolding.enabled=true; '
                              'hifi_only retains its contig endpoint')
        if joint:
            if cores < 4:
                raise ConfigError(f'{sid}: joint mapping requires at least 4 cores')
            if mem_mb < 32000:
                raise ConfigError(f'{sid}: joint chromosome route requires a task budget >= 32000 MB')
            if scaffold['mem_mb'] > mem_mb:
                raise ConfigError(f'{sid}: scaffolding exceeds task budget')
            scaffold['threads'] = min(scaffold['threads'], cores)
            scaffold['processes'] = min(scaffold['processes'], scaffold['threads'])
            if bool(scaffold['review_agp']) != bool(scaffold['review_manifest']):
                raise ConfigError(f'{sid}: review_agp and review_manifest must be supplied together')
            for field in ('review_agp', 'review_manifest'):
                if scaffold[field]:
                    p = (config_path.parent / scaffold[field]).resolve()
                    if not p.is_file() or p.is_relative_to(run_dir):
                        raise ConfigError(f'{sid}: {field} must be an existing external review file')
                    scaffold[field] = str(p)
        if settings['genome_profile']['enabled']:
            raise ConfigError(f'{sid}: genome_profile is not implemented')
        finish, polish, assess = (settings[k] for k in ('finishing','polishing','assessment'))
        advanced = finish['enabled'] or polish['enabled'] or any(assess[k] for k in ('busco','kmer','read_mapping'))
        if advanced and not joint:
            raise ConfigError(f'{sid}: post-assembly modules require a chromosome route; '
                              'not implemented for the basic contig statistics only endpoint')
        if not finish['enabled'] and (finish['gap_closing'] or finish['telomere_completion'] or
                                     finish['source_delivery'] or polish['enabled'] or advanced):
            raise ConfigError(f'{sid}: enable finishing for the post-assembly assessment route')
        if finish['enabled']:
            if any(c.isspace() for c in str(run_dir)):
                raise ConfigError('Post-assembly native shell tools require a whitespace-free run directory')
            if not finish['telomere_motif']:
                raise ConfigError(f'{sid}: finishing requires an explicit telomere_motif')
            if finish['local_window'] < finish['min_anchor']:
                raise ConfigError(f'{sid}: local_window must be >= min_anchor')
            for group in (finish,polish,assess):
                if group['mem_mb'] > mem_mb:
                    raise ConfigError(f'{sid}: post-assembly module exceeds task budget')
                group['threads'] = min(group['threads'],cores)
            if assess['terminal_max_distance'] >= assess['terminal_window']:
                raise ConfigError(f'{sid}: terminal_max_distance must be smaller than terminal_window')
            for group,field in [(finish,'source_delivery'),(assess,'busco_lineage')]:
                if group[field]:
                    p = (config_path.parent/group[field]).resolve()
                    if not p.exists() or p.is_relative_to(run_dir):
                        raise ConfigError(f'{sid}: {field} must be an existing external path')
                    group[field] = str(p)
            if finish['source_delivery'] and not Path(finish['source_delivery']).is_file():
                raise ConfigError(f'{sid}: source_delivery must point to manifest.json')
            if assess['busco'] and (not assess['busco_lineage'] or
                                   not (Path(assess['busco_lineage'])/'dataset.cfg').is_file()):
                raise ConfigError(f'{sid}: BUSCO needs an existing offline lineage directory with dataset.cfg')
        if asm['mem_mb'] > mem_mb:
            raise ConfigError(f'{sid}: assembly requires {asm["mem_mb"]} MB; '
                              f'task budget is {mem_mb} MB')
        # Bloom storage alone is 2**f bits; retain at least 512 MB for other state.
        bloom_mb = (2 ** asm['bloom_filter_bits']) / 8 / 1_000_000 if asm['bloom_filter_bits'] else 0
        if bloom_mb + 512 > asm['mem_mb']:
            raise ConfigError(f'{sid}: bloom_filter_bits exceeds the assembly memory budget')
        asm['threads'] = min(asm['threads'], cores)
        active = {}
        technologies = {'hifi', 'hic'} | ({'ont'} if asm['mode'] == 'hifi_ul_hic' else set()) if joint else {'hifi'}
        required_technologies = set(technologies)
        if finish['enabled'] and (polish['enabled'] or assess['kmer']):
            technologies.add('wgs')
        for read in by_sample.get(sid, []):
            if read['technology'] not in technologies:
                continue
            rid = read['readset_id']
            if read['technology'] in {'hic','wgs'}:
                if read['format'] != 'fastq' or not read.get('file2'):
                    raise ConfigError(f'{sid}/{rid}: Hi-C/WGS requires paired FASTQ file1 and file2')
            elif read['format'] not in {'fastq', 'bam'} or read.get('file2') or (read['technology'] == 'ont' and read['format'] != 'fastq'):
                raise ConfigError(f'{sid}/{rid}: HiFi input must be single-read FASTQ or unpaired BAM')
            path = (readset_path.parent / read['file1']).resolve()
            if not path.is_file() or path.stat().st_size == 0:
                raise ConfigError(f'{sid}/{rid}: missing or empty input {path}')
            if path.is_relative_to(run_dir):
                raise ConfigError(f'{sid}/{rid}: source input must be outside the run directory')
            if path in used_paths:
                raise ConfigError(f'Duplicate active source path (would duplicate reads): {path}')
            used_paths.add(path)
            md5 = read.get('md5_1', '').lower()
            if md5 and not re.fullmatch('[0-9a-f]{32}', md5):
                raise ConfigError(f'{sid}/{rid}: invalid md5_1')
            active[rid] = {**read, 'file1': str(path), 'expected_md5': md5 or None,
                           'source_signature': source_signature(path)}
            if read.get('prepared_manifest'):
                if not raw.get('prepared_input_cache'):
                    raise ConfigError(f'{sid}/{rid}: prepared_manifest requires an external prepared_input_cache')
                try:
                    manifest=(readset_path.parent/read['prepared_manifest']).resolve()
                    cache=(config_path.parent/raw['prepared_input_cache']).resolve()
                    active[rid]['prepared_manifest']=str(manifest)
                    active[rid]['prepared_reuse']=describe_prepared_input(active[rid],manifest,cache,run_dir,repository)
                except (OSError,ValueError,KeyError,TypeError) as exc:
                    raise ConfigError(f'{sid}/{rid}: {exc}') from exc
            if read['technology'] in {'hic','wgs'}:
                mate = (readset_path.parent / read['file2']).resolve()
                md5_2 = read.get('md5_2', '').lower()
                if not mate.is_file() or not mate.stat().st_size or mate.is_relative_to(run_dir):
                    raise ConfigError(f'{sid}/{rid}: missing or empty external mate {mate}')
                if mate in used_paths:
                    raise ConfigError(f'Duplicate active source path (would duplicate reads): {mate}')
                if md5_2 and not re.fullmatch('[0-9a-f]{32}', md5_2):
                    raise ConfigError(f'{sid}/{rid}: invalid md5_2')
                used_paths.add(mate)
                active[rid].update(file2=str(mate), expected_md5_2=md5_2 or None,
                                   source_signature_2=source_signature(mate))
        if not active:
            raise ConfigError(f'{sid}: hifi_only requires at least one HiFi readset')
        missing = required_technologies - {r['technology'] for r in active.values()}
        if missing:
            raise ConfigError(f'{sid}: required read technologies missing: {sorted(missing)}')
        if finish['enabled'] and polish['enabled'] and polish['on_missing_wgs']=='error' and not any(
                r['technology']=='wgs' for r in active.values()):
            raise ConfigError(f'{sid}: NextPolish2 requires ordinary WGS short reads, not Hi-C')
        physical_inputs={}
        for rid,spec in active.items():
            paths={spec['file1']}
            if spec.get('file2'):paths.add(spec['file2'])
            if spec.get('prepared_reuse'):paths.add(spec['prepared_reuse']['identity']['prepared_fastq'])
            for path in paths:
                if path in physical_inputs and physical_inputs[path]!=rid:
                    raise ConfigError(f'{sid}: duplicate raw/prepared source path across readsets')
                physical_inputs[path]=rid
        hic_libraries = [r['library_id'] for r in active.values() if r['technology'] == 'hic']
        if len(hic_libraries) != len(set(hic_libraries)):
            raise ConfigError(f'{sid}: combine lanes of each Hi-C library before registering; '
                              'duplicate marking must cover the complete physical library')
        result['samples'][sid] = {'sample_id': sid, 'species': row['species'],
                                  'profile': profile, 'settings': settings, 'readsets': active,
                                  'unused_readsets': [r['readset_id'] for r in by_sample.get(sid, [])
                                                     if r['technology'] not in technologies]}
    routes = {s['settings']['scaffolding']['backend'] == 'gphase' for s in result['samples'].values()}
    if len(routes) > 1:
        raise ConfigError('Diploid and polyploid routes require separate run directories')
    if routes == {True}:
        if any(c.isspace() for c in str(run_dir)):
            raise ConfigError('GPhase upstream shell scripts require a whitespace-free run directory')
        result['workflow_route'] = 'polyploid'
    return result
