#!/usr/bin/env python3
"""Run the complete assembly DAG locally; an HPC wrapper may invoke this once."""
import argparse
import ast
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True


def main(argv=None, prog=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__, prog=prog)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--run-dir', required=True, type=Path)
    parser.add_argument('--samples', nargs='+', help='Default: all samples in samples.tsv')
    parser.add_argument('--cores', required=True, type=int)
    parser.add_argument('--mem-mb', required=True, type=int)
    parser.add_argument('--target', choices=['contigs', 'chromosomes', 't2t', 'report'], default='report')
    parser.add_argument('--env-prefix', type=Path, help='Optional unified environment prefix')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check', action='store_true', help='Metadata/dependency preflight only')
    args = parser.parse_args(argv)
    repository = Path(__file__).resolve().parents[1]
    if args.env_prefix:
        prefix = args.env_prefix.resolve()
        executable = prefix / 'bin/python'
        if not executable.is_file():
            parser.error(f'Environment Python not found: {executable}')
        env = dict(os.environ, PATH=str(prefix / 'bin') + os.pathsep + os.environ.get('PATH', ''),
                   CONDA_PREFIX=str(prefix), PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
        if Path(sys.prefix).resolve() != prefix or os.environ.get('CONDA_PREFIX') != str(prefix):
            os.execve(str(executable), [str(executable), str(repository / 'run'), 'assembly', *argv], env)
        os.environ.update(env)
    sys.path.insert(0, str(repository))
    import importlib.metadata
    from assembly.configuration import ConfigError, canonical, digest, resolve_config, write_stable
    try:
        plan = resolve_config(args.config, repository, args.run_dir, args.samples,
                              args.cores, args.mem_mb)
        tools = {}
        for name, command in {'hifiasm': ['--version'], 'samtools': ['--version'],
                              'seqkit': ['version']}.items():
            path = shutil.which(name)
            if not path or not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve()):
                raise ConfigError(f'{name} must be in the unified active environment; '
                                  'activate it or supply --env-prefix')
            version = subprocess.check_output([path, *command], stderr=subprocess.STDOUT,
                                              timeout=30).decode(errors='replace').splitlines()[0]
            tools[name] = {'path': path, 'version': version}
        snakemake = Path(sys.prefix) / 'bin/snakemake'
        if not snakemake.is_file():
            raise ConfigError('snakemake missing from the unified active environment')
        tools['snakemake'] = {'path': str(snakemake),
                              'version': importlib.metadata.version('snakemake')}
        joint = any(s['settings']['scaffolding']['enabled'] for s in plan['samples'].values())
        polyploid = plan.get('workflow_route') == 'polyploid'
        postassembly = args.target in ('t2t','report') and any(
            s['settings']['finishing']['enabled'] for s in plan['samples'].values())
        if args.target=='t2t' and not all(s['settings']['finishing']['enabled'] for s in plan['samples'].values()):
            raise ConfigError('--target t2t requires finishing.enabled=true for every selected sample')
        if args.target == 'chromosomes' and not all(s['settings']['scaffolding']['enabled']
                                                    for s in plan['samples'].values()):
            raise ConfigError('--target chromosomes requires a scaffolding-enabled route for every selected sample')
        if joint:
            import importlib
            for module in ('pysam', 'scipy', 'sklearn', 'portion', 'networkx', 'matplotlib'):
                importlib.import_module(module)
            for name, options in {'bwa': [], 'samblaster': ['--help'], 'haphic': ['--version']}.items():
                path = shutil.which(name)
                if not path or not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve()):
                    raise ConfigError(f'{name} missing from the unified environment')
                output = subprocess.run([path, *options], stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, timeout=30).stdout.decode(errors='replace')
                tools[name] = {'path': path, 'version': output.strip()[:2000]}
            haphic_root = Path(tools['haphic']['path']).resolve().parent
            for name in ('filter_bam', 'juicer'):
                path = haphic_root / 'utils' / name
                if not path.is_file() or not os.access(path, os.X_OK):
                    raise ConfigError(f'HapHiC bundled utility missing: {path}')
                tools[name] = {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
            tools['haphic']['source_lock'] = json.loads((haphic_root / 'SOURCE_LOCK.json').read_text())
        if polyploid:
            addon = Path(sys.prefix) / 'opt/GPhase'
            # Added modules live under the same Python prefix; existing packages stay untouched.
            sys.path.insert(0, str(addon / 'python'))
            for module in ('statsmodels.api','igraph','kneed','community','pathos','argcomplete','networkit'):
                importlib.import_module(module)
            os.environ['PYTHONPATH'] = str(addon / 'python') + os.pathsep + os.environ.get('PYTHONPATH','')
            os.environ['PATH'] = str(addon/'bin') + os.pathsep + str(addon/'python/bin') + os.pathsep + os.environ['PATH']
            for name in ('gphase','mosdepth','yahs'):
                path = shutil.which(name)
                if not path or not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve()):
                    raise ConfigError(f'{name} missing from unified environment')
                tools[name] = {'path':path,'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest()}
            tools['gphase']['source_lock'] = json.loads((addon/'SOURCE_LOCK.json').read_text())
            for relative, expected in tools['gphase']['source_lock']['source_files'].items():
                if hashlib.sha256((addon/'source'/relative).read_bytes()).hexdigest() != expected:
                    raise ConfigError(f'Locked GPhase source changed: {relative}')
        if postassembly:
            addon=Path(sys.prefix)/'opt/T2T'
            lock=addon/'SOURCE_LOCK.json'
            if not lock.is_file():raise ConfigError('Install the locked optional T2T tools in the unified prefix first')
            locked=json.loads(lock.read_text())
            for relative,expected in locked['files'].items():
                if hashlib.sha256((addon/relative).read_bytes()).hexdigest()!=expected:
                    raise ConfigError('Locked T2T tool changed: '+relative)
            required={'tidk','minimap2'}
            for sample in plan['samples'].values():
                settings=sample['settings']
                if settings['finishing']['gap_closing']:required.add('quartet')
                if settings['polishing']['enabled'] and any(r['technology']=='wgs' for r in sample['readsets'].values()):
                    required.update(['nextPolish2','yak','fastp'])
                if settings['assessment']['read_mapping']:required.add('inspector.py')
                if settings['assessment']['kmer']:required.update(['merqury.sh','meryl','fastp'])
                if settings['assessment']['busco']:required.add('busco')
            for name in sorted(required):
                path=shutil.which(name)
                if not path or not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve()):
                    raise ConfigError(name+' missing from unified environment')
                tools[name]={'path':path,'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                             't2t_source_lock_sha256':hashlib.sha256(lock.read_bytes()).hexdigest()}
    except (ConfigError, OSError, ImportError, subprocess.SubprocessError) as exc:
        parser.exit(2, f'PRECHECK FAILED: {exc}\n')
    code_paths = [repository / 'run', *sorted((repository / 'common').rglob('*.py')),
                  *sorted((repository / 'assembly').rglob('*.py')),
                  *sorted((repository / 'assembly').rglob('*.smk')), repository / 'assembly/Snakefile',
                  repository / 'assembly/schema.yaml', repository / 'assembly/defaults.yaml']
    plan['code'] = {str(p.relative_to(repository)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in code_paths}

    def stage_code(stage, *scripts):
        """Keep full provenance in the plan; invalidate only a stage's actual code."""
        selected = {f'assembly/rules/{stage}.smk'}
        pending = list(scripts)
        while pending:
            name = pending.pop()
            relative = f'assembly/scripts/{name}.py'
            if relative in selected:
                continue
            selected.add(relative)
            for node in ast.walk(ast.parse((repository / relative).read_text())):
                if isinstance(node, ast.ImportFrom) and node.module:
                    if (repository / f'assembly/scripts/{node.module}.py').is_file():
                        pending.append(node.module)
        return {key: plan['code'][key] for key in sorted(selected)}

    stage_sources = {
        'joint': stage_code('04_assembly', 'assemble_phased_contigs'),
        'linkage': stage_code('05_linkage', 'map_filter_hic', 'merge_validate_hic_bams'),
        'scaffolding': stage_code('06_scaffolding', 'scaffold_haplotypes'),
        'delivery': stage_code('06_scaffolding', 'export_chromosome_fasta'),
        'plot': stage_code('06_scaffolding', 'plot_hic_contacts'),
    }
    plan['tools'] = tools
    plan['target'] = args.target
    plan['plan_id'] = digest(plan)
    print(f'PRECHECK PASS: samples={",".join(plan["sample_order"])} '
          f'cores={args.cores} mem_mb={args.mem_mb} '
          f'modes={",".join(plan["samples"][s]["settings"]["assembly"]["mode"] for s in plan["sample_order"])} '
          'executor=local', flush=True)
    if args.check:
        return 0
    run_dir = Path(plan['run_dir'])
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / '.launcher.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.exit(2, f'Run directory already has an active launcher: {run_dir}\n')
        for sid, sample in plan['samples'].items():
            for rid, read in sample['readsets'].items():
                if read['technology']=='wgs':continue
                script = 'validate_hic_pairs' if read['technology'] == 'hic' else 'prepare_long_reads'
                if read.get('prepared_reuse'):script='reuse_prepared_reads'
                spec = {**read, 'samtools': tools['samtools'], 'code': stage_code('01_inputs', script)}
                if read.get('prepared_reuse'):spec['seqkit']=tools['seqkit']
                write_stable(run_dir / f'config/readsets/{sid}/{rid}.json', spec)
            write_stable(run_dir / f'config/samples/{sid}.assembly.json',
                         {'sample_id': sid, 'assembly': sample['settings']['assembly'],
                          'readsets': list(sample['readsets']), 'hifiasm': tools['hifiasm'],
                          'code': stage_code('04_assembly', 'assemble_hifi_contigs')})
            write_stable(run_dir / f'config/samples/{sid}.metadata.json',
                         {**sample, 'tools': tools, 'code': plan['code']})
            if sample['settings']['scaffolding']['enabled']:
                # Keep review-only changes downstream of expensive assembly and mapping.
                stages = {
                    'joint': {'assembly': sample['settings']['assembly']},
                    'linkage': {'scaffolding': {k: v for k, v in sample['settings']['scaffolding'].items()
                                               if k in ('threads', 'mapq', 'nm_exclusive', 'mem_mb',
                                                        'hic_reference', 'min_retained_fraction')}},
                    'scaffolding': {'chromosomes_per_haplotype': sample['settings']['chromosomes_per_haplotype'],
                                    'scaffolding': {k: v for k, v in sample['settings']['scaffolding'].items()
                                                   if k not in ('review_agp', 'review_manifest', 'plot_bin_kb')}},
                    'delivery': {'chromosomes_per_haplotype': sample['settings']['chromosomes_per_haplotype'],
                                 'scaffolding': {k: v for k, v in sample['settings']['scaffolding'].items()
                                                if k in ('review_agp', 'review_manifest')}},
                    'plot': {'scaffolding': {'plot_bin_kb': sample['settings']['scaffolding']['plot_bin_kb']}},
                }
                stage_tools = {
                    'joint': ('hifiasm', 'samtools'),
                    'linkage': ('bwa', 'samblaster', 'samtools', 'filter_bam'),
                    'scaffolding': ('haphic',),
                    'delivery': ('samtools',),
                    'plot': ('haphic',),
                }
                for stage, settings in stages.items():
                    write_stable(run_dir / f'config/samples/{sid}.{stage}.json',
                                 {'sample_id': sid, 'settings': settings,
                                  'tools': {key: tools[key] for key in stage_tools[stage]},
                                  'code': stage_sources[stage]})
            if polyploid:
                stages = {
                    'poly_assembly': ('04_assembly', ['assemble_polyploid_unitigs'],
                                      {'assembly': sample['settings']['assembly'], 'ploidy': sample['settings']['ploidy']}),
                    'poly_scaffolding': ('06_scaffolding', ['scaffold_polyploid_gphase'],
                        {'ploidy':sample['settings']['ploidy'],
                         'chromosomes_per_haplotype':sample['settings']['chromosomes_per_haplotype'],
                         'scaffolding':{k:sample['settings']['scaffolding'][k]
                                        for k in ('restriction_sites','processes','threads','mem_mb')}}),
                    'poly_delivery': ('06_scaffolding', ['export_polyploid_fasta'],
                        {k:sample['settings'][k] for k in ('ploidy','chromosomes_per_haplotype')}),
                    'poly_plot': ('06_scaffolding', ['plot_polyploid_contacts'],
                                  {'scaffolding': {k:sample['settings']['scaffolding'][k]
                                                   for k in ('plot_bin_kb','mapq','nm_exclusive')}}),
                }
                poly_tools = {'poly_assembly':('hifiasm','samtools'),
                              'poly_scaffolding':('gphase','yahs','samtools'),
                              'poly_delivery':('samtools',),
                              'poly_plot':('bwa','samblaster','samtools','filter_bam','haphic')}
                for stage, (number, scripts, settings) in stages.items():
                    write_stable(run_dir / f'config/samples/{sid}.{stage}.json',
                                 {'sample_id':sid,'settings':settings,
                                  'tools':{key:tools[key] for key in poly_tools[stage]},'code':stage_code(number,*scripts)})
            if postassembly and sample['settings']['finishing']['enabled']:
                settings=sample['settings']
                terminal={k:settings['assessment'][k] for k in ('terminal_window','terminal_max_distance','terminal_min_repeats')}
                stages={
                    'finishing':('07_finishing','finish_chromosomes',
                        {'finishing':settings['finishing'],'assessment':terminal,
                         'ploidy':settings['ploidy'],'chromosomes_per_haplotype':settings['chromosomes_per_haplotype']},
                        ('hifi','ont'),('hifiasm','minimap2','samtools','quartet','tidk')),
                    'polishing':('08_polishing','polish_chromosomes',{'polishing':settings['polishing']},
                        ('wgs',),('nextPolish2','minimap2','samtools','yak','fastp')),
                    't2t_assessment':('09_assessment','assess_t2t',
                        {'assessment':settings['assessment'],'finishing':{'telomere_motif':settings['finishing']['telomere_motif']},
                         'ploidy':settings['ploidy'],'chromosomes_per_haplotype':settings['chromosomes_per_haplotype']},
                        ('wgs',),('tidk','minimap2','samtools','inspector.py','merqury.sh','meryl','fastp','busco'))}
                for name,(number,script,stage_settings,technologies,selected_tools) in stages.items():
                    write_stable(run_dir/f'config/samples/{sid}.{name}.json',
                        {'sample_id':sid,'settings':stage_settings,
                         'readsets':{r:spec for r,spec in sample['readsets'].items() if spec['technology'] in technologies},
                         'tools':{t:tools[t] for t in selected_tools if t in tools},'code':stage_code(number,script)})
        write_stable(run_dir / 'config/plan.json', plan)
        command = [str(snakemake), '--snakefile', str(repository / 'assembly/Snakefile'),
                   '--directory', str(run_dir), '--configfile', str(run_dir / 'config/plan.json'),
                   '--executor', 'local', '--cores', str(args.cores),
                   '--resources', f'mem_mb={args.mem_mb}', '--rerun-incomplete',
                   '--printshellcmds', 'assembly_' + args.target]
        if args.dry_run:
            command.append('--dry-run')
        status = {'plan_id': plan['plan_id'], 'samples': plan['sample_order'],
                  'target': args.target, 'executor': 'local', 'command': command,
                  'started_utc': datetime.now(timezone.utc).isoformat()}
        status_path = run_dir / ('dryrun_status.json' if args.dry_run else 'run_status.json')
        write_stable(status_path, {**status, 'status': 'PLANNING' if args.dry_run else 'RUNNING'})
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONNOUSERSITE='1',
                   PYTHONDONTWRITEBYTECODE='1')
        handlers, interrupted = {}, []
        try:
            process = subprocess.Popen(command, env=env, start_new_session=True)

            def forward_signal(signum, frame):
                interrupted.append(signum)
                try:
                    os.killpg(process.pid, signum)
                except ProcessLookupError:
                    pass

            for signum in (signal.SIGINT, signal.SIGTERM):
                handlers[signum] = signal.signal(signum, forward_signal)
            interruption_started = None
            while True:
                try:
                    rc = process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    if interrupted:
                        if interruption_started is None:
                            interruption_started = time.monotonic()
                        elif time.monotonic() - interruption_started > 15:
                            # Some Snakemake worker shutdown races otherwise leave
                            # the outer allocation and launcher lock alive forever.
                            import psutil
                            try:
                                descendants = psutil.Process(process.pid).children(recursive=True)
                            except psutil.NoSuchProcess:
                                descendants = []
                            for child in reversed(descendants):
                                try:
                                    child.kill()
                                except psutil.NoSuchProcess:
                                    pass
                            try:
                                os.killpg(process.pid, signal.SIGKILL)
                            except ProcessLookupError:
                                pass
                            rc = process.wait()
                            break
            if interrupted:
                rc = 128 + interrupted[-1]
        except OSError:
            rc = 1
            raise
        finally:
            for signum, handler in handlers.items():
                signal.signal(signum, handler)
            write_stable(status_path, {**status, 'exit_code': locals().get('rc', 1),
                                       'status': ('PLANNED' if args.dry_run else 'COMPLETED')
                                       if locals().get('rc', 1) == 0 else 'FAILED',
                                       'finished_utc': datetime.now(timezone.utc).isoformat()})
        return rc


if __name__ == '__main__':
    raise SystemExit(main())
