#!/usr/bin/env python3
"""Run HiFi assembly transactionally and normalize GFA/FASTA outputs."""
import argparse
import hashlib
import json
import shutil
import tempfile
from fasta_utils import gfa_to_fasta
from io_utils import execute, load, save
from pathlib import Path


def assemble(args):
    spec = load(args.spec)
    options = spec['assembly']
    qcs = [load(p) for p in args.qcs]
    if any(q['status'] != 'PASS' for q in qcs):
        raise ValueError('Read QC did not pass')
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Every recomputation gets a fresh cache; old hifiasm .bin files are never silently reused.
    with tempfile.TemporaryDirectory(prefix='.hifiasm-', dir=destination.parent) as temp:
        stage = Path(temp) / 'assembly'
        stage.mkdir()
        prefix = stage / 'asm'
        command = ['hifiasm', '-o', str(prefix), '-t', str(args.threads), '--primary',
                   '-f', str(options['bloom_filter_bits'])]
        if options['purge_level'] is not None:
            command += ['-l', str(options['purge_level'])]
        command += args.fastqs
        execute(command)
        print('HIFIASM_OUTPUTS', json.dumps(sorted(p.name for p in stage.iterdir())), flush=True)
        sets = []
        for label, suffix in [('primary', 'p_ctg'), ('alternate', 'a_ctg')]:
            candidates = [stage / f'asm.bp.{suffix}.gfa', stage / f'asm.{suffix}.gfa']
            existing = [p for p in candidates if p.is_file()]
            if len(existing) > 1:
                raise ValueError(f'Ambiguous hifiasm {label} GFA output: {existing}')
            if label == 'primary' and not existing:
                raise ValueError(f'Missing primary GFA; expected one of {candidates}')
            gfa = existing[0] if existing else candidates[0]
            normalized_gfa = stage / f'{label}.gfa'
            if gfa.exists():
                shutil.copyfile(gfa, normalized_gfa)
            else:
                normalized_gfa.write_text('')
            fasta = stage / f'{label}.fa'
            count = gfa_to_fasta(gfa, fasta, label == 'primary')
            sets.append({'set_id': label, 'representation': label, 'num_contigs': count,
                         'has_sequence': bool(count), 'fasta': fasta.name,
                         'fai': fasta.name + '.fai', 'normalized_gfa': normalized_gfa.name,
                         'source_gfa': gfa.name if gfa.exists() else None})
        inputs = [load(p) for p in args.manifests]
        identity = {'spec': spec, 'input_sha256': [p['source_sha256'] for p in inputs]}
        assembly_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(stage / 'manifest.json', {'sample_id': spec['sample_id'], 'assembly_id': assembly_id,
                                      'mode': 'hifi_only', 'representation': 'primary_alternate',
                                      'phase_level': 'NOT_CERTIFIED', 'sets': sets,
                                      'command': command, 'settings': options,
                                      'input_manifests': inputs, 'tools': {'hifiasm': spec['hifiasm']},
                                      'code': spec['code'], 'cache_policy': 'fresh_on_recompute'})
        (stage / '.omics-managed').write_text('omics-workflows assembly output\n')
        if destination.exists():
            # Snakemake creates output parents before executing a file-output rule.
            if any(destination.iterdir()) and not (destination / '.omics-managed').is_file():
                raise ValueError(f'Refusing to replace unmanaged assembly directory: {destination}')
            destination.replace(Path(temp) / 'previous')
        stage.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--fastqs", nargs="+", required=True)
    parser.add_argument("--manifests", nargs="+", required=True)
    parser.add_argument("--qcs", nargs="+", required=True)
    assemble(parser.parse_args())


if __name__ == "__main__":
    main()
