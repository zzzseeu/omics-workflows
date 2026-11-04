#!/usr/bin/env python3
"""Assemble phased contigs, prefix IDs and preserve source provenance."""
import argparse
import hashlib
import json
import shutil
from fasta_utils import gfa_to_fasta, prefixed_gfa
from io_utils import assembly_scratch, checksum, execute, load, publish, save, source_check
from pathlib import Path


def joint_assemble(args):
    metadata = load(args.metadata)
    options = metadata['settings']['assembly']
    for q in args.qcs:
        if load(q)['status'] != 'PASS':
            raise ValueError('Read QC did not pass')
    hic = [load(p) for p in args.paired]
    for item in hic:
        source_check(item['source_spec'])
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temp:
        stage = Path(temp) / 'assembly'
        stage.mkdir()
        command = ['hifiasm', '-o', stage / 'asm', '-t', args.threads,
                   '-f', options['bloom_filter_bits'],
                   '--h1', ','.join(x['source_spec']['file1'] for x in hic),
                   '--h2', ','.join(x['source_spec']['file2'] for x in hic)]
        if args.ul:
            command += ['--ul', ','.join(args.ul)]
        if options['purge_level'] is not None:
            command += ['-l', options['purge_level']]
        command += args.hifi
        execute(command)
        sets = []
        with (stage / 'combined.fa').open('wb') as combined, (stage / 'contig_haplotypes.tsv').open('w') as mapping:
            mapping.write('contig_id\thaplotype\toriginal_id\n')
            for hap in ('hap1', 'hap2'):
                source = stage / f'asm.hic.{hap}.p_ctg.gfa'
                if not source.is_file():
                    raise ValueError(f'Missing phased hifiasm output: {source}')
                target = stage / f'{hap}.gfa'
                prefixed_gfa(source, target, hap + '_')
                fasta = stage / f'{hap}.fa'
                count = gfa_to_fasta(target, fasta, True)
                with fasta.open('rb') as f:
                    shutil.copyfileobj(f, combined)
                with fasta.open() as f:
                    for line in f:
                        if line.startswith('>'):
                            cid = line[1:].strip()
                            mapping.write(f'{cid}\t{hap}\t{cid[len(hap)+1:]}\n')
                sets.append({'set_id': hap, 'fasta': fasta.name, 'normalized_gfa': target.name,
                             'source_gfa': source.name, 'num_contigs': count})
        execute(['samtools', 'faidx', stage / 'combined.fa'])
        for item in hic:
            source_check(item['source_spec'])
        inputs = [load(p) for p in args.manifests] + hic
        identity = {'options': options, 'inputs': inputs, 'tools': metadata['tools'], 'code': metadata['code']}
        assembly_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(stage / 'manifest.json', {'sample_id': metadata['sample_id'], 'assembly_id': assembly_id,
             'mode': options['mode'], 'representation': 'phased', 'phase_level': 'HIC_PHASED_UNVALIDATED',
             'sets': sets, 'command': [str(x) for x in command], 'metadata': metadata,
             'input_manifests': inputs, 'combined_sha256': checksum(stage / 'combined.fa'),
             'limitations': ['Haplotype labels are not parental identities; chromosome-wide phase switches remain possible']})
        publish(stage, destination, temp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--hifi", nargs="+", required=True)
    parser.add_argument("--ul", nargs="*", required=False, default=[])
    parser.add_argument("--paired", nargs="+", required=True)
    parser.add_argument("--qcs", nargs="+", required=True)
    parser.add_argument("--manifests", nargs="+", required=True)
    joint_assemble(parser.parse_args())


if __name__ == "__main__":
    main()
