#!/usr/bin/env python3
"""Collect contig statistics and explicit unperformed evaluation states."""
import argparse
import csv
from fasta_utils import fasta_statistics
from io_utils import load, save
from pathlib import Path


def assess(args):
    assembly = Path(args.assembly_dir)
    manifest = load(assembly / 'manifest.json')
    stats = {s['set_id']: fasta_statistics(assembly / s['fasta']) for s in manifest['sets']}
    qcs = [load(p) for p in args.qcs]
    disabled = {key: {'status': 'NOT_RUN', 'reason': 'Not implemented in this MVP'}
                for key in ['BUSCO', 'kmer_QV', 'read_mapping', 'structural_validation',
                            'telomere_validation', 'phasing_validation', 'contamination_screen']}
    result = {'sample_id': manifest['sample_id'], 'assembly_id': manifest['assembly_id'],
              'quality_level': 'CONTIG_ONLY_UNVALIDATED', 'sets': stats,
              'input_reads': sum(q['num_reads'] for q in qcs),
              'input_bases': sum(q['num_bases'] for q in qcs), 'evaluations': disabled,
              'T2T_status': 'NOT_ASSESSED'}
    Path(args.table).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.table).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['set_id', 'num_contigs', 'total_bases', 'N50', 'L50', 'longest_contig', 'N_bases'])
        for label, values in stats.items():
            writer.writerow([label, *(values[k] for k in
                                     ['num_contigs', 'total_bases', 'N50', 'L50', 'longest_contig', 'N_bases'])])
    save(args.json, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--qcs", nargs="+", required=True)
    assess(parser.parse_args())


if __name__ == "__main__":
    main()
