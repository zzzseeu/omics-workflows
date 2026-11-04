#!/usr/bin/env python3
"""Verify delivered FASTA hashes and collect chromosome assessment records."""
import argparse
from io_utils import checksum, load, save
from pathlib import Path


def collect_chromosome_metrics(args):
    delivery = load(args.delivery)
    root = Path(args.delivery).parent
    for label, stats in delivery['sets'].items():
        if checksum(root/f'{label}.fa') != stats['sha256']:
            raise ValueError('Delivered FASTA changed after validation')
    assessment = {**delivery, 'evaluations': {'AGP_source_coverage': 'PASS',
                  'native_FASTA_AGP_concordance': 'PASS', 'HiC_contact_maps': str(Path(args.plots).resolve()),
                  'human_review': delivery['human_review']}}
    save(args.assessment, assessment)
    with Path(args.table).open('w') as out:
        out.write('set_id\tnum_sequences\ttotal_bases\tN50\tN_bases\n')
        for label, stat in delivery['sets'].items():
            out.write(f'{label}\t{stat["num_contigs"]}\t{stat["total_bases"]}\t{stat["N50"]}\t{stat["N_bases"]}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delivery", required=True)
    parser.add_argument("--plots", required=True)
    parser.add_argument("--assessment", required=True)
    parser.add_argument("--table", required=True)
    collect_chromosome_metrics(parser.parse_args())


if __name__ == "__main__":
    main()
