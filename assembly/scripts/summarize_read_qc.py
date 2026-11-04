#!/usr/bin/env python3
"""Validate seqkit statistics and write the read QC decision."""
import argparse
import csv
from io_utils import save
from pathlib import Path


def read_qc(args):
    with Path(args.stats).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 1 or rows[0]['format'].upper() != 'FASTQ':
        raise ValueError('Expected one FASTQ statistics row')
    row = rows[0]
    count = int(row['num_seqs'].replace(',', ''))
    bases = int(row['sum_len'].replace(',', ''))
    if count < 1 or bases < 1:
        raise ValueError('Empty FASTQ is not a valid assembly input')
    save(args.json, {'status': 'PASS', 'validation': 'seqkit full-file parsing and statistics',
                     'num_reads': count, 'num_bases': bases, 'statistics': row,
                     'limitations': ['No cross-library duplicate-read-ID audit or contamination screen']})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--json", required=True)
    read_qc(parser.parse_args())


if __name__ == "__main__":
    main()
