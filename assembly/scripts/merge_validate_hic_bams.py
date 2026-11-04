#!/usr/bin/env python3
"""Merge Hi-C libraries and verify adjacent pairs and count conservation."""
import argparse
from hic_bam_utils import bam_pairs
from io_utils import checksum, execute, load, save
from pathlib import Path


def merge_hic(args):
    # Concatenation preserves adjacent pairs. Name sorting across libraries is unnecessary.
    out = Path(args.bam)
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_suffix('.partial.bam')
    try:
        import pysam
        headers = []
        for path in args.inputs:
            with pysam.AlignmentFile(str(path), 'rb') as bam:
                headers.append(list(zip(bam.references, bam.lengths)))
        if any(h != headers[0] for h in headers[1:]):
            raise ValueError('Hi-C library BAM reference dictionaries differ; refusing unsafe concatenation')
        execute(['samtools', 'cat', '-o', temporary, *args.inputs])
        execute(['samtools', 'quickcheck', '-v', temporary])
        expected = sum(load(p)['retained_pairs'] for p in args.metrics)
        actual = bam_pairs(temporary)
        if actual != expected:
            raise ValueError('Merged Hi-C pair count differs from library totals')
        temporary.replace(out)
        libraries = [load(p) for p in args.metrics]
        modes = {r.get('hic_reference', 'combined') for r in libraries}
        if len(modes) != 1:
            raise ValueError('Hi-C libraries use inconsistent reference modes')
        mode = modes.pop()
        save(args.out, {'status': 'PASS', 'retained_pairs': actual, 'bam_sha256': checksum(out),
                        'hic_reference': mode, 'phase_validation': 'NOT_ASSESSED',
                        'count_unit': 'reference_pair_observations' if mode == 'per_haplotype' else 'physical_pairs',
                        'libraries': libraries})
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--metrics", nargs="+", required=True)
    merge_hic(parser.parse_args())


if __name__ == "__main__":
    main()
