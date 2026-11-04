#!/usr/bin/env python3
"""Verify checksums and every Hi-C FASTQ mate ID, count and length."""
import argparse
import hashlib
import itertools
from io_utils import load, save, source_check
from pathlib import Path


def paired_input(args):
    import pysam
    spec = load(args.spec)
    source_check(spec)
    hashes = []
    for i in (1, 2):
        path = Path(spec[f'file{i}'])
        md5, sha = hashlib.md5(), hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                md5.update(block)
                sha.update(block)
        expected = spec['expected_md5' if i == 1 else 'expected_md5_2']
        if expected and expected != md5.hexdigest():
            raise ValueError(f'MD5 mismatch: {path}')
        hashes.append({'file': str(path), 'md5': md5.hexdigest(), 'sha256': sha.hexdigest(),
                       'provided_md5_status': 'PASS' if expected else 'NOT_PROVIDED'})
    count, bases = 0, [0, 0]
    with pysam.FastxFile(spec['file1']) as f1, pysam.FastxFile(spec['file2']) as f2:
        for a, b in itertools.zip_longest(f1, f2):
            if a is None or b is None:
                raise ValueError('Paired FASTQ files have unequal record counts')
            names = [r.name[:-2] if r.name.endswith(('/1', '/2')) else r.name for r in (a, b)]
            if names[0] != names[1]:
                raise ValueError(f'Paired FASTQ IDs differ: {a.name}, {b.name}')
            for i, read in enumerate((a, b)):
                if not read.sequence or read.quality is None or len(read.sequence) != len(read.quality):
                    raise ValueError('Hi-C FASTQ sequence/quality is empty or malformed')
                bases[i] += len(read.sequence)
            count += 1
    if not count:
        raise ValueError('Empty paired Hi-C input')
    source_check(spec)
    save(args.out, {'status': 'PASS', 'sample_id': spec['sample_id'],
                    'readset_id': spec['readset_id'], 'library_id': spec['library_id'],
                    'num_pairs': count, 'num_reads': count * 2, 'num_bases': sum(bases),
                    'sources': hashes, 'source_spec': spec,
                    'validation': 'Full MD5/SHA256, full paired FASTQ parse, order/ID/length checks'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--out", required=True)
    paired_input(parser.parse_args())


if __name__ == "__main__":
    main()
