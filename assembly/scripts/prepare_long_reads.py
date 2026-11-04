#!/usr/bin/env python3
"""Verify long-read sources and prepare lossless compressed FASTQ."""
import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
from io_utils import execute, load, save, signature
from pathlib import Path


def prepare(args):
    spec = load(args.spec)
    source = Path(spec['file1'])
    before = signature(source)
    if before != spec['source_signature']:
        raise ValueError('Source changed after preflight; restart the launcher to create a new plan')
    md5, sha = hashlib.md5(), hashlib.sha256()
    with source.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            md5.update(block)
            sha.update(block)
    if spec['expected_md5'] and md5.hexdigest() != spec['expected_md5']:
        raise ValueError(f'MD5 mismatch: {source}; observed {md5.hexdigest()}')
    output = Path(args.fastq)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.partial.fastq.gz')
    try:
        if spec['format'] == 'bam':
            execute(['samtools', 'quickcheck', '-u', '-v', source])
            paired = execute(['samtools', 'view', '-c', '-f', '1', source],
                             stdout=subprocess.PIPE).stdout.decode().strip()
            if int(paired):
                raise ValueError('HiFi BAM contains paired reads; only unpaired HiFi BAM is supported')
            # Unpaired primary reads go to category 0. Native compression also works
            # with a one-core budget without a second Python compressor process.
            execute(['samtools', 'fastq', '-n', '-@', str(max(0, args.threads - 1)),
                     '-0', temporary, '-s', os.devnull, source], stdout=subprocess.DEVNULL)
            conversion = 'samtools fastq; unpaired primary reads; original IDs retained'
        else:
            with source.open('rb') as probe:
                compressed = probe.read(2) == b'\x1f\x8b'
            if compressed:
                shutil.copyfile(source, temporary)
            else:
                with source.open('rb') as src, gzip.open(temporary, 'wb', compresslevel=3) as target:
                    shutil.copyfileobj(src, target, 8 * 1024 * 1024)
            conversion = 'FASTQ byte copy (gzip) or lossless gzip compression (plain text)'
        if signature(source) != before:
            raise ValueError('Source changed while preparing reads')
        # Compressed FASTQ byte copies already have exactly the source hash.
        # Converted/compressed outputs get a producer hash for future reuse.
        if spec['format']=='fastq' and compressed:
            prepared_sha=sha.hexdigest()
        else:
            prepared=hashlib.sha256()
            with temporary.open('rb') as handle:
                for block in iter(lambda:handle.read(8*1024*1024),b''):prepared.update(block)
            prepared_sha=prepared.hexdigest()
        temporary.replace(output)
        save(args.manifest, {'sample_id': spec['sample_id'], 'readset_id': spec['readset_id'],
                             'library_id': spec['library_id'], 'source': str(source),
                             'source_signature': before, 'source_md5': md5.hexdigest(),
                             'source_sha256': sha.hexdigest(),
                             'provided_md5_status': 'PASS' if spec['expected_md5'] else 'NOT_PROVIDED',
                             'preparation_format_version': 2, 'technology':spec['technology'],
                             'source_format':spec['format'],
                             'prepared_fastq': str(output.absolute()), 'conversion': conversion,
                             'prepared_sha256':prepared_sha,'prepared_signature':signature(output),
                             'filtering': 'No quality trimming or deduplication performed'})
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--threads", required=True, type=int)
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
