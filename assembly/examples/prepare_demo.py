#!/usr/bin/env python3
"""Create a small synthetic two-sample project outside the workflow repository."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]


def fastq(path, records):
    # Sequence/read sampling and gzip timestamps are deterministic.
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='ascii') as out:
                for name, sequence in records:
                    out.write(f'@{name}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', required=True, type=Path,
                        help='New project directory, outside the workflow repository')
    parser.add_argument('--seed', type=int, default=20261030)
    args = parser.parse_args()
    project = args.project_dir.resolve()
    if project == ROOT or project.is_relative_to(ROOT):
        parser.error('Synthetic data and results must be outside the code repository')
    if project.exists():
        parser.error('Project directory already exists; choose a new one (nothing overwritten)')
    samtools = Path(sys.prefix) / 'bin/samtools'
    if not samtools.is_file():
        parser.error('Run with the unified environment Python; samtools is required')
    project.mkdir(parents=True)
    data, config_dir = project / 'data', project / 'config'
    data.mkdir()
    shutil.copytree(ROOT / 'assembly/examples/hifi_only', config_dir)
    rng = random.Random(args.seed)
    for sid in ['sample_A', 'sample_B']:
        reference = ''.join(rng.choices('ACGT', k=100000))
        (data / f'{sid}.synthetic_reference.fa').write_text(f'>{sid}_synthetic\n{reference}\n')
        records = []
        for i in range(400):
            start = rng.randrange(0, 90001)
            records.append((f'{sid}_read_{i}', reference[start:start + 10000]))
        if sid == 'sample_A':
            sam = data / 'sample_A.cell01.sam'
            with sam.open('w') as out:
                out.write('@HD\tVN:1.6\tSO:unknown\n')
                for name, sequence in records[:200]:
                    out.write(f'{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t'
                              + 'I' * len(sequence) + '\n')
            subprocess.run([str(samtools), 'view', '-b', '-o', str(data / 'sample_A.cell01.bam'),
                            str(sam)], check=True)
            fastq(data / 'sample_A.cell02.fastq.gz', records[200:])
        else:
            fastq(data / 'sample_B.cell01.fastq.gz', records)
    table = config_dir / 'readsets.tsv'
    with table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        row['md5_1'] = hashlib.md5((config_dir / row['file1']).read_bytes()).hexdigest()
    with table.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=[*rows[0].keys()], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    config_path = config_dir / 'config.yaml'
    config = yaml.safe_load(config_path.read_text())
    config['analysis_profiles']['hifi_contigs'] = {'assembly': {
        'threads': 2, 'mem_mb': 2048, 'bloom_filter_bits': 0, 'purge_level': 0}}
    for sid in ['sample_A', 'sample_B']:
        config['samples'][sid]['ploidy'] = 1
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    (config_dir / 'samples.tsv').write_text(
        'sample_id\tspecies\tprofile\n'
        'sample_A\tSynthetic sequence A\thifi_contigs\n'
        'sample_B\tSynthetic sequence B\thifi_contigs\n')
    (project / 'DEMO_ONLY.json').write_text(json.dumps({
        'purpose': 'Execution example only; not biological or T2T validation',
        'seed': args.seed, 'samples': ['sample_A', 'sample_B'],
        'reference_bases_per_sample': 100000, 'reads_per_sample': 400,
        'bases_per_read': 10000, 'profile': 'small synthetic data only',
        'note': 'BAM bytes may differ across samtools versions and project paths',
    }, indent=2) + '\n')
    print(f'SYNTHETIC_DEMO_PROJECT={project}')
    print('Config: ' + str(config_path))
    print('Synthetic-only profile: 2 threads / 2048 MB per assembly; NOT for real genomes')


if __name__ == '__main__':
    main()
