#!/usr/bin/env python3
"""Validate executables, dependencies and a small assembly outside the code repo."""
import argparse
from contextlib import ExitStack
import gzip
import importlib.metadata
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    checks = []

    def run(name, command, accepted=(0,), stdout_path=None):
        with ExitStack() as stack:
            log = stack.enter_context((out / (name + '.log')).open('w'))
            target = stack.enter_context(stdout_path.open('w')) if stdout_path else log
            result = subprocess.run(command, cwd=out, stdout=target, stderr=log, timeout=300)
        if result.returncode not in accepted:
            raise RuntimeError(f'{name} failed ({result.returncode}); see {out / (name + ".log")}')
        checks.append(name)

    executables = ['snakemake', 'hifiasm', 'samtools', 'seqkit', 'pigz', 'gfatools',
                   'minimap2', 'yak', 'fastp', 'busco', 'quast.py', 'meryl',
                   'merqury.sh', 'multiqc']
    paths = {}
    for name in executables:
        path = shutil.which(name)
        if not path or not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve()):
            raise RuntimeError(f'{name} missing from environment: {path}')
        paths[name] = path
    (out / 'executables.json').write_text(json.dumps(paths, indent=2) + '\n')
    import yaml
    import jsonschema
    plugin = importlib.metadata.version('snakemake-executor-plugin-cluster-generic')
    (out / 'cluster_plugin.version.txt').write_text(plugin + '\n')
    run('pip_check', [sys.executable, '-m', 'pip', 'check'])
    for name, command in {
        'snakemake_version': ['snakemake', '--version'],
        'hifiasm_version': ['hifiasm', '--version'],
        'samtools_version': ['samtools', '--version'],
        'seqkit_version': ['seqkit', 'version'],
        'pigz_version': ['pigz', '--version'],
        'minimap2_version': ['minimap2', '--version'],
        'fastp_version': ['fastp', '--version'],
        'busco_version': ['busco', '--version'],
        'quast_version': ['quast.py', '--version'],
        'meryl_version': ['meryl', '--version'],
        'multiqc_version': ['multiqc', '--version'],
    }.items():
        run(name, command)
    run('yak_help', ['yak'], accepted=(0, 1))
    run('gfatools_help', ['gfatools'], accepted=(0, 1))
    run('merqury_help', ['merqury.sh'], accepted=(0, 1))

    rng = random.Random(20261030)
    sequence = ''.join(rng.choices('ACGT', k=100000))
    (out / 'reference.fa').write_text('>synthetic_reference\n' + sequence + '\n')
    with gzip.open(out / 'reads.fq.gz', 'wt') as handle:
        for i in range(400):
            start = rng.randrange(0, 90001)
            read = sequence[start:start + 10000]
            handle.write(f'@synthetic_{i}\n{read}\n+\n' + 'I' * len(read) + '\n')
    run('hifiasm_smoke', ['hifiasm', '-t', '2', '-f0', '-l0', '-o', 'smoke', 'reads.fq.gz'])
    gfa = out / 'smoke.bp.p_ctg.gfa'
    if not gfa.is_file():
        raise RuntimeError('hifiasm did not produce primary GFA')
    run('gfa_to_fasta', ['gfatools', 'gfa2fa', str(gfa)], stdout_path=out / 'assembly.fa')
    if not (out / 'assembly.fa').read_text().startswith('>'):
        raise RuntimeError('No primary contigs in GFA')
    run('fasta_index', ['samtools', 'faidx', 'assembly.fa'])
    assembled_bases = sum(int(line.split('\t')[1]) for line in (out / 'assembly.fa.fai').read_text().splitlines())
    if assembled_bases < 80000:
        raise RuntimeError(f'Synthetic assembly is unexpectedly short: {assembled_bases}')
    run('assembly_stats', ['seqkit', 'stats', '-T', 'assembly.fa'])
    run('minimap2_smoke', ['minimap2', '-t', '2', '-x', 'map-hifi', 'assembly.fa', 'reads.fq.gz'])
    alignments = [line for line in (out / 'minimap2_smoke.log').read_text().splitlines()
                  if line.startswith('synthetic_') and len(line.split('\t')) >= 12]
    if len(alignments) < 300:
        raise RuntimeError(f'Too few synthetic reads aligned: {len(alignments)}')
    run('meryl_smoke', ['meryl', 'k=21', 'threads=2', 'count', 'output', 'reads.meryl', 'reads.fq.gz'])
    if not (out / 'reads.meryl' / 'merylIndex').is_file():
        raise RuntimeError('meryl did not produce a database index')
    run('quast_smoke', ['quast.py', 'assembly.fa', '-o', 'quast', '-t', '2'])
    if not (out / 'quast' / 'report.tsv').is_file():
        raise RuntimeError('QUAST report missing')
    (out / 'Snakefile').write_text('''rule all:
    input: "engine_stats.tsv"
rule stats:
    input: "assembly.fa"
    output: "engine_stats.tsv"
    shell: "seqkit stats -T {input:q} > {output:q}"
''')
    run('snakemake_dryrun', ['snakemake', '-n', '--cores', '1'])
    run('snakemake_execute', ['snakemake', '--cores', '1'])
    run('snakemake_restart', ['snakemake', '--cores', '1'])
    if 'Nothing to be done' not in (out / 'snakemake_restart.log').read_text():
        raise RuntimeError('Snakemake did not recognize completed outputs')
    run('cluster_plugin', ['snakemake', '--help'])
    if '--cluster-generic-submit-cmd' not in (out / 'cluster_plugin.log').read_text():
        raise RuntimeError('cluster-generic plugin is not registered')
    (out / 'validation.json').write_text(json.dumps({
        'status': 'PASS', 'prefix': sys.prefix, 'checks': checks,
        'synthetic_assembly_bases': assembled_bases, 'synthetic_alignments': len(alignments),
        'limitations': ['BUSCO lineage databases not downloaded',
                        'Merqury full evaluation not executed',
                        'CSUB executor adapter not yet implemented',
                        'Synthetic smoke assembly is not biological validation'],
    }, indent=2) + '\n')
    print('ASSEMBLY_ENVIRONMENT_VALIDATION=PASS', flush=True)


if __name__ == '__main__':
    main()
