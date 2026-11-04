#!/usr/bin/env python3
"""End-to-end portable DAG validation. Run once inside an outer compute job."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
import yaml

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    if out == ROOT or out.is_relative_to(ROOT):
        raise ValueError('Test output must be outside the code repository')
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    process_env = dict(os.environ, PATH=str(Path(sys.prefix) / 'bin') + os.pathsep + os.defpath)
    if shutil.which('csub', path=process_env['PATH']):
        raise AssertionError('Portable validation PATH unexpectedly contains csub')
    checks.append('scheduler_free_PATH')

    def run(name, command, success=True):
        with (out / f'{name}.log').open('w') as log:
            completed = subprocess.run([str(c) for c in command], stdout=log,
                                       stderr=subprocess.STDOUT, timeout=900, env=process_env)
        if (completed.returncode == 0) != success:
            raise RuntimeError(f'{name}: unexpected exit {completed.returncode}; inspect {out / (name + ".log")}')
        checks.append(name)
        print('PASS', name, flush=True)

    run('config_contracts', [sys.executable, '-B', '-m', 'unittest', 'discover',
                             '-s', ROOT / 'tests/assembly', '-p', 'test_*.py', '-v'])
    fixture = out / 'fixture with spaces'
    fixture.mkdir()
    rng = random.Random(20261030)
    rows = []
    for sid in ['A', 'B']:
        sequence = ''.join(rng.choices('ACGT', k=100000))
        records = []
        for i in range(400):
            start = rng.randrange(0, 90001)
            records.append((f'{sid}_read_{i}', sequence[start:start + 10000]))
        groups = [records[:200], records[200:]] if sid == 'A' else [records]
        for i, group in enumerate(groups, 1):
            if sid == 'A' and i == 1:
                sam = fixture / 'A cell1.sam'
                with sam.open('w') as f:
                    f.write('@HD\tVN:1.6\tSO:unknown\n')
                    for name, read in group:
                        f.write(f'{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{read}\t' + 'I' * len(read) + '\n')
                path = fixture / 'A cell1.bam'
                run('make_unpaired_bam', ['samtools', 'view', '-b', '-o', path, sam])
                fmt = 'bam'
            else:
                path = fixture / f'{sid} cell{i}.fq'  # Plain FASTQ conversion path.
                with path.open('w') as f:
                    for name, read in group:
                        f.write(f'@{name}\n{read}\n+\n' + 'I' * len(read) + '\n')
                if sid == 'B':
                    zipped = Path(str(path) + '.gz')
                    with path.open('rb') as src, gzip.open(zipped, 'wb') as dst:
                        dst.write(src.read())
                    path = zipped
                fmt = 'fastq'
            rows.append([sid, f'cell{i}', 'library1', 'hifi', fmt, path.name, '',
                         hashlib.md5(path.read_bytes()).hexdigest(), ''])
    (fixture / 'samples.tsv').write_text('sample_id\tspecies\tprofile\nA\tSynthetic A\tsmall\nB\tSynthetic B\tsmall\n')
    with (fixture / 'readsets.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'readset_id', 'library_id', 'technology', 'format',
                         'file1', 'file2', 'md5_1', 'md5_2'])
        writer.writerows(rows)
    cfg = {'schema_version': 1, 'samples_table': 'samples.tsv', 'readsets_table': 'readsets.tsv',
           'analysis_profiles': {'small': {'assembly': {'threads': 2, 'mem_mb': 2048,
                                                       'bloom_filter_bits': 0, 'purge_level': 0}}}}
    config_path = fixture / 'config.yaml'
    config_path.write_text(yaml.safe_dump(cfg))
    run_dir = out / 'portable run'
    base = [sys.executable, ROOT / 'run', 'assembly', '--config', config_path,
            '--run-dir', run_dir, '--cores', '4', '--mem-mb', '4096']
    if Path('/usr/bin/python3').exists():
        run('environment_bootstrap', ['/usr/bin/python3', ROOT / 'run', 'assembly',
                                      '--env-prefix', sys.prefix, '--config', config_path,
                                      '--run-dir', out / 'check only run', '--samples', 'A',
                                      '--cores', '2', '--mem-mb', '4096', '--check'])
        if (out / 'check only run').exists():
            raise AssertionError('--check created a run directory')
    run('dryrun', [*base, '--samples', 'A', '--dry-run'])
    if (run_dir / 'samples/A/04_assembly/hifiasm').exists():
        raise AssertionError('Dry-run executed the assembler')
    run('selected_A', [*base, '--samples', 'A'])
    a = run_dir / 'samples/A/04_assembly/hifiasm/primary.fa'
    if (run_dir / 'samples/B').exists():
        raise AssertionError('Unselected B ran')
    a_time = a.stat().st_mtime_ns
    run('add_sample_B', base)
    if a.stat().st_mtime_ns != a_time:
        raise AssertionError('Adding B unnecessarily recomputed A assembly')
    for sid in ['A', 'B']:
        report = json.loads((run_dir / f'samples/{sid}/10_report/manifest.json').read_text())
        stat = report['assessment']['sets']['primary']
        if stat['total_bases'] < 80000 or report['assessment']['input_reads'] != 400:
            raise AssertionError(f'{sid}: incorrect synthetic output/input count')
        if report['T2T_status'] != 'NOT_ASSESSED':
            raise AssertionError('Synthetic result incorrectly certified T2T')
    before = {sid: (run_dir / f'samples/{sid}/04_assembly/hifiasm/primary.fa').stat().st_mtime_ns
              for sid in ['A', 'B']}
    run('resume_no_compute', base)
    if 'Nothing to be done' not in (out / 'resume_no_compute.log').read_text():
        raise AssertionError('Resume did not recognize completed DAG')
    cfg['samples'] = {'B': {'assembly': {'threads': 1}}}
    config_path.write_text(yaml.safe_dump(cfg))
    run('change_B_parameters', base)
    after = {sid: (run_dir / f'samples/{sid}/04_assembly/hifiasm/primary.fa').stat().st_mtime_ns
             for sid in ['A', 'B']}
    if after['A'] != before['A'] or after['B'] == before['B']:
        raise AssertionError('Sample parameter invalidation/isolation failed')
    # Remove one final product and verify dependency-aware recovery.
    (run_dir / 'samples/B/10_report/report.html').unlink()
    run('recover_missing_report', base)
    if not (run_dir / 'samples/B/10_report/report.html').is_file():
        raise AssertionError('Missing HTML report was not recovered')
    if (run_dir / 'samples/B/04_assembly/hifiasm/primary.fa').stat().st_mtime_ns != after['B']:
        raise AssertionError('Recovering a report unnecessarily reran assembly')
    (run_dir / 'samples/B/04_assembly/hifiasm/primary.fa.fai').unlink()
    run('recover_missing_index', base)
    if (run_dir / 'samples/B/04_assembly/hifiasm/primary.fa').stat().st_mtime_ns != after['B']:
        raise AssertionError('Recovering an index unnecessarily reran assembly')
    (run_dir / 'samples/B/04_assembly/hifiasm/primary.fa').unlink()
    run('recover_missing_fasta', base)
    if not (run_dir / 'samples/B/04_assembly/hifiasm/primary.fa').is_file():
        raise AssertionError('Missing FASTA was not recovered')
    if a.stat().st_mtime_ns != after['A']:
        raise AssertionError('Recovery of B unexpectedly recomputed A')
    # A valid file with an incorrect delivered checksum must stop before assembly.
    wrong = [r[:] for r in rows]
    wrong[0][7] = '0' * 32
    with (fixture / 'bad_md5.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'readset_id', 'library_id', 'technology', 'format',
                         'file1', 'file2', 'md5_1', 'md5_2'])
        writer.writerows(wrong)
    bad = {**cfg, 'readsets_table': 'bad_md5.tsv'}
    bad_path = fixture / 'bad.yaml'
    bad_path.write_text(yaml.safe_dump(bad))
    run('checksum_failure', [sys.executable, ROOT / 'run', 'assembly', '--config', bad_path,
                             '--run-dir', out / 'bad checksum run', '--samples', 'A',
                             '--cores', '2', '--mem-mb', '4096'], success=False)
    if 'MD5 mismatch' not in (out / 'bad checksum run/samples/A/logs/01_inputs/cell1.log').read_text():
        raise AssertionError('Wrong failure reason in checksum test')
    # Terminate the launcher during a real task; it must forward the signal and release
    # the run lock only after its local Snakemake process exits, then resume successfully.
    interrupted_run = out / 'interrupted run'
    interrupted_command = [sys.executable, ROOT / 'run', 'assembly', '--config', config_path,
                           '--run-dir', interrupted_run, '--samples', 'A',
                           '--cores', '2', '--mem-mb', '4096']
    with (out / 'terminate.log').open('w') as log:
        process = subprocess.Popen([str(x) for x in interrupted_command], env=process_env,
                                   stdout=log, stderr=subprocess.STDOUT)
        deadline = time.monotonic() + 90
        while not (interrupted_run / 'samples/A/01_inputs/cell1.fastq.gz').exists():
            if process.poll() is not None or time.monotonic() > deadline:
                process.kill()
                process.wait()
                raise AssertionError('Could not reach a running input task for interruption test')
            time.sleep(0.1)
        process.terminate()
        if process.wait(timeout=60) == 0:
            raise AssertionError('Terminated launcher returned success')
    state = json.loads((interrupted_run / 'run_status.json').read_text())
    if state['status'] != 'FAILED':
        raise AssertionError('Terminated run status was not recorded')
    checks.append('terminate_and_record_status')
    # A forced worker shutdown can leave a Snakemake lock. Exercise the documented
    # recovery procedure only after verifying that this test run has no live workers.
    if list((interrupted_run / '.snakemake/locks').glob('*.lock')):
        import psutil
        for worker in psutil.process_iter():
            try:
                if worker.status() != psutil.STATUS_ZOMBIE and worker.cwd() == str(interrupted_run):
                    raise AssertionError(f'Interrupted workflow still has a live worker: {worker.pid}')
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        run('unlock_after_forced_shutdown', [Path(sys.prefix) / 'bin/snakemake',
            '--snakefile', ROOT / 'assembly/Snakefile', '--directory', interrupted_run,
            '--configfile', interrupted_run / 'config/plan.json', '--unlock'])
    run('resume_interrupted', interrupted_command)
    # Invalid FASTQ must not reach the assembler, even without a supplied checksum.
    invalid = fixture / 'invalid.fq'
    invalid.write_text('@bad\nACGTACGT\n+\nII\n')
    (fixture / 'invalid.tsv').write_text(
        'sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\n'
        'A\tbad\tlib\thifi\tfastq\tinvalid.fq\n')
    bad_path.write_text(yaml.safe_dump({**cfg, 'readsets_table': 'invalid.tsv'}))
    run('invalid_fastq', [sys.executable, ROOT / 'run', 'assembly', '--config', bad_path,
                          '--run-dir', out / 'invalid fastq run', '--samples', 'A',
                          '--cores', '2', '--mem-mb', '4096'], success=False)
    if (out / 'invalid fastq run/samples/A/04_assembly/hifiasm').exists():
        raise AssertionError('Malformed FASTQ reached assembly')
    result = {'status': 'PASS', 'checks': checks, 'samples': ['A', 'B'],
              'limitations': ['Synthetic execution validation only; not biological or T2T validation'],
              'run_directory': str(run_dir)}
    (out / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print('ASSEMBLY_MVP_VALIDATION=PASS', flush=True)


if __name__ == '__main__':
    main()
