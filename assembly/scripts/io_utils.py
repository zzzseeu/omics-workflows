#!/usr/bin/env python3
"""Shared JSON, source integrity, process and atomic output helpers."""
import hashlib
import json
import shutil
import subprocess
import tempfile
from contextlib import ExitStack, contextmanager
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text())



def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + '\n')
    tmp.replace(path)



def signature(path):
    stat = Path(path).stat()
    return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
            'ctime_ns': stat.st_ctime_ns, 'inode': stat.st_ino}



def checksum(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()



def execute(command, **kwargs):
    print('COMMAND', json.dumps([str(x) for x in command]), flush=True)
    return subprocess.run([str(x) for x in command], check=True, **kwargs)



def source_check(spec):
    for i in (1, 2):
        key = 'source_signature' if i == 1 else 'source_signature_2'
        if spec.get(f'file{i}') and signature(spec[f'file{i}']) != spec[key]:
            raise ValueError(f'Source changed after preflight: {spec[f"file{i}"]}')



def publish(stage, destination, scratch):
    destination = Path(destination)
    if destination.exists():
        # Snakemake pre-creates every output parent, including nested empty
        # directories. Only actual files/symlinks require an ownership marker.
        if not (destination / '.omics-managed').is_file() and any(
                p.is_file() or p.is_symlink() for p in destination.rglob('*')):
            raise ValueError(f'Refusing to replace unmanaged directory: {destination}')
        destination.replace(Path(scratch) / 'previous')
    (stage / '.omics-managed').write_text('omics-workflows chromosome module\n')
    stage.replace(destination)



@contextmanager
def assembly_scratch(destination):
    temporary = Path(tempfile.mkdtemp(prefix='.joint-', dir=destination.parent))
    try:
        yield str(temporary)
    except BaseException:
        failed = destination.parent / 'failed_attempts'
        failed.mkdir(exist_ok=True)
        retained = failed / temporary.name.lstrip('.')
        temporary.replace(retained)
        print(f'FAILED_ASSEMBLY_RETAINED {retained.resolve()}', flush=True)
        raise
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)



def pipeline(commands, target):
    """No shell interpolation; fail when any producer fails, including SIGPIPE."""
    print('PIPELINE', json.dumps([[str(x) for x in c] for c in commands]), flush=True)
    processes = []
    try:
        with Path(target).open('wb') as out:
            previous = None
            for i, command in enumerate(commands):
                p = subprocess.Popen([str(x) for x in command], stdin=previous,
                                     stdout=out if i == len(commands)-1 else subprocess.PIPE)
                if previous is not None:
                    previous.close()
                previous = p.stdout
                processes.append(p)
            codes = [p.wait() for p in processes]
            if any(codes):
                raise RuntimeError(f'Pipeline subprocess failure: {codes}')
    finally:
        for p in processes:
            if p.poll() is None:
                p.terminate()
        for p in processes:
            p.wait()
