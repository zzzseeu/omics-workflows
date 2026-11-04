#!/usr/bin/env python3
"""Install pinned official HapHiC source into an existing unified environment."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--env-prefix',type=Path,required=True)
    p.add_argument('--source-dir',type=Path,help='Existing official Git checkout; otherwise clone upstream')
    args=p.parse_args()
    prefix=args.env_prefix.resolve()
    if not (prefix/'conda-meta').is_dir():
        p.error('Supply the existing unified micromamba environment')
    lock=json.loads((ROOT/'env/haphic.lock.json').read_text())
    destination=prefix/'share'/('HapHiC-'+lock['commit'])
    executable=prefix/'bin/haphic'
    if executable.exists() and executable.resolve()!=destination/'haphic':
        p.error('Refusing to replace another HapHiC installation')
    if not destination.exists():
        with tempfile.TemporaryDirectory(prefix='omics-haphic-') as temp:
            source=args.source_dir
            if not source:
                source=Path(temp)/'source'
                subprocess.run(['git','clone',lock['project'],str(source)],check=True)
                subprocess.run(['git','checkout',lock['commit']],cwd=source,check=True)
            commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()
            if commit!=lock['commit']:
                p.error('HapHiC checkout does not match the source lock')
            if subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True).strip():
                p.error('HapHiC checkout contains unrecorded changes')
            shutil.copytree(source,destination,ignore=shutil.ignore_patterns('.git','__pycache__'))
    (destination/'SOURCE_LOCK.json').write_text(json.dumps(lock,indent=2)+'\n')
    if not executable.exists():executable.symlink_to(destination/'haphic')
    print(executable)


if __name__=='__main__':main()
