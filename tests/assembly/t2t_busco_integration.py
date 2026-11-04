#!/usr/bin/env python3
"""Validate native offline BUSCO execution and a deliberately failed gene gate."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

ROOT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    expected='08fb6d94a1699643c257653b77e05af5'
    if hashlib.md5(a.archive.read_bytes()).hexdigest()!=expected:raise ValueError('Official lineage archive MD5 mismatch')
    with tarfile.open(a.archive) as tar:tar.extractall(a.output/'database',filter='data')
    lineage=a.output/'database/eukaryota_odb12'
    if not (lineage/'dataset.cfg').exists():raise ValueError('Offline lineage missing')
    meta=json.loads((a.fixture/'run/config/samples/S.t2t_assessment.json').read_text())
    meta['settings']['assessment'].update(busco=True,busco_lineage=str(lineage),read_mapping=False,kmer=False)
    meta['tools']['busco']={'path':shutil.which('busco'),'sha256':hashlib.sha256(Path(shutil.which('busco')).read_bytes()).hexdigest()}
    meta['code']={rel:hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() for rel in meta['code']}
    metadata=a.output/'metadata.json';metadata.write_text(json.dumps(meta,indent=2)+'\n')
    command=[sys.executable,ROOT/'assembly/scripts/assess_t2t.py','--metadata',metadata,
             '--polishing',a.fixture/'run/samples/S/08_polishing/manifest.json',
             '--hifi',a.fixture/'run/samples/S/01_inputs/h.fastq.gz','--threads','8','--out-dir',a.output/'result']
    with (a.output/'execute.log').open('w') as log:subprocess.run(list(map(str,command)),stdout=log,stderr=log,check=True,timeout=3600)
    result=json.loads((a.output/'result/assessment.json').read_text())
    if result['evaluations']['busco']['status']!='COMPLETED' or result['busco_screen_pass'] is not False:
        raise ValueError('Random DNA must not pass the BUSCO completeness gate')
    if result['T2T_status']!='INCOMPLETE_QUALITY':raise ValueError('Failed BUSCO gate not propagated')
    (a.output/'validation.json').write_text(json.dumps({'status':'PASS','offline_busco':'PASS',
        'expected_gene_completeness_failure':'PASS','archive_md5':expected,'copies':2},indent=2)+'\n')


if __name__=='__main__':main()
