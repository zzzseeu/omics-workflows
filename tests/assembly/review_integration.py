#!/usr/bin/env python3
"""Verify edited-AGP application and isolation of review changes from assembly."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import yaml

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'assembly/scripts'))
from io_utils import load,save,checksum


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,help='New review-test directory; default: fixture/review_fixture')
    args=p.parse_args()
    validation=args.validation_dir.resolve()
    fixture=validation/'known_structure'
    review=(args.output or fixture/'review_fixture').resolve()
    if review == ROOT or review.is_relative_to(ROOT):
        raise ValueError('Review validation output must be outside the code repository')
    review.mkdir(exist_ok=False)
    env=dict(os.environ,PATH=str(Path(sys.prefix)/'bin')+os.pathsep+os.defpath,
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    native=fixture/'haphic/04.build/scaffolds.raw.agp'
    rows=native.read_text().splitlines()
    for i,line in enumerate(rows):
        c=line.split('\t')
        if c[4]=='W':
            c[8]='-' if c[8]=='+' else '+'
            rows[i]='\t'.join(c)
            break
    edited=review/'edited.raw.agp'
    edited.write_text('\n'.join(rows)+'\n')
    request=load(fixture/'delivery/review.request.json')
    request.update(status='APPROVED',reviewer='SYNTHETIC_TEST_ONLY',reviewed_at='2026-11-03',
                   review_agp_sha256=checksum(edited))
    save(review/'review.json',request)
    metadata=load(fixture/'metadata.json')
    metadata['settings']['scaffolding'].update(review_agp=str(edited),review_manifest=str(review/'review.json'))
    save(review/'metadata.json',metadata)
    with (review/'apply.log').open('w') as log:
        subprocess.run([sys.executable,str(ROOT/'assembly/scripts/export_chromosome_fasta.py'),
            '--metadata',str(review/'metadata.json'),'--assembly',str(fixture/'assembly.json'),
            '--fasta',str(fixture/'combined.fa'),'--mapping',str(fixture/'mapping.tsv'),
            '--haphic',str(fixture/'haphic'),'--out-dir',str(review/'delivery')],
            env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    applied=load(review/'delivery/manifest.json')
    assert applied['human_review']=='APPROVED' and applied['quality_level']=='CHROMOSOME_REVIEWED'
    assert checksum(review/'delivery/chromosomes.fa')!=checksum(fixture/'delivery/chromosomes.fa')
    assert applied['source_coverage_validation']=='PASS_EXACTLY_ONCE_ALL_BASES'
    # Planning a review of the original joint run must leave its expensive stage specs unchanged.
    run=validation/'joint_run'
    data=validation/'joint fixture'
    # Refresh the normal plan first, so a code update between validation jobs
    # cannot be mistaken for invalidation caused by review parameters.
    with (review/'baseline_dag.log').open('w') as log:
        subprocess.run([str(ROOT/'run'), 'assembly','--env-prefix',sys.prefix,'--config',str(data/'config.yaml'),
            '--run-dir',str(run),'--cores','8','--mem-mb','64000','--target','chromosomes','--dry-run'],
            env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    specs=[run/f'config/samples/S.{stage}.json' for stage in ('joint','linkage','scaffolding')]
    before={str(path):checksum(path) for path in specs}
    config=yaml.safe_load((data/'config.yaml').read_text())
    config['defaults']['scaffolding'].update(review_agp=str(edited),review_manifest=str(review/'review.json'))
    (data/'review.config.yaml').write_text(yaml.safe_dump(config))
    with (review/'review_dag.log').open('w') as log:
        subprocess.run([str(ROOT/'run'), 'assembly','--env-prefix',sys.prefix,'--config',str(data/'review.config.yaml'),
            '--run-dir',str(run),'--cores','8','--mem-mb','64000','--target','chromosomes','--dry-run'],
            env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    assert before=={str(path):checksum(path) for path in specs}, 'Review changed assembly/linkage/scaffold specs'
    assert 'rule joint_hifiasm:' not in (review/'review_dag.log').read_text(), 'Review planned reassembly'
    save(review/'validation.json',{'status':'PASS','checks':['edited_AGP_applied','hashes_bound_to_review',
        'FASTA_changes_with_orientation','source_coverage_preserved','review_does_not_reassemble'],
        'scope':'Synthetic interface test; no real biological result approved'})
    print('REVIEW_INTEGRATION=PASS',flush=True)


if __name__=='__main__':main()
