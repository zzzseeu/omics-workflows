#!/usr/bin/env python3
"""Native zero-copy prepared-input, shared-cache, QC and recovery validation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import yaml

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'assembly/scripts'))
from io_utils import checksum,save,signature


def call(command,log,success=True):
    with log.open('w') as out:r=subprocess.run(list(map(str,command)),stdout=out,stderr=out,timeout=1200)
    if (r.returncode==0)!=success:raise RuntimeError(f'Unexpected exit {r.returncode}: {log}')
    return log.read_text()


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);old=out/'old/samples/A'
    (old/'01_inputs').mkdir(parents=True);(old/'02_read_qc').mkdir()
    rng=random.Random(20261104);sequence=''.join(rng.choices('ACGT',k=100000));rows=[];protected={};modes={}
    for rid,mode in [('legacy','legacy'),('producer','producer'),('touched','touched')]:
        raw=out/(rid+('.bam' if mode=='legacy' else '.fq.gz'))
        reads=[]
        for i in range(200):
            start=rng.randrange(90001);reads.append((rid+'_'+str(i),sequence[start:start+10000]))
        if mode=='legacy':
            sam=out/'legacy.sam'
            with sam.open('w') as f:
                f.write('@HD\tVN:1.6\tSO:unknown\n')
                for name,seq in reads:f.write(f'{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t'+('I'*len(seq))+'\n')
            call(['samtools','view','-b','-o',raw,sam],out/'make_bam.log');fmt='bam'
        else:
            with gzip.open(raw,'wt') as f:
                for name,seq in reads:f.write('@'+name+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
            fmt='fastq'
        spec={'sample_id':'A','readset_id':rid,'library_id':rid,'technology':'hifi','format':fmt,
              'file1':str(raw),'source_signature':signature(raw),'expected_md5':checksum(raw,'md5')}
        source_spec=out/(rid+'.spec.json');save(source_spec,spec)
        fq=old/'01_inputs'/(rid+'.fastq.gz');manifest=old/'01_inputs'/(rid+'.json')
        call([sys.executable,ROOT/'assembly/scripts/prepare_long_reads.py','--spec',source_spec,
              '--fastq',fq,'--manifest',manifest,'--threads','2'],out/(rid+'.prepare.log'))
        qc=old/'02_read_qc'/(rid+'.json');stats=old/'02_read_qc'/(rid+'.tsv')
        call([sys.executable,ROOT/'assembly/scripts/check_read_quality.py','--fastq',fq,
              '--prepared-manifest',manifest,'--stats',stats,'--json',qc,'--seqkit','seqkit','--threads','2'],out/(rid+'.qc.log'))
        if mode=='legacy':
            data=json.loads(manifest.read_text())
            for k in ['prepared_sha256','prepared_signature','preparation_format_version']:data.pop(k,None)
            data['prepared_fastq']=f'samples/A/01_inputs/{rid}.fastq.gz';save(manifest,data)
        if mode=='touched':
            stat=fq.stat();os.utime(fq,ns=(stat.st_atime_ns,stat.st_mtime_ns-60_000_000_000))
        protected[rid]={'sha256':checksum(fq),'signature':signature(fq)};modes[rid]=str(fq)
        rows.append(['A',rid,rid,'hifi',fmt,str(raw),'',spec['expected_md5'],'',str(manifest)])
    (out/'samples.tsv').write_text('sample_id\tspecies\tprofile\nA\tSynthetic\t\n')
    with (out/'readsets.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample_id','readset_id','library_id',
            'technology','format','file1','file2','md5_1','md5_2','prepared_manifest']);w.writerows(rows)
    config={'schema_version':1,'samples_table':'samples.tsv','readsets_table':'readsets.tsv',
            'prepared_input_cache':str(out/'cache'),
            'defaults':{'assembly':{'threads':4,'mem_mb':4096,'bloom_filter_bits':0}}}
    cfg=out/'config.yaml';cfg.write_text(yaml.safe_dump(config))
    def command(run):return [sys.executable,ROOT/'run','assembly','--config',cfg,'--run-dir',out/run,
                             '--env-prefix',sys.prefix,'--cores','4','--mem-mb','8192','--target','contigs']
    call([*command('run1'),'--dry-run'],out/'dryrun.log')
    # Two consumers share one attestation per source, protected by flock.
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(call,command(name),out/(name+'.log')) for name in ['run1','run2']]
        for f in futures:f.result()
    for rid,path in modes.items():
        f=Path(path)
        if signature(f)!=protected[rid]['signature'] or checksum(f)!=protected[rid]['sha256']:
            raise ValueError('Original prepared FASTQ was touched or changed')
        manifests=[]
        for name in ['run1','run2']:
            output=out/name/'samples/A/01_inputs'/(rid+'.fastq.gz')
            if not output.is_symlink() or output.resolve()!=f:raise ValueError('Prepared FASTQ was copied')
            m=json.loads(output.with_name(rid+'.json').read_text());manifests.append(m)
            qc=json.loads((out/name/'samples/A/02_read_qc'/(rid+'.json')).read_text())
            if 'reuse' not in qc:raise ValueError('Read QC was needlessly recomputed')
        if sum(not m['reuse']['cache_hit'] for m in manifests)!=1:raise ValueError('Shared cache scanned a source twice')
    methods=[json.loads(f.read_text())['validation'] for f in (out/'cache').glob('*.json')]
    for expected in ['LEGACY_ONE_PASS_COMPRESSED_SHA256_AND_FASTQ_QC_COUNTS',
                     'PRODUCER_SHA256_SIGNATURE_AND_PASS_QC','PRODUCER_SHA256_REVERIFIED_AFTER_TIMESTAMP_CHANGE']:
        if methods.count(expected)!=1:raise ValueError('Wrong attestation mode: '+str(methods))
    if 'Nothing to be done' not in call(command('run1'),out/'resume.log'):raise ValueError('Resume recomputed inputs')
    # A third independent consumer reuses every cache entry without scanning.
    call(command('run3'),out/'run3.log')
    for log in (out/'run3/samples/A/logs/01_inputs').glob('*.log'):
        if 'VERIFY_LEGACY_PREPARED_FASTQ' in log.read_text():raise ValueError('Cached source was rescanned')
    # Signature drift after certification must fail, never become a new baseline.
    Path(modes['legacy']).touch()
    call(command('run4'),out/'changed_input.log',success=False)
    errors=''.join(f.read_text() for f in (out/'run4/samples/A/logs/01_inputs').glob('*.log'))
    if 'Cached prepared FASTQ changed' not in errors:raise ValueError('Changed-input failure was not provenance-bound')
    save(out/'validation.json',{'status':'PASS','native_bam_and_fastq':'PASS','symlink_no_copy':'PASS',
         'original_fastq_signature_unchanged_during_reuse':'PASS','shared_cache_single_attestation':'PASS',
         'three_attestation_modes':methods,'QC_reuse':'PASS','resume':'PASS','changed_input_blocks':'PASS'})


if __name__=='__main__':main()
