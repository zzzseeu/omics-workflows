#!/usr/bin/env python3
"""Attest once, then link immutable prepared FASTQ without reconversion or copying."""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
from pathlib import Path
import re
import subprocess
import tempfile
from io_utils import checksum, load, save, signature, source_check
from prepared_input_utils import POLICY


def unchanged(spec):
    source_check(spec);r=spec['prepared_reuse'];identity=r['identity']
    for key, sigkey, hashkey in [('manifest','manifest_signature','manifest_sha256'),('qc','qc_signature','qc_sha256')]:
        path=identity[key]
        if signature(path)!=r[sigkey] or checksum(path)!=identity[hashkey]:
            raise ValueError('Prepared-input metadata changed after preflight: '+key)
    if signature(identity['prepared_fastq'])!=r['prepared_signature']:
        raise ValueError('Prepared FASTQ changed after preflight')


def scan_legacy(fastq, seqkit, expected_reads, expected_bases, folder):
    """Hash compressed bytes and stream once through native seqkit FASTQ parsing."""
    sha=hashlib.sha256()
    with tempfile.TemporaryFile(mode='w+b',dir=folder) as stats:
        command=[seqkit,'stats','-T','-j','1','-']
        print('VERIFY_LEGACY_PREPARED_FASTQ',fastq,flush=True)
        process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=stats)
        try:
            with Path(fastq).open('rb') as src:
                for block in iter(lambda:src.read(8*1024*1024),b''):
                    sha.update(block);process.stdin.write(block)
            process.stdin.close()
            if process.wait():raise ValueError('Prepared FASTQ parsing failed')
        finally:
            if process.stdin and not process.stdin.closed:process.stdin.close()
            if process.poll() is None:process.terminate();process.wait()
        stats.seek(0);rows=list(csv.DictReader(stats.read().decode().splitlines(),delimiter='\t'))
        if len(rows)!=1 or rows[0]['format'].upper()!='FASTQ':raise ValueError('Prepared input is not FASTQ')
        reads=int(rows[0]['num_seqs'].replace(',',''));bases=int(rows[0]['sum_len'].replace(',',''))
        if (reads,bases)!=(expected_reads,expected_bases):raise ValueError('Prepared FASTQ does not match prior QC counts')
    return sha.hexdigest()


def attest(spec):
    unchanged(spec);r=spec['prepared_reuse'];identity=r['identity'];cache=Path(r['cache_file'])
    cache.parent.mkdir(parents=True,exist_ok=True)
    with cache.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        unchanged(spec)
        hit=cache.exists()
        if hit:
            if cache.is_symlink():raise ValueError('Prepared-input cache record must not be a symlink')
            record=load(cache)
            if record.get('identity')!=identity or record.get('status')!='PASS' or record.get('prepared_signature')!=r['prepared_signature']:
                raise ValueError('Cached prepared FASTQ changed or cache provenance is invalid')
            if not re.fullmatch('[0-9a-f]{64}',record.get('prepared_sha256','')):
                raise ValueError('Prepared-input cache checksum is invalid')
            if r['producer_sha256'] and record['prepared_sha256']!=r['producer_sha256']:
                raise ValueError('Prepared-input cache checksum disagrees with producer')
        else:
            if r['producer_sha256'] and r['producer_signature_matches']:
                sha=r['producer_sha256'];method='PRODUCER_SHA256_SIGNATURE_AND_PASS_QC'
            else:
                sha=scan_legacy(identity['prepared_fastq'],spec['seqkit']['path'],
                                r['expected_num_reads'],r['expected_num_bases'],cache.parent)
                if r['producer_sha256'] and sha!=r['producer_sha256']:
                    raise ValueError('Prepared FASTQ producer SHA256 mismatch')
                method=('PRODUCER_SHA256_REVERIFIED_AFTER_TIMESTAMP_CHANGE' if r['producer_sha256'] else
                        'LEGACY_ONE_PASS_COMPRESSED_SHA256_AND_FASTQ_QC_COUNTS')
            unchanged(spec)
            record={'status':'PASS','identity':identity,'prepared_signature':r['prepared_signature'],
                    'prepared_sha256':sha,'validation':method,
                    'num_reads':r['expected_num_reads'],'num_bases':r['expected_num_bases'],
                    'created_utc':datetime.now(timezone.utc).isoformat()}
            save(cache,record)
    unchanged(spec)
    return record,hit,checksum(cache)


def reuse(args):
    spec=load(args.spec);record,hit,cache_sha=attest(spec);r=spec['prepared_reuse'];source=Path(r['identity']['prepared_fastq'])
    output=Path(args.fastq);output.parent.mkdir(parents=True,exist_ok=True)
    if output.resolve()==source:raise ValueError('Reuse destination aliases the original FASTQ')
    old=load(r['identity']['manifest'])
    temporary=output.with_name(output.name+'.reuse.partial')
    try:
        if temporary.exists() or temporary.is_symlink():temporary.unlink()
        temporary.symlink_to(source)
        unchanged(spec);temporary.replace(output)
        # Preserve the original raw-source provenance; never replace the raw
        # source hash with the compressed FASTQ hash.
        save(args.manifest,{**{k:old[k] for k in ('sample_id','readset_id','library_id','source',
            'source_signature','source_md5','source_sha256','provided_md5_status','conversion','filtering')},
            'preparation_format_version':2,'technology':spec['technology'],'source_format':spec['format'],
            'prepared_fastq':str(output.absolute()),'prepared_sha256':record['prepared_sha256'],
            'prepared_signature':record['prepared_signature'],
            'reuse':{'policy':POLICY,'original_manifest':r['identity']['manifest'],
                     'original_manifest_sha256':r['identity']['manifest_sha256'],
                     'qc_manifest':r['identity']['qc'],'qc_manifest_sha256':r['identity']['qc_sha256'],
                     'prepared_target':str(source),'cache_file':r['cache_file'],'cache_sha256':cache_sha,
                     'validation':record['validation'],'cache_hit':hit,'storage':'SYMLINK_NO_COPY'}})
    finally:
        if temporary.exists() or temporary.is_symlink():temporary.unlink()
    print('PREPARED_INPUT_REUSED',spec['sample_id'],spec['readset_id'],'CACHE_HIT' if hit else 'ATTESTED',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('spec','fastq','manifest'):p.add_argument('--'+key,required=True)
    p.add_argument('--threads',type=int,required=True)
    reuse(p.parse_args())


if __name__=='__main__':main()
