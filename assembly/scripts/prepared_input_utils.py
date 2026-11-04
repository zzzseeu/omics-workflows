"""Metadata-only identity checks for reusing numbered long-read preparation outputs."""
import hashlib
import json
from pathlib import Path
import re

POLICY = 'lossless-prepared-fastq-reuse-v1'


def signature(path):
    stat = Path(path).stat()
    return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
            'ctime_ns': stat.st_ctime_ns, 'inode': stat.st_ino}


def metadata(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > 1024*1024:
        raise ValueError(f'Missing or oversized prepared-input metadata: {path}')
    content = path.read_bytes()
    return json.loads(content), hashlib.sha256(content).hexdigest()


def describe_prepared_input(spec, manifest, cache_dir, run_dir, repository):
    """No reads, decompression or large-file hashing during preflight."""
    sid, rid = spec['sample_id'], spec['readset_id']
    if spec['technology'] not in ('hifi','ont'):
        raise ValueError('Prepared FASTQ reuse is supported only for HiFi/ONT')
    manifest, cache_dir, run_dir, repository = map(lambda p: Path(p).resolve(),
                                                (manifest,cache_dir,run_dir,repository))
    for path in (manifest,cache_dir):
        if path == run_dir or path.is_relative_to(run_dir) or path == repository or path.is_relative_to(repository):
            raise ValueError('Prepared-input metadata and cache must be external to the new run and workflow repository')
    if cache_dir.exists() and not cache_dir.is_dir():
        raise ValueError('Prepared-input cache is not a directory')
    old, manifest_sha = metadata(manifest)
    for key in ('sample_id','readset_id','library_id'):
        if old.get(key) != spec[key]:raise ValueError('Prepared-input identity mismatch: '+key)
    if old.get('technology',spec['technology']) != spec['technology']:
        raise ValueError('Prepared-input technology mismatch')
    if Path(old.get('source','')).resolve() != Path(spec['file1']).resolve():
        raise ValueError('Prepared-input raw source path mismatch')
    if old.get('source_signature') != spec['source_signature']:
        raise ValueError('Prepared-input raw source changed after preparation')
    if old.get('filtering') != 'No quality trimming or deduplication performed':
        raise ValueError('Prepared-input filtering policy is incompatible')
    for key,size in [('source_md5',32),('source_sha256',64)]:
        if not re.fullmatch('[0-9a-f]{'+str(size)+'}',old.get(key,'')):
            raise ValueError('Prepared-input source checksum missing: '+key)
    if spec.get('expected_md5') and old['source_md5'] != spec['expected_md5']:
        raise ValueError('Prepared-input raw source MD5 mismatch')
    if old.get('source_format',spec['format']) != spec['format']:
        raise ValueError('Prepared-input source format mismatch')
    conversion=old.get('conversion','')
    if not old.get('reuse') and not conversion.startswith('samtools fastq' if spec['format']=='bam' else 'FASTQ byte copy'):
        raise ValueError('Prepared-input conversion policy is incompatible')
    parent=manifest.parent
    if parent.name != '01_inputs' or parent.parent.name != sid or parent.parent.parent.name != 'samples':
        raise ValueError('Prepared-input manifest must belong to a numbered samples/<sample>/01_inputs directory')
    expected=parent/(rid+'.fastq.gz')
    declared=Path(old.get('prepared_fastq',''))
    if declared.is_absolute():
        if declared.resolve()!=expected.resolve():raise ValueError('Prepared FASTQ path disagrees with its manifest')
    elif declared.as_posix()!=f'samples/{sid}/01_inputs/{rid}.fastq.gz':
        raise ValueError('Unsupported relative prepared FASTQ path')
    fastq=expected.resolve()
    if not fastq.is_file() or not fastq.stat().st_size or fastq.is_relative_to(run_dir) or fastq.is_relative_to(repository):
        raise ValueError('Prepared FASTQ is missing, empty or not external')
    qc_path=Path(old['reuse']['qc_manifest']) if old.get('reuse') else parent.parent/'02_read_qc'/(rid+'.json')
    qc,qc_sha=metadata(qc_path)
    if old.get('reuse') and old['reuse'].get('qc_manifest_sha256')!=qc_sha:
        raise ValueError('Inherited prepared-input QC checksum mismatch')
    if qc.get('status')!='PASS' or qc.get('num_reads',0)<1 or qc.get('num_bases',0)<1:
        raise ValueError('Prepared-input QC is not complete/PASS')
    recorded=qc.get('statistics',{}).get('file','')
    expected_relative=f'samples/{sid}/01_inputs/{rid}.fastq.gz'
    if not recorded or (Path(recorded).is_absolute() and Path(recorded).resolve()!=fastq) or (
            not Path(recorded).is_absolute() and recorded!=expected_relative):
        raise ValueError('Prepared-input QC refers to a different FASTQ')
    observed=signature(fastq)
    producer_sha=old.get('prepared_sha256')
    if producer_sha is not None:
        if not re.fullmatch('[0-9a-f]{64}',producer_sha) or not isinstance(old.get('prepared_signature'),dict):
            raise ValueError('Prepared FASTQ producer checksum/signature record is invalid')
        if old['prepared_signature'].get('size')!=observed['size']:
            raise ValueError('Prepared FASTQ size differs from its producer record')
    identity={'policy':POLICY,'sample_id':sid,'readset_id':rid,'library_id':spec['library_id'],
              'technology':spec['technology'],'prepared_fastq':str(fastq),
              'manifest':str(manifest),'manifest_sha256':manifest_sha,
              'qc':str(qc_path),'qc_sha256':qc_sha,'source_sha256':old['source_sha256']}
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    return {'identity':identity,'prepared_signature':observed,
            'manifest_signature':signature(manifest),'qc_signature':signature(qc_path),
            'producer_sha256':producer_sha,'expected_num_reads':qc['num_reads'],'expected_num_bases':qc['num_bases'],
            'producer_signature_matches':old.get('prepared_signature')==observed,
            'cache_file':str(cache_dir/(key+'.json'))}
