#!/usr/bin/env python3
"""Native post-assembly validation on known two- and four-copy fixtures.

All computation belongs on a compute node. Sources and outputs are outside the
workflow repository. The fixture verifies functionality, not biological quality.
"""
import argparse
from collections import OrderedDict
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
from fasta_utils import fasta_statistics
from io_utils import checksum,save
from postassembly_utils import write_sequences


def invoke(command,log,accepted=(0,)):
    with log.open('w') as out:
        p=subprocess.run(list(map(str,command)),stdout=out,stderr=subprocess.STDOUT,timeout=3600)
    if p.returncode not in accepted:raise RuntimeError(f'{command[0]} exit {p.returncode}: {log}')
    return log.read_text()


def fixture(out,ploidy,wgs):
    rng=random.Random(20261104+ploidy);out.mkdir(parents=True,exist_ok=True)
    delivery=out/'source';delivery.mkdir(exist_ok=True)
    truth=OrderedDict();draft=OrderedDict()
    for h in range(1,ploidy+1):
        name=f'hap{h}_Chr01';seq='TTTAGGG'*30+''.join(rng.choices('ACGT',k=59580))+'CCCTAAA'*30
        truth[name]=seq
        if wgs and h==1:
            # A known consensus error away from the finishing anchors must be
            # corrected by real NextPolish2, not hidden by an identity smoke test.
            seq=seq[:10000]+next(b for b in 'ACGT' if b!=seq[10000])+seq[10001:]
        # One gap plus one terminal truncation. Other copies remain intact.
        draft[name]=seq if h!=1 else seq[2000:30000]+'N'*100+seq[31000:]
    write_sequences(out/'truth.fa',truth);write_sequences(delivery/'chromosomes.fa',draft)
    (delivery/'unplaced.fa').write_text('')
    for h in range(1,ploidy+1):
        name=f'hap{h}_Chr01';write_sequences(delivery/f'hap{h}.chromosomes.fa',{name:draft[name]})
    labels=['chromosomes',*[f'hap{h}.chromosomes' for h in range(1,ploidy+1)],'unplaced']
    save(delivery/'manifest.json',{'sample_id':'S','assembly_id':'known-fixture','quality_level':'CHROMOSOME_CANDIDATE',
         'sets':{label:fasta_statistics(delivery/(label+'.fa')) for label in labels}})
    def fastq(path,records):
        with gzip.open(path,'wt') as f:
            for name,seq in records:f.write('@'+name+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
    records=[]
    for h,(name,seq) in enumerate(truth.items(),1):
        starts=[0]*5+[len(seq)-8000]*5+[rng.randrange(0,len(seq)-8000) for _ in range(350)]
        records.extend((f'h{h}_{i}',seq[s:s+8000]) for i,s in enumerate(starts))
    fastq(out/'hifi.fq.gz',records)
    fastq(out/'ont.fq.gz',[(f'ul{h}_{i}',seq) for h,seq in enumerate(truth.values()) for i in range(5)])
    fastq(out/'hic.R1.fq.gz',[('hic/1',next(iter(truth.values()))[1000:1150])])
    fastq(out/'hic.R2.fq.gz',[('hic/2',next(iter(truth.values()))[30000:30150])])
    reads=[['S','h','hifi','hifi','fastq','hifi.fq.gz',''],['S','u','ont','ont','fastq','ont.fq.gz',''],
           ['S','c','hic','hic','fastq','hic.R1.fq.gz','hic.R2.fq.gz']]
    if wgs:
        r1=[];r2=[]
        from postassembly_utils import reverse_complement
        for h,seq in enumerate(truth.values()):
            for i in range(7000):
                start=rng.randrange(0,len(seq)-500);name=f'w{h}_{i}'
                r1.append((name+'/1',seq[start:start+150]));r2.append((name+'/2',reverse_complement(seq[start+350:start+500])))
        fastq(out/'wgs.R1.fq.gz',r1);fastq(out/'wgs.R2.fq.gz',r2)
        reads.append(['S','w','wgs','wgs','fastq','wgs.R1.fq.gz','wgs.R2.fq.gz'])
    (out/'readsets.tsv').write_text('sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\tfile2\n'+
                                  ''.join('\t'.join(r)+'\n' for r in reads))
    (out/'samples.tsv').write_text('sample_id\tspecies\tprofile\nS\tSynthetic\t\n')
    config={'schema_version':1,'samples_table':'samples.tsv','readsets_table':'readsets.tsv',
        'defaults':{'ploidy':ploidy,'chromosomes_per_haplotype':1,
            'assembly':{'mode':'hifi_ul_hic','representation':'phased' if ploidy==2 else 'unitigs',
                        'bloom_filter_bits':0,'mem_mb':32000,'purge_level':None if ploidy==2 else 0},
            'scaffolding':{'enabled':True,'backend':'haphic' if ploidy==2 else 'gphase','mem_mb':32000},
            'finishing':{'enabled':True,'source_delivery':str(delivery/'manifest.json'),
                         'gap_closing':True,'telomere_completion':True,'telomere_motif':'TTTAGGG',
                         'threads':8,'mem_mb':32000},
            'polishing':{'enabled':True,'on_missing_wgs':'skip','threads':8,'mem_mb':32000},
            'assessment':{'read_mapping':True,'kmer':True,'busco':False,'threads':8,'mem_mb':32000}}}
    (out/'config.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
    return config


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--ploidy',type=int,choices=[2,4],required=True);p.add_argument('--wgs',action='store_true')
    p.add_argument('--reuse-prepared',action='store_true')
    args=p.parse_args();out=args.output.resolve();fixture(out,args.ploidy,args.wgs)
    if args.reuse_prepared:
        import csv
        from io_utils import signature
        with (out/'readsets.tsv').open() as f:reader=csv.DictReader(f,delimiter='\t');fields=list(reader.fieldnames);rows=list(reader)
        fields.append('prepared_manifest')
        for row in rows:
            row['prepared_manifest']=''
            if row['technology'] not in ('hifi','ont'):continue
            rid=row['readset_id'];raw=out/row['file1'];old=out/'prepared_source/samples/S'
            (old/'01_inputs').mkdir(parents=True,exist_ok=True);(old/'02_read_qc').mkdir(exist_ok=True)
            spec={'sample_id':'S','readset_id':rid,'library_id':row['library_id'],'technology':row['technology'],
                  'format':'fastq','file1':str(raw),'source_signature':signature(raw),'expected_md5':None}
            source=out/(rid+'.spec.json');save(source,spec);fq=old/'01_inputs'/(rid+'.fastq.gz');manifest=fq.with_name(rid+'.json')
            invoke([sys.executable,ROOT/'assembly/scripts/prepare_long_reads.py','--spec',source,
                    '--fastq',fq,'--manifest',manifest,'--threads','2'],out/(rid+'.seed.log'))
            invoke([sys.executable,ROOT/'assembly/scripts/check_read_quality.py','--fastq',fq,'--prepared-manifest',manifest,
                    '--stats',old/'02_read_qc'/(rid+'.tsv'),'--json',old/'02_read_qc'/(rid+'.json'),
                    '--seqkit','seqkit','--threads','2'],out/(rid+'.qc.log'))
            row['prepared_manifest']=str(manifest)
        with (out/'readsets.tsv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
        cfg=yaml.safe_load((out/'config.yaml').read_text());cfg['prepared_input_cache']=str(out/'prepared_cache')
        (out/'config.yaml').write_text(yaml.safe_dump(cfg))
    command=[sys.executable,ROOT/'run','assembly','--config',out/'config.yaml','--run-dir',out/'run',
             '--target','t2t','--cores','8','--mem-mb','64000','--env-prefix',sys.prefix]
    invoke([*command,'--dry-run'],out/'dryrun.log')
    invoke(command,out/'execute.log')
    sample=out/'run/samples/S';finish=json.loads((sample/'07_finishing/manifest.json').read_text())
    if finish['accepted_changes']<2:raise ValueError('Known bridge/terminal extension was not recovered')
    if any(row['N_bases'] for row in finish['after']):raise ValueError('Known gap not closed')
    expected='COMPLETED' if args.wgs else 'SKIPPED_MISSING_WGS'
    polish=json.loads((sample/'08_polishing/manifest.json').read_text())
    if polish['status']!=expected:raise ValueError('Incorrect missing-data or polishing status')
    if args.wgs:
        from fasta_utils import fasta_dict
        truth=fasta_dict(out/'truth.fa');before=fasta_dict(sample/'07_finishing/genome.fa');after=fasta_dict(sample/'08_polishing/genome.fa')
        if before['hap1_Chr01'][10000]==truth['hap1_Chr01'][10000]:raise ValueError('Known error missing from polishing input')
        if after['hap1_Chr01'][10000]!=truth['hap1_Chr01'][10000]:raise ValueError('NextPolish2 did not correct the known error')
    assessment=json.loads((sample/'09_assessment/t2t/assessment.json').read_text())
    if assessment['gap_free_two_telomere_chromosomes']!=args.ploidy:raise ValueError('Terminal/gap screening mismatch')
    if assessment['evaluations']['inspector']['status']!='COMPLETED':raise ValueError('Inspector did not complete')
    text=invoke(command,out/'resume.log')
    if 'Nothing to be done' not in text:raise ValueError('Resume recomputed completed work')
    for path in sample.rglob('*'):
        if path.is_symlink() and not path.exists():raise ValueError('Broken published evidence link: '+str(path))
    if args.reuse_prepared:
        for rid in ('h','u'):
            if not (sample/'01_inputs'/(rid+'.fastq.gz')).is_symlink():raise ValueError('Postassembly copied prepared reads')
    save(out/'validation.json',{'status':'PASS','ploidy':args.ploidy,'wgs':args.wgs,'accepted_changes':finish['accepted_changes'],
         'screened_chromosomes':assessment['gap_free_two_telomere_chromosomes'],'polishing':expected,
         'T2T_certification':'NOT_VALIDATED_SYNTHETIC_FIXTURE','resume':'PASS','prepared_input_reuse':args.reuse_prepared})


if __name__=='__main__':main()
