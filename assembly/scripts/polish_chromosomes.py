#!/usr/bin/env python3
"""NextPolish2 with ordinary WGS k-mers and explicit missing-data records."""
import argparse
from collections import defaultdict
from pathlib import Path
import shutil
import tempfile
import subprocess
from fasta_utils import fasta_dict
from io_utils import assembly_scratch, checksum, execute, load, publish, save, source_check
from postassembly_utils import map_reads, rebase_symlinks, write_bundle, write_checksums


def prepare_wgs(meta, destination, threads):
    """QC ordinary paired WGS only; Hi-C is deliberately not an alternative."""
    destination.mkdir(parents=True,exist_ok=True);outputs=[];records=[]
    for rid,spec in meta['readsets'].items():
        if spec['technology']!='wgs':continue
        source_check(spec)
        for i in (1,2):
            if spec.get('expected_md5' if i==1 else 'expected_md5_2'):
                expected=spec['expected_md5' if i==1 else 'expected_md5_2']
                if checksum(spec[f'file{i}'],'md5')!=expected:raise ValueError('WGS source MD5 mismatch')
        r1,r2=[destination/(rid+suffix) for suffix in ('.R1.fq.gz','.R2.fq.gz')]
        execute(['fastp','-i',spec['file1'],'-I',spec['file2'],'-o',r1,'-O',r2,
                 '-5','-3','-n','0','-f','5','-F','5','-t','5','-T','5','-q','20',
                 '-w',min(threads,16),'-j',destination/(rid+'.json'),'-h',destination/(rid+'.html')])
        source_check(spec);outputs.extend([r1,r2])
        records.append({'readset_id':rid,'sources':{spec[f'file{i}']:checksum(spec[f'file{i}']) for i in (1,2)},
                        'cleaned':{p.name:checksum(p) for p in (r1,r2)}})
    save(destination/'manifest.json',{'libraries':records,'technology':'ordinary_wgs','used_for_polishing':True})
    return outputs


def audit_polished(before,after,work,threads):
    old,new=fasta_dict(before),fasta_dict(after)
    if set(old)!=set(new):raise ValueError('Polishing changed chromosome IDs or copy count')
    for name,sequence in new.items():
        if abs(len(sequence)-len(old[name]))>max(1000,len(old[name])*0.05):
            raise ValueError('Unexpected length change after polishing: '+name)
    paf=work/'polishing_changes.paf'
    with paf.open('w') as out:
        execute(['minimap2','-c','--cs=long','-x','asm5','-t',threads,before,after],stdout=out)
    intervals=defaultdict(list)
    for line in paf.read_text().splitlines():
        f=line.split('\t')
        if 'tp:A:P' not in f[12:]:continue
        if f[0]!=f[5] or f[4]!='+':
            raise ValueError('Polished sequence changes copy/strand assignment; review required')
        intervals[f[0]].append((int(f[2]),int(f[3])))
    for name,seq in new.items():
        end=covered=0
        for start,stop in sorted(intervals[name]):
            covered+=max(0,stop-max(start,end));end=max(end,stop)
        if covered<len(seq)*0.95:raise ValueError('Insufficient same-copy polishing concordance: '+name)
    return new


def polish(args):
    meta=load(args.metadata);source=load(args.finishing);options=meta['settings']['polishing']
    base=Path(args.finishing).parent
    if checksum(base/'genome.fa')!=source['sets']['genome']['sha256']:
        raise ValueError('Finishing FASTA changed after validation')
    dest=Path(args.out_dir).resolve();dest.parent.mkdir(parents=True,exist_ok=True)
    wgs=[r for r in meta['readsets'].values() if r['technology']=='wgs']
    if options['enabled'] and not wgs and options['on_missing_wgs']=='error':
        raise ValueError('NextPolish2 needs ordinary WGS reads; Hi-C is not accepted')
    with assembly_scratch(dest) as scratch:
        stage=Path(scratch)/'result';stage.mkdir();shutil.copyfile(base/'unplaced.fa',stage/'unplaced.fa')
        sequence=fasta_dict(base/'genome.fa');status='DISABLED'
        if options['enabled'] and not wgs:status='SKIPPED_MISSING_WGS'
        elif options['enabled']:
            clean=prepare_wgs(meta,stage/'wgs',args.threads)
            combined=stage/'wgs/combined.fq.gz'
            with combined.open('wb') as out:
                for p in clean:
                    with p.open('rb') as handle:shutil.copyfileobj(handle,out,8*1024*1024)
            for k in (21,31):
                # The same reads are supplied for bloom filtering and counting.
                execute(['yak','count','-k',k,'-b','30','-t',args.threads,'-o',stage/f'k{k}.yak',combined,combined])
            bam=stage/'hifi.bam';map_reads(base/'genome.fa',args.hifi,bam,args.threads)
            # Exclude ambiguous cross-copy evidence before NextPolish2's own phasing.
            unique=stage/'hifi.unique.bam'
            execute(['samtools','view','-b','-q','30','-F','2308','-o',unique,bam])
            execute(['samtools','index',unique])
            with (stage/'polished.fa').open('w') as out:
                execute(['nextPolish2','-t',args.threads,'-m','ref','-L',options['min_contig_length'],
                         '--uppercase',unique,base/'genome.fa',stage/'k21.yak',stage/'k31.yak'],stdout=out)
            sequence=audit_polished(base/'genome.fa',stage/'polished.fa',stage,args.threads)
            status='COMPLETED'
        stats=write_bundle(stage,sequence,source['groups'])
        save(stage/'manifest.json',{'sample_id':meta['sample_id'],'groups':source['groups'],'sets':stats,
             'status':status,'T2T_status':'NOT_ASSESSED','source_finishing_sha256':checksum(args.finishing),
             'source_finishing':str(Path(args.finishing).resolve()),'input_genome_sha256':source['sets']['genome']['sha256'],
             'tools':meta['tools'],'code':meta['code'],
             'parameters':options,
             'reason': 'Ordinary WGS reads are absent; Hi-C was not substituted' if status=='SKIPPED_MISSING_WGS' else None,
             'phase_policy':'MAPQ>=30 HiFi evidence; same-copy alignment screening does not independently validate phase'})
        rebase_symlinks(stage,dest);write_checksums(stage);publish(stage,dest,scratch)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('metadata','finishing','out-dir'):p.add_argument('--'+name,required=True)
    p.add_argument('--hifi',nargs='+',required=True);p.add_argument('--threads',type=int,required=True)
    polish(p.parse_args())


if __name__=='__main__':main()
