#!/usr/bin/env python3
"""Estimate unitig copy number from budgeted HiFi mapping and pinned popCNV."""
import argparse
from pathlib import Path
import sys
from io_utils import assembly_scratch, checksum, execute, pipeline, publish, save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True)
    p.add_argument('--hifi',nargs='+',required=True)
    p.add_argument('--out-dir',required=True)
    p.add_argument('--threads',type=int,required=True)
    args=p.parse_args();dest=Path(args.out_dir);dest.parent.mkdir(parents=True,exist_ok=True)
    fasta=Path(args.fasta).resolve()
    gphase=Path(sys.prefix)/'opt/GPhase/source'
    with assembly_scratch(dest) as temp:
        stage=Path(temp)/'coverage';stage.mkdir();bam=stage/'reads.bam'
        mapper=max(1,(args.threads-2)//2);sorter=max(1,args.threads-mapper-2)
        pipeline([['minimap2','-t',mapper,'-I','16G','-ax','map-hifi',fasta,*[Path(x).resolve() for x in args.hifi]],
                  ['samtools','sort','-@',sorter,'-m','512M','-o',bam,'-']],stage/'mapping.stdout')
        execute(['samtools','quickcheck','-v',bam]);execute(['samtools','index','-@',args.threads,bam])
        (stage/'bam_files').mkdir();(stage/'bam_files/reads.bam').symlink_to('../reads.bam')
        (stage/'bam_files/reads.bam.bai').symlink_to('../reads.bam.bai')
        (stage/'group.list').write_text('reads\treads\n')
        names=set()
        with Path(str(fasta)+'.fai').open() as src,(stage/'gene.list').open('w') as out:
            for line in src:
                fields=line.split('\t');names.add(fields[0]);out.write(f'{fields[0]}\t1\t{fields[1]}\t{fields[0]}\n')
        command=[sys.executable,gphase/'src/popCNV/popCNV.py','-g',fasta,'-s','1000','-b','bam_files',
                 '-l','gene.list','-w','popcnv','--group','group.list','--wild','0','-t',args.threads]
        execute(command,cwd=stage)
        cn=stage/'popcnv/06.genes.round.cn';observed=set()
        with cn.open() as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():continue
                fields=line.split();observed.add(fields[0])
                value=float(fields[1])
                if not (value >= 0 and value.is_integer()):raise ValueError('Invalid or missing unitig copy-number estimate')
        if names!=observed:raise ValueError('popCNV does not cover every input unitig')
        save(stage/'manifest.json',{'status':'PASS','input_fasta_sha256':checksum(fasta),
             'copy_number_sha256':checksum(cn),'unitig_count':len(names),'command':list(map(str,command)),
             'limitations':['Read-depth copy numbers are estimates and do not validate chromosome phasing']})
        publish(stage,dest,temp)


if __name__=='__main__':main()
