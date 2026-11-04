#!/usr/bin/env python3
"""Remap Hi-C to graph-merged components and draw coordinate-consistent contact maps."""
import argparse
from pathlib import Path
from types import SimpleNamespace
from io_utils import assembly_scratch, execute, load, publish, save
from map_filter_hic import map_hic
from merge_validate_hic_bams import merge_hic


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('metadata','fasta','agp','out-dir'):p.add_argument('--'+key,required=True)
    p.add_argument('--paired',nargs='+',required=True)
    p.add_argument('--threads',type=int,required=True)
    args=p.parse_args();meta=load(args.metadata);dest=Path(args.out_dir);dest.parent.mkdir(parents=True,exist_ok=True)
    with assembly_scratch(dest) as temp:
        stage=Path(temp)/'contact_maps';stage.mkdir();(stage/'bwa').mkdir()
        prefix=(stage/'bwa/components').resolve()
        execute(['bwa','index','-p',prefix,Path(args.fasta).resolve()])
        libraries=[];metrics=[]
        for i,paired in enumerate(args.paired,1):
            bam=stage/f'library{i}.bam';record=stage/f'library{i}.json'
            map_hic(SimpleNamespace(metadata=args.metadata,paired=paired,prefix=str(prefix),
                                    bam=str(bam),out=str(record),threads=args.threads))
            libraries.append(str(bam));metrics.append(str(record))
        bam=stage/'HiC.remapped.bam'
        merge_hic(SimpleNamespace(inputs=libraries,metrics=metrics,bam=str(bam),out=str(stage/'mapping.json')))
        execute(['haphic','plot',Path(args.agp).resolve(),bam.resolve(),
                 '--bin_size',meta['settings']['scaffolding']['plot_bin_kb'],'--min_len','1',
                 '--separate_plots','--normalization','log10','--threads',args.threads],cwd=stage)
        for name in ('contact_map.pdf','separate_plots.pdf'):
            if not (stage/name).is_file():raise ValueError(f'Missing contact map: {name}')
        save(stage/'manifest.json',{'status':'PASS','mapping':'Hi-C remapped to GPhase final component FASTA; MAPQ>=1 and NM<3 pairwise',
             'coordinate_system':'Delivery contig assembly AGP and components.fa',
             'human_review':'NOT_REVIEWED','retained_mapping':'HiC.remapped.bam'})
        publish(stage,dest,temp)


if __name__=='__main__':main()
