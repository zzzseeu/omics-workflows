#!/usr/bin/env python3
"""Map a physical Hi-C library to raw unitigs while retaining MAPQ0 for GPhase."""
import argparse
from pathlib import Path
import tempfile
from hic_bam_utils import bam_pairs
from io_utils import checksum, execute, load, pipeline, save, source_check


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('paired','prefix','bam','out'):p.add_argument('--'+key,required=True)
    p.add_argument('--threads',type=int,required=True)
    args=p.parse_args();paired=load(args.paired);source=paired['source_spec'];source_check(source)
    destination=Path(args.bam);destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.map-',dir=destination.parent) as temp:
        bam=Path(temp)/'mapped.bam'
        pipeline([['bwa','mem','-5SP','-t',max(1,args.threads-3),'-R',
                  f'@RG\\tID:{source["readset_id"]}\\tSM:{source["sample_id"]}\\tLB:{source["library_id"]}',
                  args.prefix,source['file1'],source['file2']],['samblaster'],
                  ['samtools','view','-F','3340','-b','-']],bam)
        execute(['samtools','quickcheck','-v',bam]);pairs=bam_pairs(bam);source_check(source)
        bam.replace(destination)
        save(args.out,{'status':'PASS','retained_pairs':pairs,'input_pairs':paired['num_pairs'],
                      'bam_sha256':checksum(destination),'readset_id':source['readset_id'],
                      'mapping':'bwa mem -5SP; samblaster per physical library; -F3340; MAPQ0 retained; no NM filter',
                      'sort_order':'query-adjacent'})


if __name__=='__main__':main()
