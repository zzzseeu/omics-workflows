#!/usr/bin/env python3
"""Run graph-assisted polyploid scaffolding with locked inputs and explicit ploidy."""
import argparse
from pathlib import Path
from io_utils import assembly_scratch, checksum, execute, load, publish, save


def relocate_internal_links(stage, destination):
    """Rebase absolute internal links while preserving external input links."""
    stage=Path(stage);destination=Path(destination).resolve()
    for path in stage.rglob('*'):
        if path.is_symlink():
            target=path.readlink()
            if target.is_absolute() and target.is_relative_to(stage):
                new=destination/target.relative_to(stage)
                path.unlink();path.symlink_to(new)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('metadata','assembly','fasta','gfa','copies','bam','out-dir'):p.add_argument('--'+key,required=True)
    p.add_argument('--threads',type=int,required=True)
    args=p.parse_args();meta=load(args.metadata);assembly=load(args.assembly)
    if checksum(args.fasta)!=assembly['fasta_sha256'] or checksum(args.gfa)!=assembly['gfa_sha256']:
        raise ValueError('Primary unitig graph was modified after assembly')
    dest=Path(args.out_dir);dest.parent.mkdir(parents=True,exist_ok=True)
    with assembly_scratch(dest) as temp:
        stage=Path(temp)/'gphase';stage.mkdir()
        command=['gphase','pipeline','-f',Path(args.fasta).resolve(),'-g',Path(args.gfa).resolve(),
                 '-c',Path(args.copies).resolve(),'-m',Path(args.bam).resolve(),
                 '--n_chr',meta['settings']['chromosomes_per_haplotype'],
                 '--n_hap',meta['settings']['ploidy'],'-p','assembly','--rescue',
                 '--thread',max(1,args.threads-2),'--processes',min(args.threads,meta['settings']['scaffolding']['processes']),
                 '-e',meta['settings']['scaffolding']['restriction_sites']]
        execute(command,cwd=stage)
        required=['gphase_final.unitig.fasta','gphase_final.unitig.scaffold.agp',
                  'gphase_final.unitig.scaffold.fasta','gphase_final.contig.fasta',
                  'gphase_final.contig.scaffold.agp','gphase_final.contig.scaffold.fasta','gphase_final.ctg2utg.txt']
        for name in required:
            path=stage/'gphase_output/scaffold_hap'/name
            if not path.is_file() or not path.stat().st_size:raise ValueError(f'Missing GPhase output: {name}')
        if checksum(args.fasta)!=assembly['fasta_sha256'] or checksum(args.gfa)!=assembly['gfa_sha256']:
            raise ValueError('GPhase modified its raw unitig inputs')
        save(stage/'manifest.json',{'status':'AUTOMATIC_POLYPLOID_CANDIDATE','assembly_id':assembly['assembly_id'],
             'command':list(map(str,command)),'input_fasta_sha256':checksum(args.fasta),
             'input_gfa_sha256':checksum(args.gfa),'copy_number_sha256':checksum(args.copies),
             'human_review':'NOT_REVIEWED','T2T_status':'NOT_ASSESSED'})
        # GPhase makes absolute links between its stages. Preserve their meaning
        # when publishing the complete working directory from a temporary path.
        relocate_internal_links(stage,dest)
        publish(stage,dest,temp)


if __name__=='__main__':main()
