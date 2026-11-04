#!/usr/bin/env python3
"""Write a polyploid stage report and require one scaffold per chromosome copy."""
import argparse
from pathlib import Path
from io_utils import checksum, load, save
from report_utils import write_reports


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('metadata','assessment','delivery','manifest','markdown','html'):p.add_argument('--'+key,required=True)
    args=p.parse_args();meta=load(args.metadata);assessment=load(args.assessment);delivery=load(args.delivery)
    if assessment['assembly_id']!=delivery['assembly_id'] or assessment['sets']!=delivery['sets']:
        raise ValueError('Assessment differs from delivery')
    base=Path(args.delivery).parent;sets=[]
    for label,stat in delivery['sets'].items():
        fasta=base/f'{label}.fa';fai=Path(str(fasta)+'.fai')
        if not fai.is_file() or checksum(fasta)!=stat['sha256']:raise ValueError('Delivered FASTA/FAI changed or missing')
        sets.append({'set_id':label,'fasta':str(fasta.resolve()),'fai':str(fai.resolve()),'sha256':stat['sha256']})
    save(args.manifest,{**delivery,'species':meta['species'],'metadata':meta,'sets':sets,
         'assessment':assessment,'bundle_type':'POLYPLOID_CHROMOSOME_STAGE_REPORT'})
    text=f'# {meta["sample_id"]} — polyploid chromosome candidate\n\n'
    text+='Route: HiFi/UL raw unitigs → HiFi depth/popCNV → MAPQ0-preserving Hi-C mapping → GPhase → graph-validated FASTA.\n\n'
    text+=f'Quality: {delivery["quality_level"]}; T2T: NOT_ASSESSED; human review: NOT_REVIEWED\n\n'
    text+=f'Chromosome copies: {delivery["chromosome_counts"]}\n\n'
    for label,stat in delivery['sets'].items():text+=f'- {label}: {stat["num_contigs"]} sequences; {stat["total_bases"]} bp; N50 {stat["N50"]}; Ns {stat["N_bases"]}\n'
    text+='\n'+ '\n'.join('- '+s for s in delivery['limitations'])+'\n'
    write_reports(args.markdown,args.html,text)
    if delivery['quality_level']=='INCOMPLETE_SCAFFOLDING':
        raise ValueError('Not exactly one scaffold for each configured chromosome/haplotype group; partial products retained')


if __name__=='__main__':main()
