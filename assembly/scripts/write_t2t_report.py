#!/usr/bin/env python3
"""Report completed post-assembly operations separately from T2T validation."""
import argparse
from io_utils import checksum,load,save
from report_utils import write_reports


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('assessment','finishing','polishing','markdown','html','manifest'):p.add_argument('--'+name,required=True)
    args=p.parse_args();a=load(args.assessment);f=load(args.finishing);pol=load(args.polishing)
    text=f'# {a["sample_id"]} post-assembly report\n\n'
    text+=f'T2T status: {a["T2T_status"]}; this is not an independently certified T2T genome.\n\n'
    text+=f'Final FASTA: {a["genome"]}\nSHA256: {a["genome_sha256"]}\n\n'
    text+=f'07 finishing: {f["accepted_changes"]} supported edits; {len(f["attempts"])} attempted loci.\n'
    text+=f'08 polishing: {pol["status"]}.\n'
    text+=f'09 assessment: {a["observed_chromosomes"]}/{a["expected_chromosomes"]} chromosomes; '
    text+=f'{a["gap_free_two_telomere_chromosomes"]} gap-free chromosomes with two terminal arrays.\n\n'
    text+='| Copy | Sequences | Bases | N50 | N bases |\n|---|---:|---:|---:|---:|\n'
    for name,stat in a['sets'].items():text+=f'| {name} | {stat["num_contigs"]} | {stat["total_bases"]} | {stat["N50"]} | {stat["N_bases"]} |\n'
    text+='\nEvaluations:\n'
    for key,value in a['evaluations'].items():
        text+=f'- {key}: {value.get("status") if isinstance(value,dict) else value}\n'
    text+='\n'+ '\n'.join('- '+s for s in a['limitations'])+'\n'
    write_reports(args.markdown,args.html,text)
    save(args.manifest,{'sample_id':a['sample_id'],'status':'COMPLETED','T2T_status':a['T2T_status'],
         'genome_sha256':a['genome_sha256'],'assessment_sha256':checksum(args.assessment),
         'finishing_sha256':checksum(args.finishing),'polishing_sha256':checksum(args.polishing)})


if __name__=='__main__':main()
