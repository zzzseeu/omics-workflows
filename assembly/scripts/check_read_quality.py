#!/usr/bin/env python3
"""Run seqkit QC or reuse identical prepared-input QC with provenance checks."""
import argparse
import csv
from pathlib import Path
from io_utils import checksum,execute,load,save,signature
from summarize_read_qc import read_qc


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('fastq','prepared-manifest','stats','json','seqkit'):p.add_argument('--'+name,required=True)
    p.add_argument('--threads',type=int,required=True);args=p.parse_args()
    prepared=load(args.prepared_manifest);reused=prepared.get('reuse')
    if reused:
        if signature(args.fastq)!=prepared['prepared_signature']:
            raise ValueError('Prepared FASTQ changed before QC reuse')
        source=reused['qc_manifest']
        if checksum(source)!=reused['qc_manifest_sha256']:raise ValueError('Source QC changed before reuse')
        qc=load(source)
        if qc['status']!='PASS':raise ValueError('Only PASS QC can be reused')
        stats={**qc['statistics'],'file':args.fastq}
        with Path(args.stats).open('w') as out:
            writer=csv.DictWriter(out,fieldnames=list(stats),delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerow(stats)
        read_qc(args)
        result=load(args.json)
        result['validation']='Reused PASS seqkit statistics for identity-checked prepared FASTQ'
        result['reuse']={'source_qc':source,'source_qc_sha256':reused['qc_manifest_sha256'],
                         'prepared_sha256':prepared['prepared_sha256']}
        save(args.json,result)
        print('READ_QC_REUSED',flush=True)
    else:
        with Path(args.stats).open('w') as out:
            execute([args.seqkit,'stats','-a','-T','-j',args.threads,args.fastq],stdout=out)
        read_qc(args)


if __name__=='__main__':main()
