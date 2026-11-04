#!/usr/bin/env python3
"""Check polyploid delivery hashes and record graph-specific assessment boundaries."""
import argparse
from pathlib import Path
from io_utils import checksum, load, save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('delivery','plots','json','table'):p.add_argument('--'+key,required=True)
    args=p.parse_args();delivery=load(args.delivery);base=Path(args.delivery).parent
    for label,stat in delivery['sets'].items():
        if checksum(base/f'{label}.fa')!=stat['sha256']:raise ValueError('Delivered FASTA changed after validation')
    assessment={**delivery,'evaluations':{'native_AGP_FASTA':'PASS','raw_unitig_copy_sequences':'PASS',
                'HiC_contact_maps':str(Path(args.plots).resolve()),'human_review':'NOT_REVIEWED',
                'BUSCO':'NOT_RUN','kmer_QV':'NOT_RUN','independent_phasing':'NOT_RUN','telomere_validation':'NOT_RUN'}}
    save(args.json,assessment)
    with Path(args.table).open('w') as out:
        out.write('set_id\tnum_sequences\ttotal_bases\tN50\tN_bases\n')
        for label,stat in delivery['sets'].items():
            out.write(f'{label}\t{stat["num_contigs"]}\t{stat["total_bases"]}\t{stat["N50"]}\t{stat["N_bases"]}\n')


if __name__=='__main__':main()
