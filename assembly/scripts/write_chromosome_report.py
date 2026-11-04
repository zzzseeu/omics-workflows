#!/usr/bin/env python3
"""Write chromosome reports and block completion for incomplete scaffolding."""
import argparse
from io_utils import checksum, load, save
from pathlib import Path
from report_utils import write_reports


def chromosome_report(args):
    metadata, delivery, assembly = load(args.metadata), load(args.delivery), load(args.assembly)
    assessment = load(args.assessment)
    root = Path(args.delivery).parent
    if assessment["assembly_id"] != delivery["assembly_id"] or assessment["sets"] != delivery["sets"]:
        raise ValueError("Assessment does not match chromosome delivery")
    for label, stats in delivery["sets"].items():
        if checksum(root / f"{label}.fa") != stats["sha256"]:
            raise ValueError("Delivered FASTA changed after assessment")
    sets = [{'set_id': label, 'fasta': str((root/f'{label}.fa').resolve()),
             'fai': str((root/f'{label}.fa.fai').resolve()), 'sha256': stat['sha256']}
            for label, stat in delivery['sets'].items()]
    save(args.manifest, {**delivery, 'species': metadata['species'], 'metadata': metadata,
         'bundle_type': 'CHROMOSOME_STAGE_REPORT', 'assessment': assessment, 'sets': sets,
         'assembly_manifest': str(Path(args.assembly).resolve()), 'phase_level': assembly['phase_level'],
         'delivery_directory': str(root.resolve())})
    text = f'# {metadata["sample_id"]} — chromosome assembly\n\n'
    text += f'Route: {assembly["mode"]} → bwa/samblaster → HapHiC → AGP-validated FASTA\n\n'
    text += f'Quality: {delivery["quality_level"]}; human review: {delivery["human_review"]}; T2T: NOT_ASSESSED\n\n'
    text += f'Chromosomes per haplotype: {delivery["chromosome_counts"]}; expected {delivery["expected_chromosomes_per_haplotype"]}\n\n'
    for label, stat in delivery['sets'].items():
        text += f'- {label}: {stat["num_contigs"]} sequences; {stat["total_bases"]:,} bp; N50 {stat["N50"]:,}; Ns {stat["N_bases"]:,}\n'
    text += f'\nDelivery: {root.resolve()}\n\nContact maps: {Path(args.plots).resolve()}\n\n'
    text += '\n'.join('- '+x for x in delivery['limitations'])+'\n'
    write_reports(args.markdown, args.html, text)
    if delivery['quality_level'] == 'INCOMPLETE_SCAFFOLDING':
        raise ValueError('Chromosome counts do not match the configured target. '
                         'Partial delivery and contact maps are retained for review; '
                         'the chromosome target must not be marked complete.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--delivery", required=True)
    parser.add_argument("--plots", required=True)
    parser.add_argument("--assessment", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--markdown", required=True)
    parser.add_argument("--html", required=True)
    chromosome_report(parser.parse_args())


if __name__ == "__main__":
    main()
