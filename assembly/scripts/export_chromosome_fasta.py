#!/usr/bin/env python3
"""Validate AGP coverage and review provenance, then export chromosome FASTA."""
import argparse
import csv
import shutil
import tempfile
from agp_utils import parse_agp, reconstruct
from contextlib import ExitStack, contextmanager
from fasta_utils import fasta_dict, fasta_statistics, write_fasta
from io_utils import checksum, execute, load, publish, save
from pathlib import Path


def deliver(args):
    metadata, assembly = load(args.metadata), load(args.assembly)
    options = metadata['settings']['scaffolding']
    sequences = fasta_dict(args.fasta)
    if checksum(args.fasta) != assembly['combined_sha256']:
        raise ValueError('Joint assembly FASTA differs from its source manifest')
    haps = {}
    with Path(args.mapping).open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            haps[row['contig_id']] = row['haplotype']
    if set(haps) != set(sequences):
        raise ValueError('Haplotype table does not match source FASTA')
    raw_agp = Path(args.haphic) / '04.build/scaffolds.raw.agp'
    automatic = parse_agp(raw_agp, sequences)
    native = fasta_dict(Path(args.haphic) / '04.build/scaffolds.fa')
    if set(native) != set(automatic):
        raise ValueError('Native HapHiC FASTA and AGP objects differ')
    for obj, rows in automatic.items():
        if reconstruct(rows, sequences) != native[obj]:
            raise ValueError(f'Native HapHiC FASTA disagrees with AGP: {obj}')
    # Only objects derived from explicit HapHiC tours are anchored chromosome candidates.
    tours = list((Path(args.haphic) / '03.sort/final_tours').glob('*.tour'))
    anchored = {p.stem.rsplit('_', 1)[0] for p in tours}
    if not anchored or not anchored <= set(automatic):
        raise ValueError('Cannot trace HapHiC chromosome objects to sorted tours')
    reviewed = False
    chosen_agp, objects = raw_agp, automatic
    if options['review_agp']:
        review = load(options['review_manifest'])
        required = {'status', 'reviewer', 'reviewed_at', 'input_fasta_sha256', 'automatic_agp_sha256',
                    'review_agp_sha256', 'chromosome_objects'}
        if not required <= set(review) or review['status'] != 'APPROVED' or not review['reviewer'] or not review['reviewed_at']:
            raise ValueError('Manual review manifest is incomplete or not APPROVED')
        if (review['input_fasta_sha256'] != checksum(args.fasta) or
            review['automatic_agp_sha256'] != checksum(raw_agp) or
            review['review_agp_sha256'] != checksum(options['review_agp'])):
            raise ValueError('Manual review is stale or refers to a different assembly/AGP')
        chosen_agp = Path(options['review_agp'])
        objects = parse_agp(chosen_agp, sequences)
        anchored = set(review['chromosome_objects'])
        if not anchored or not anchored <= set(objects):
            raise ValueError('Reviewed chromosome object list is invalid')
        reviewed = True
    assigned, unresolved = {'hap1': [], 'hap2': []}, []
    for obj, rows in objects.items():
        hapset = {haps[r[5]] for r in rows if r[4] == 'W'}
        if obj in anchored and len(hapset) == 1:
            assigned[next(iter(hapset))].append(obj)
        else:
            unresolved.append(obj)
    expected = metadata['settings']['chromosomes_per_haplotype']
    counts = {hap: len(v) for hap, v in assigned.items()}
    quality = ('CHROMOSOME_REVIEWED' if reviewed else 'CHROMOSOME_CANDIDATE') if all(
        count == expected for count in counts.values()) else 'INCOMPLETE_SCAFFOLDING'
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.delivery-', dir=destination.parent) as temp:
        stage = Path(temp) / 'delivery'
        stage.mkdir()
        names = {}
        for hap, group in assigned.items():
            # Local names reflect length order, not syntenic correspondence to reference chromosomes.
            group.sort(key=lambda o: int(objects[o][-1][2]), reverse=True)
            names.update({o: f'{hap}_Chr{i:02d}' for i, o in enumerate(group, 1)})
        names.update({o: f'unplaced_{i:06d}' for i, o in enumerate(unresolved, 1)})
        with ExitStack() as stack:
            handles = {k: stack.enter_context((stage / f'{k}.fa').open('w'))
                       for k in ('chromosomes', 'hap1.chromosomes', 'hap2.chromosomes', 'unplaced')}
            agp = stack.enter_context((stage / 'assembly.agp').open('w'))
            chr_agp = stack.enter_context((stage / 'chromosomes.agp').open('w'))
            table = stack.enter_context((stage / 'chromosome_names.tsv').open('w'))
            table.write('sequence_id\tsource_object\thaplotype\tstatus\tlength\n')
            for obj, name in names.items():
                rows = objects[obj]
                sequence = reconstruct(rows, sequences)
                hap = next((h for h in assigned if obj in assigned[h]), None)
                if hap:
                    write_fasta(handles['chromosomes'], name, sequence)
                    write_fasta(handles[hap+'.chromosomes'], name, sequence)
                else:
                    write_fasta(handles['unplaced'], name, sequence)
                for row in rows:
                    line = '\t'.join([name, *row[1:]]) + '\n'
                    agp.write(line)
                    if hap:
                        chr_agp.write(line)
                table.write(f'{name}\t{obj}\t{hap or "unresolved"}\t{"chromosome_candidate" if hap else "unplaced"}\t{len(sequence)}\n')
        for name in ('chromosomes', 'hap1.chromosomes', 'hap2.chromosomes', 'unplaced'):
            path = stage / f'{name}.fa'
            if path.stat().st_size:
                execute(['samtools', 'faidx', path])
            else:
                Path(str(path)+'.fai').write_text('')
            with (stage / f'{name}.chrom.sizes').open('w') as out:
                for line in Path(str(path)+'.fai').read_text().splitlines():
                    out.write('\t'.join(line.split('\t')[:2])+'\n')
        save(stage / 'review.request.json', {'status': 'APPROVED' if reviewed else 'PENDING',
             'reviewer': None, 'reviewed_at': None, 'input_fasta_sha256': checksum(args.fasta),
             'automatic_agp_sha256': checksum(raw_agp), 'review_agp_sha256': None,
             'chromosome_objects': sorted(anchored),
             'instructions': 'Review contact maps and original-coordinate AGP; supply external review_agp and a completed manifest. Do not approve by changing status alone.'})
        if reviewed:
            shutil.copyfile(options['review_manifest'], stage / 'review.applied.json')
        result = {'sample_id': metadata['sample_id'], 'assembly_id': assembly['assembly_id'],
             'quality_level': quality, 'human_review': 'APPROVED' if reviewed else 'NOT_REVIEWED',
             'T2T_status': 'NOT_ASSESSED', 'chromosome_counts': counts,
             'hic_reference': load(Path(args.haphic)/'manifest.json').get('hic_reference', 'combined'),
             'phase_validation': 'NOT_ASSESSED',
             'expected_chromosomes_per_haplotype': expected, 'unplaced_objects': len(unresolved),
             'mixed_haplotype_objects': [o for o in anchored if len({haps[r[5]] for r in objects[o] if r[4]=='W'}) > 1],
             'input_fasta_sha256': checksum(args.fasta), 'applied_agp_sha256': checksum(chosen_agp),
             'source_coverage_validation': 'PASS_EXACTLY_ONCE_ALL_BASES',
             'native_fasta_agp_concordance': 'PASS',
             'sets': {name: fasta_statistics(stage/f'{name}.fa') for name in
                      ('chromosomes', 'hap1.chromosomes', 'hap2.chromosomes', 'unplaced')},
             'limitations': ['Chromosome numbers are length ranks, not reference chromosome identities',
                            'Hi-C gaps contain Ns; this is not a T2T or gap-free assembly',
                            'Hi-C phasing is not independent phase validation; hap labels are not parents',
                            'BUSCO, k-mer QV, contamination and telomere completeness are not assessed']}
        if result['hic_reference'] == 'per_haplotype':
            result['limitations'].append('Hi-C molecules were reused within each haplotype reference; '
                                        'absence of inter-haplotype contacts is by construction, '
                                        'not evidence of correct phasing')
        save(stage/'manifest.json', result)
        with (stage/'SHA256SUMS').open('w') as out:
            for p in sorted(stage.iterdir()):
                if p.is_file() and p.name != 'SHA256SUMS':
                    out.write(f'{checksum(p)}  {p.name}\n')
        publish(stage, destination, temp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--haphic", required=True)
    parser.add_argument("--out-dir", required=True)
    deliver(parser.parse_args())


if __name__ == "__main__":
    main()
