#!/usr/bin/env python3
"""Map and filter one physical Hi-C library with paired-BAM validation."""
import argparse
from hic_bam_utils import bam_pairs, combine_haplotype_bams
from io_utils import checksum, execute, load, pipeline, save, source_check
from failure_scratch import native_scratch
from pathlib import Path


def map_hic(args):
    metadata, paired = load(args.metadata), load(args.paired)
    source = paired['source_spec']
    source_check(source)
    options = metadata['settings']['scaffolding']
    mode = options.get('hic_reference', 'combined')
    prefixes = args.prefix if isinstance(args.prefix, list) else [args.prefix]
    if mode not in ('combined', 'per_haplotype') or len(prefixes) != (2 if mode == 'per_haplotype' else 1):
        raise ValueError('Hi-C reference mode and BWA index prefixes disagree')
    out = Path(args.bam)
    out.parent.mkdir(parents=True, exist_ok=True)
    with native_scratch(out, prefix='.hic-') as temp:
        references, filtered_inputs = [], []
        for index, prefix in enumerate(prefixes, 1):
            label = f'hap{index}' if mode == 'per_haplotype' else 'combined'
            raw = Path(temp)/f'{label}.mapped.bam'
            filtered = Path(temp)/f'{label}.filtered.bam'
            pipeline([['bwa', 'mem', '-5SP', '-t', max(1, args.threads-3),
                       '-R', f'@RG\\tID:{source["readset_id"]}\\tSM:{metadata["sample_id"]}\\tLB:{source["library_id"]}',
                       prefix, source['file1'], source['file2']],
                      ['samblaster'], ['samtools', 'view', '-F', '3340', '-b', '-']], raw)
            pipeline([[metadata['tools']['filter_bam']['path'], raw, options['mapq'],
                       '--nm', options['nm_exclusive'], '--remove-singletons',
                       '--threads', max(1, args.threads-2)],
                      ['samtools', 'view', '-b', '-']], filtered)
            execute(['samtools', 'quickcheck', '-v', filtered])
            count = bam_pairs(filtered)
            fraction = count / paired['num_pairs']
            result = {'reference': label, 'retained_pairs': count, 'retained_fraction': fraction}
            references.append(result)
            save(Path(temp)/'retention.json', {'input_pairs': paired['num_pairs'], 'references': references})
            if fraction < options.get('min_retained_fraction', 0.001):
                raise ValueError(f'Hi-C retained only {count}/{paired["num_pairs"]} pairs '
                                 f'({fraction:.6%}) against {label}; below min_retained_fraction. '
                                 'Inspect ambiguous MAPQ=0 alignments and reference duplication before scaffolding.')
            filtered_inputs.append(filtered)
        filtered = filtered_inputs[0]
        if mode == 'per_haplotype':
            filtered = Path(temp)/'union.filtered.bam'
            combine_haplotype_bams(filtered_inputs, filtered)
        execute(['samtools', 'quickcheck', '-v', filtered])
        pairs = bam_pairs(filtered)
        if pairs != sum(r['retained_pairs'] for r in references):
            raise ValueError('Hi-C reference observation counts were not conserved')
        source_check(source)
        filtered.replace(out)
        save(args.out, {'status': 'PASS', 'readset_id': source['readset_id'],
             'input_pairs': paired['num_pairs'], 'retained_pairs': pairs,
             'retained_fraction': pairs/(paired['num_pairs']*len(prefixes)),
             'reference_pair_opportunities': paired['num_pairs']*len(prefixes),
             'bam_sha256': checksum(out),
             'hic_reference': mode, 'references': references,
             'count_unit': 'reference_pair_observations' if mode == 'per_haplotype' else 'physical_pairs',
             'shared_molecules_across_references': mode == 'per_haplotype',
             'phase_validation': 'NOT_ASSESSED',
             'mapping': 'bwa mem -5SP; samblaster per physical library; -F3340; pairwise MAPQ/NM',
             'thresholds': {'mapq_min': options['mapq'], 'NM_exclusive': options['nm_exclusive']},
             'sort_order': 'query-adjacent, not coordinate-sorted'})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--paired", required=True)
    parser.add_argument("--prefix", nargs='+', required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threads", required=True, type=int)
    map_hic(parser.parse_args())


if __name__ == "__main__":
    main()
