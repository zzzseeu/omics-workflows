#!/usr/bin/env python3
"""Verify assessed contigs and write the sample bundle and reports."""
import argparse
from io_utils import checksum, load, save
from pathlib import Path
from report_utils import write_reports


def report(args):
    metadata, assessment = load(args.metadata), load(args.assessment)
    assembly = Path(args.assembly_dir)
    native = load(assembly / 'manifest.json')
    outputs = []
    for item in native['sets']:
        fasta = assembly / item['fasta']
        fai = assembly / item['fai']
        if checksum(fasta) != assessment['sets'][item['set_id']]['sha256']:
            raise ValueError('Assembly changed after assessment')
        if not fai.is_file():
            raise ValueError('Missing FASTA index')
        outputs.append({**item, 'fasta': str(fasta.resolve()), 'fai': str(fai.resolve()),
                        'sha256': checksum(fasta), 'fai_sha256': checksum(fai)})
    bundle = {'sample_id': metadata['sample_id'], 'species': metadata['species'],
              'assembly_id': native['assembly_id'], 'bundle_type': 'CONTIG_STAGE_REPORT',
              'quality_level': assessment['quality_level'], 'phase_level': native['phase_level'],
              'T2T_status': assessment['T2T_status'], 'sets': outputs,
              'assessment': assessment, 'metadata': metadata,
              'assembly_manifest': str((assembly / 'manifest.json').resolve()),
              'limitations': ['Primary/alternate is not two certified complete haplotypes',
                              'No chromosome scaffolding or T2T certification',
                              'Basic statistics only; no independent accuracy/contamination validation',
                              'This is a stage bundle, not a frozen final reference release']}
    save(args.manifest, bundle)
    text = f'# {metadata["sample_id"]} — HiFi-only contig report\n\n'
    text += f'Species: {metadata["species"] or "not specified"}\n\n'
    text += f'Assembly ID: {native["assembly_id"]}\n\n'
    text += f'Quality: {assessment["quality_level"]}; T2T: NOT_ASSESSED; phasing: NOT_CERTIFIED\n\n'
    text += f'HiFi input: {assessment["input_reads"]:,} reads / {assessment["input_bases"]:,} bp\n\n'
    for label, values in assessment['sets'].items():
        text += f'- {label}: {values["num_contigs"]} contigs, {values["total_bases"]:,} bp, '
        text += f'N50 {values["N50"]:,} bp, N bases {values["N_bases"]:,}\n'
    text += '\n## Limits\n\n' + '\n'.join('- ' + s for s in bundle['limitations']) + '\n'
    text += '\nDisabled evaluations: ' + ', '.join(assessment['evaluations']) + '\n'
    write_reports(args.markdown, args.html, text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--assessment", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--markdown", required=True)
    parser.add_argument("--html", required=True)
    report(parser.parse_args())


if __name__ == "__main__":
    main()
