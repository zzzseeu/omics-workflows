#!/usr/bin/env python3
"""Summarize selected sample reports with quality and T2T status."""
import argparse
from io_utils import load, save
from pathlib import Path
from report_utils import write_reports


def batch(args):
    plan = load(args.plan)
    bundles = [load(p) for p in args.manifests]
    if [b['sample_id'] for b in bundles] != plan['sample_order']:
        raise ValueError('Batch sample list differs from selected plan')
    result = {'plan_id': plan['plan_id'], 'executor': 'local', 'status': 'STAGE_REPORTS_COMPLETE',
              'samples': [{'sample_id': b['sample_id'], 'assembly_id': b['assembly_id'],
                           'quality_level': b['quality_level'], 'T2T_status': b['T2T_status'],
                           'primary': b['assessment']['sets'].get('primary', b['assessment']['sets'].get('chromosomes')),
                           'manifest': str(Path(p).resolve())}
                          for b, p in zip(bundles, args.manifests)]}
    save(args.json, result)
    text = '# Assembly batch report\n\nSee each sample quality level and review status. T2T is not certified.\n\n'
    for sample in result['samples']:
        stat = sample['primary']
        text += f'- {sample["sample_id"]}: {sample["quality_level"]}; {stat["num_contigs"]} sequences, '
        text += f'{stat["total_bases"]:,} bp, N50 {stat["N50"]:,} bp; T2T NOT_ASSESSED\n'
    write_reports(args.markdown, args.html, text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--markdown", required=True)
    parser.add_argument("--html", required=True)
    parser.add_argument("--manifests", nargs="+", required=True)
    batch(parser.parse_args())


if __name__ == "__main__":
    main()
