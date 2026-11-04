#!/usr/bin/env python3
"""Select chromosome objects and publish global and per-chromosome contact maps."""
import argparse
import tempfile
from io_utils import execute, load, publish
from pathlib import Path


def contact_plot(args):
    options = load(args.metadata)['settings']['scaffolding']
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.plot-', dir=destination.parent) as temp:
        stage = Path(temp)/'contact_maps'
        stage.mkdir()
        chromosome_objects = []
        with Path(args.agp).open() as f:
            for line in f:
                obj = line.split('\t', 1)[0]
                if obj.startswith(('hap1_Chr', 'hap2_Chr')) and obj not in chromosome_objects:
                    chromosome_objects.append(obj)
        # Original contig coordinates must match the BAM, including correction splits.
        command = ['haphic', 'plot', Path(args.agp).resolve(), Path(args.bam).resolve(),
                   '--bin_size', options['plot_bin_kb'], '--min_len', '1',
                   '--separate_plots', '--normalization', 'log10', '--threads', args.threads]
        if chromosome_objects:
            command += ['--specified_scaffolds', ','.join(chromosome_objects)]
        execute(command, cwd=stage)
        for name in ('contact_map.pdf', 'separate_plots.pdf'):
            if not (stage/name).is_file():
                raise ValueError(f'Contact plot missing: {name}')
        publish(stage, destination, temp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--agp", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    contact_plot(parser.parse_args())


if __name__ == "__main__":
    main()
