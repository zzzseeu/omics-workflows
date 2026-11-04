#!/usr/bin/env python3
"""Run HapHiC transactionally, validate outputs and repair published helper paths."""
import argparse
import os
from io_utils import checksum, execute, load, publish, save
from failure_scratch import native_scratch
from pathlib import Path
from scaffold_haplotypes_separately import separate_scaffolding


def scaffold(args):
    metadata = load(args.metadata)
    options = metadata['settings']['scaffolding']
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with native_scratch(destination, prefix='.haphic-') as temp:
        stage = Path(temp) / 'haphic'
        mode = load(args.linkage).get('hic_reference', 'combined')
        command = ['haphic', 'pipeline', Path(args.fasta).resolve(), Path(args.bam).resolve(),
                   2 * metadata['settings']['chromosomes_per_haplotype'],
                   '--gfa', ','.join(str(Path(p).resolve()) for p in args.gfas),
                   '--phasing_weight', '1', '--RE', options['restriction_sites'],
                   '--correct_nrounds', options['correct_nrounds'], '--Ns', options['gap_size'],
                   '--threads', args.threads, '--processes', min(options['processes'], args.threads),
                   '--outdir', stage]
        if mode == 'per_haplotype':
            commands = separate_scaffolding(args, metadata, stage)
        else:
            execute(command)
            commands = [[str(x) for x in command]]
        for name in ('scaffolds.fa', 'scaffolds.raw.agp', 'scaffolds.agp', 'juicebox.sh'):
            if not (stage / '04.build' / name).is_file():
                raise ValueError(f'Missing HapHiC build artifact: {name}')
        save(stage / 'manifest.json', {'status': 'AUTOMATIC_CANDIDATE',
             'command': commands[0] if len(commands)==1 else None, 'settings': options,
             'commands_executed': commands,
             'input_fasta_sha256': checksum(args.fasta),
             'input_bam_sha256': load(args.linkage)['bam_sha256'],
             'hic_reference': load(args.linkage).get('hic_reference', 'combined'),
             'phase_validation': 'NOT_ASSESSED',
             'automatic_correction_rounds_requested': options['correct_nrounds'],
             'human_review': 'NOT_REVIEWED'})
        # HapHiC embeds temporary paths in Juicebox helper. Replace with published paths.
        for helper in stage.rglob('*.sh'):
            if not helper.is_symlink():
                helper.write_text(helper.read_text().replace(str(stage), str(destination.resolve())))
        for path in stage.rglob('*'):
            if path.is_symlink():
                target = os.readlink(path)
                if target.startswith(str(stage)):
                    path.unlink()
                    path.symlink_to(target.replace(str(stage), str(destination.resolve()),1))
        publish(stage, destination, temp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--linkage", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--gfas", nargs="+", required=True)
    scaffold(parser.parse_args())


if __name__ == "__main__":
    main()
