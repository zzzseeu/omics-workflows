#!/usr/bin/env python3
"""Adapt two independent native HapHiC runs to the chromosome delivery contract."""
from pathlib import Path
import shutil
import pysam
from fasta_utils import fasta_dict, write_fasta
from hic_bam_utils import bam_pairs
from io_utils import execute, save


def separate_scaffolding(args, metadata, stage):
    options = metadata['settings']['scaffolding']
    inputs = stage/'inputs'
    inputs.mkdir(parents=True)
    fastas = [inputs/f'hap{i}.fa' for i in (1,2)]
    with fastas[0].open('w') as one, fastas[1].open('w') as two:
        with pysam.FastxFile(str(args.fasta)) as source:
            for record in source:
                hap = record.name.split('_',1)[0]
                if hap not in ('hap1','hap2'):
                    raise ValueError('Per-haplotype scaffolding requires hap1_/hap2_ contig IDs')
                write_fasta(one if hap=='hap1' else two, record.name, record.sequence)
    bams = [inputs/f'hap{i}.bam' for i in (1,2)]
    with pysam.AlignmentFile(str(args.bam), 'rb') as source:
        headers = []
        for i in (1,2):
            header = source.header.to_dict()
            header['SQ'] = [sq for sq in header['SQ'] if sq['SN'].startswith(f'hap{i}_')]
            headers.append(pysam.AlignmentHeader.from_dict(header))
        with pysam.AlignmentFile(str(bams[0]),'wb',header=headers[0]) as one, \
             pysam.AlignmentFile(str(bams[1]),'wb',header=headers[1]) as two:
            for read in source:
                hap = read.reference_name.split('_',1)[0]
                if hap not in ('hap1','hap2') or not read.next_reference_name.startswith(hap+'_'):
                    raise ValueError('Per-haplotype linkage contains an inter-haplotype or unmapped mate')
                index = 0 if hap=='hap1' else 1
                (one if index==0 else two).write(pysam.AlignedSegment.fromstring(read.to_string(),headers[index]))
    commands = []
    for i in (1,2):
        hap = f'hap{i}'
        bam_pairs(bams[i-1])
        gfa = [p for p in args.gfas if Path(p).name == hap+'.gfa']
        if len(gfa)!=1:
            raise ValueError(f'Expected one matching {hap}.gfa')
        command = ['haphic','pipeline',fastas[i-1].resolve(),bams[i-1].resolve(),
                   metadata['settings']['chromosomes_per_haplotype'],
                   '--gfa',Path(gfa[0]).resolve(),'--RE',options['restriction_sites'],
                   '--correct_nrounds',options['correct_nrounds'],'--Ns',options['gap_size'],
                   '--threads',args.threads,'--processes',min(options['processes'],args.threads),
                   '--outdir',stage/hap]
        execute(command)
        commands.append([str(x) for x in command])
    build = stage/'04.build'
    build.mkdir()
    tours = stage/'03.sort/final_tours'
    tours.mkdir(parents=True)
    for filename in ('scaffolds.fa','scaffolds.raw.agp','scaffolds.agp'):
        with (build/filename).open('w') as output:
            for i in (1,2):
                hap = f'hap{i}'
                path = stage/hap/'04.build'/filename
                if filename.endswith('.fa'):
                    for name,seq in fasta_dict(path).items():
                        write_fasta(output,hap+'_'+name,seq)
                else:
                    for line in path.read_text().splitlines():
                        if not line or line.startswith('#'):
                            continue
                        fields = line.split('\t')
                        fields[0] = hap+'_'+fields[0]
                        output.write('\t'.join(fields)+'\n')
    for i in (1,2):
        hap = f'hap{i}'
        for tour in (stage/hap/'03.sort/final_tours').glob('*.tour'):
            shutil.copyfile(tour,tours/(hap+'_'+tour.name))
    (build/'juicebox.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n'
        '# Each haplotype has its own native Juicebox helper and review coordinates.\n'
        'HERE="$(cd -- "$(dirname -- "$0")" && pwd)"\n'
        'for hap in hap1 hap2; do (cd "$HERE/../$hap/04.build" && bash juicebox.sh); done\n')
    save(stage/'separate_haplotypes.json',{'status':'AUTOMATIC_CANDIDATE','commands':commands,
         'phase_validation':'NOT_ASSESSED','shared_physical_molecules':True,
         'limitations':['Structural differences between haplotypes can bias shared Hi-C mappings',
                        'No inter-haplotype contacts were constructed; this does not prove phase correctness']})
    return commands
