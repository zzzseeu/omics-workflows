#!/usr/bin/env python3
"""Shared provenance, alignment and terminal-repeat checks for post-assembly."""
from collections import OrderedDict
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from fasta_utils import fasta_dict, fasta_statistics, write_fasta
from io_utils import checksum, execute, load, pipeline, save


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


def write_sequences(path, sequences):
    with Path(path).open('w') as handle:
        for name, seq in sequences.items():
            write_fasta(handle,name,seq)


def terminal_repeats(sequence, motif, window=10000, distance=1000, repeats=10):
    """Require a terminal array near the real sequence end, not an internal hit.

    Exact cyclic rotations in either orientation are searched. Interrupted arrays
    are intentionally unresolved, even if tidk reports a nearby repeat signal.
    """
    motifs={motif[i:]+motif[:i] for i in range(len(motif))}
    rc=reverse_complement(motif)
    motifs|={rc[i:]+rc[:i] for i in range(len(rc))}
    result={}
    for side, region in [('left',sequence[:window]),('right',sequence[-window:])]:
        hits=[]
        for m in motifs:
            for match in re.finditer('(?:'+m+'){'+str(repeats)+',}',region):
                offset=match.start() if side=='left' else len(region)-match.end()
                if offset<=distance:
                    hits.append({'distance':offset,'array_bases':match.end()-match.start(),'motif':m})
        result[side]=max(hits,key=lambda h:h['array_bases']) if hits else None
    return result


def inspect_sequences(sequences, options, motif):
    result=[]
    for name,seq in sequences.items():
        ends=terminal_repeats(seq,motif,options['terminal_window'],
                              options['terminal_max_distance'],options['terminal_min_repeats'])
        gaps=[{'start':m.start(),'end':m.end(),'length':m.end()-m.start()}
              for m in re.finditer('N+',seq)]
        ambiguity=sum(base not in 'ACGT' for base in seq)
        result.append({'sequence_id':name,'length':len(seq),'gaps':gaps,'N_bases':seq.count('N'),
                       'ambiguous_bases':ambiguity,'telomeres':ends,
                       'gap_free_two_telomeres':not ambiguity and bool(ends['left']) and bool(ends['right'])})
    return result


def verify_delivery(path, sample):
    manifest=load(path);base=Path(path).parent
    if manifest.get('sample_id')!=sample or manifest.get('quality_level')=='INCOMPLETE_SCAFFOLDING':
        raise ValueError('Source delivery has a different sample or incomplete scaffolding')
    sets=manifest.get('sets',{})
    if 'chromosomes' not in sets:
        raise ValueError('Source delivery lacks chromosome FASTA provenance')
    for label,stat in sets.items():
        if checksum(base/(label+'.fa'))!=stat['sha256']:
            raise ValueError('Source delivery FASTA checksum mismatch: '+label)
    sequences=fasta_dict(base/'chromosomes.fa')
    if not sequences:
        raise ValueError('No chromosome sequences in source delivery')
    groups={}
    for label in sets:
        if re.fullmatch(r'hap\d+\.chromosomes',label):
            names=list(fasta_dict(base/(label+'.fa')))
            for name in names:
                if name in groups or name not in sequences:
                    raise ValueError('Invalid/overlapping haplotype sets')
                groups[name]=label.split('.')[0]
    if set(groups)!=set(sequences):
        raise ValueError('Source delivery does not assign every chromosome to one copy')
    # Validate actual sequence identity as well as independent set hashes.
    for hap in set(groups.values()):
        for name,seq in fasta_dict(base/(hap+'.chromosomes.fa')).items():
            if sequences[name]!=seq:
                raise ValueError('Copy FASTA disagrees with combined FASTA')
    return manifest,sequences,groups


def map_reads(reference, reads, output, threads, preset='map-hifi', min_mapq=None):
    """Two concurrent processes share the advertised thread budget."""
    reference,output=Path(reference).resolve(),Path(output).resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    commands=[['minimap2','-a','-x',preset,'--secondary=yes','-N','10','-t',max(1,threads-2),reference,
               *[str(Path(p).resolve()) for p in reads]],
              ['samtools','sort','-@','1','-m','512M','-o',output,'-']]
    # minimap2 treats multiple FASTQ arguments as paired fragments. Long-read
    # libraries are therefore streamed as ONE concatenated input, never as mates.
    if len(reads)>1:
        stream=('import sys,gzip,shutil\n'
                'for p in sys.argv[1:]:\n'
                ' with open(p,"rb") as h: magic=h.read(2)\n'
                ' opener=gzip.open if magic==bytes([31,139]) else open\n'
                ' with opener(p,"rb") as h: shutil.copyfileobj(h,sys.stdout.buffer,8388608)\n')
        commands[0]=commands[0][:-len(reads)]+['-']
        commands[0][commands[0].index('-t')+1]=max(1,threads-3)
        commands.insert(0,[sys.executable,'-c',stream,*[str(Path(p).resolve()) for p in reads]])
    sink=output.with_suffix('.pipeline.stdout')
    pipeline(commands,sink);sink.unlink()
    execute(['samtools','quickcheck','-v',output])
    execute(['samtools','index',output])


def supported_alignments(bam, chrom, start, end, min_mapq, min_identity):
    """Primary, unique mappings with a reliable aligned base fraction."""
    seen=set()
    for record in bam.fetch(chrom,max(0,start),end):
        if record.is_unmapped or record.is_secondary or record.is_supplementary or record.mapping_quality<min_mapq:
            continue
        aligned=record.query_alignment_length or 0
        if not aligned or not record.has_tag('NM') or 1-record.get_tag('NM')/aligned<min_identity:
            continue
        if record.query_name in seen:
            raise ValueError('Duplicate primary long-read ID: '+record.query_name)
        seen.add(record.query_name)
        yield record


def write_bundle(stage,sequences,groups):
    write_sequences(stage/'genome.fa',sequences)
    execute(['samtools','faidx',stage/'genome.fa'])
    stats={'genome':fasta_statistics(stage/'genome.fa')}
    for hap in sorted(set(groups.values())):
        write_sequences(stage/(hap+'.fa'),OrderedDict((n,s) for n,s in sequences.items() if groups[n]==hap))
        execute(['samtools','faidx',stage/(hap+'.fa')])
        stats[hap]=fasta_statistics(stage/(hap+'.fa'))
    return stats


def write_checksums(stage):
    with (stage/'SHA256SUMS').open('w') as out:
        for f in sorted(stage.rglob('*')):
            if f.is_file() and not f.is_symlink() and f.name!='SHA256SUMS':
                out.write(f'{checksum(f)}  {f.relative_to(stage)}\n')


def rebase_symlinks(stage,destination):
    """Keep native evidence links valid after transactional directory publication."""
    stage=Path(stage).resolve();links={}
    for path in stage.rglob('*'):
        if not path.is_symlink():continue
        target=path.resolve()
        if target.is_relative_to(stage):
            target=Path(destination).resolve()/target.relative_to(stage)
            path.unlink();path.symlink_to(target)
        links[str(path.relative_to(stage))]=str(target)
    if links:save(stage/'symlinks.json',links)


def native_tidk(genome, outdir, motif):
    outdir=Path(outdir).resolve();outdir.mkdir(parents=True,exist_ok=True)
    execute(['tidk','search','--string',motif,'--window','1000','--output','repeats',
             '--dir',outdir,Path(genome).resolve()])
    if not list(outdir.glob('repeats*')):
        raise ValueError('tidk completed without repeat output')
