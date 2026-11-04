#!/usr/bin/env python3
"""GFA conversion, FASTA parsing, sequence writing and basic metrics."""
from collections import OrderedDict, defaultdict
from io_utils import checksum
from pathlib import Path


def gfa_to_fasta(gfa, fasta, required):
    count, seen = 0, set()
    with Path(fasta).open('w') as out:
        if Path(gfa).exists():
            with Path(gfa).open() as handle:
                for line in handle:
                    if not line.startswith('S\t'):
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 3 or fields[2] == '*' or fields[1] in seen:
                        raise ValueError(f'Invalid sequence-containing GFA segment in {gfa}')
                    seen.add(fields[1])
                    out.write(f'>{fields[1]}\n{fields[2]}\n')
                    count += 1
    if required and not count:
        raise ValueError('No primary contigs produced by hifiasm')
    return count



def fasta_statistics(path):
    lengths, seen, name, length, n_count, total_gc = [], set(), None, 0, 0, 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    if not length:
                        raise ValueError('Empty FASTA sequence')
                    lengths.append(length)
                name = line[1:].split()[0]
                if name in seen:
                    raise ValueError('Duplicate contig ID')
                seen.add(name)
                length = 0
            else:
                seq = line.strip().upper()
                if seq and name is None:
                    raise ValueError('FASTA sequence without a header')
                if any(base not in 'ACGTRYSWKMBDHVN' for base in seq):
                    raise ValueError('Unexpected FASTA alphabet')
                length += len(seq)
                n_count += seq.count('N')
                total_gc += seq.count('G') + seq.count('C')
    if name is not None:
        if not length:
            raise ValueError('Empty FASTA sequence')
        lengths.append(length)
    total, cumulative, n50, l50 = sum(lengths), 0, 0, 0
    for i, value in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += value
        if cumulative >= total / 2:
            n50, l50 = value, i
            break
    return {'num_contigs': len(lengths), 'total_bases': total, 'N50': n50, 'L50': l50,
            'longest_contig': max(lengths, default=0), 'N_bases': n_count,
            'GC_percent': 100 * total_gc / total if total else None,
            'sha256': checksum(path)}



def prefixed_gfa(source, target, prefix):
    with Path(source).open() as src, Path(target).open('w') as out:
        for line in src:
            fields = line.rstrip('\n').split('\t')
            if fields[0] in {'S', 'A', 'L'}:
                fields[1] = prefix + fields[1]
                if fields[0] == 'L':
                    fields[3] = prefix + fields[3]
            elif fields[0] != 'H':
                raise ValueError(f'Unsupported hifiasm GFA record: {fields[0]}')
            out.write('\t'.join(fields) + '\n')



def fasta_dict(path):
    import pysam
    result = OrderedDict()
    with pysam.FastxFile(str(path)) as f:
        for r in f:
            if r.name in result or not r.sequence:
                raise ValueError('Duplicate ID or empty FASTA sequence')
            result[r.name] = r.sequence.upper()
    return result



def write_fasta(handle, name, sequence):
    handle.write(f'>{name}\n')
    for i in range(0, len(sequence), 60):
        handle.write(sequence[i:i+60] + '\n')
