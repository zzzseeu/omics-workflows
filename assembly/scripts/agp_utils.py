#!/usr/bin/env python3
"""Validate original-coordinate AGP and reconstruct complete sequences."""
from collections import OrderedDict, defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans('ACGTRYMKBDHVNSW', 'TGCAYRKMVHDBNSW')


def parse_agp(path, sequences, require_complete=True):
    """Validate object tiling, component bounds, orientation, and exact source coverage."""
    objects, coverage = OrderedDict(), defaultdict(list)
    closed_objects, current = set(), None
    with Path(path).open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            c = line.rstrip().split('\t')
            if len(c) != 9:
                raise ValueError('AGP requires nine tab-separated fields')
            obj, start, end, part, kind = c[:5]
            start, end, part = int(start), int(end), int(part)
            if obj != current:
                if obj in closed_objects:
                    raise ValueError('AGP object rows must be contiguous')
                if current is not None:
                    closed_objects.add(current)
                current = obj
            rows = objects.setdefault(obj, [])
            if start != (int(rows[-1][2])+1 if rows else 1) or end < start or part != len(rows)+1:
                raise ValueError('AGP object coordinates or part numbers do not tile')
            if kind == 'W':
                component, a, b, orientation = c[5], int(c[6]), int(c[7]), c[8]
                if component not in sequences or not (1 <= a <= b <= len(sequences[component])):
                    raise ValueError('AGP component out of bounds or absent from input')
                if b-a != end-start or orientation not in {'+', '-'}:
                    raise ValueError('AGP component length/orientation is invalid')
                coverage[component].append((a, b))
            elif kind in {'N', 'U'}:
                if int(c[5]) != end-start+1 or c[7] != 'yes' or c[8] != 'proximity_ligation':
                    raise ValueError('Invalid or unsupported AGP scaffold gap')
            else:
                raise ValueError(f'Unsupported AGP component type: {kind}')
            rows.append(c)
    for cid, sequence in sequences.items():
        end = 0
        for a, b in sorted(coverage[cid]):
            if a <= end or (require_complete and a != end+1):
                raise ValueError(f'AGP duplicates or loses source sequence: {cid}')
            end = b
        if require_complete and end != len(sequence):
            raise ValueError(f'AGP loses source sequence: {cid}')
    return objects



def reconstruct(rows, sequences):
    parts = []
    for row in rows:
        if row[4] in {'N', 'U'}:
            parts.append('N' * int(row[5]))
        else:
            seq = sequences[row[5]][int(row[6])-1:int(row[7])]
            parts.append(seq if row[8] == '+' else seq.translate(COMPLEMENT)[::-1])
    return ''.join(parts)
