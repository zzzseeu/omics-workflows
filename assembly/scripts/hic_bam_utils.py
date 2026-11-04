#!/usr/bin/env python3
"""Validate adjacent primary mate pairs in Hi-C BAM files."""


def combine_haplotype_bams(inputs, output):
    """Translate both reference and mate IDs by name into a union BAM header.

    These are shared physical molecules mapped to separate copy references,
    not independent phase-informative observations. Namespace query IDs to
    avoid confusing the two reference observations of the same molecule.
    """
    import pysam
    headers = []
    for path in inputs:
        with pysam.AlignmentFile(str(path), 'rb') as bam:
            headers.append(bam.header.to_dict())
    header = dict(headers[0])
    header['SQ'] = [sq for h in headers for sq in h['SQ']]
    names = [sq['SN'] for sq in header['SQ']]
    if len(names) != len(set(names)):
        raise ValueError('Per-haplotype BAM reference names must be disjoint')
    for i, h in enumerate(headers, 1):
        if any(not sq['SN'].startswith(f'hap{i}_') for sq in h['SQ']):
            raise ValueError(f'Expected hap{i}_ contig names in per-haplotype BAM')
    header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(str(output), 'wb', header=header) as target:
        for i, path in enumerate(inputs, 1):
            with pysam.AlignmentFile(str(path), 'rb') as source:
                for read in source:
                    # Parsing SAM against the union header translates RNEXT as
                    # well as RNAME; concatenating these BAMs would corrupt IDs.
                    mapped = pysam.AlignedSegment.fromstring(read.to_string(), header)
                    name = f'hap{i}:{read.query_name}'
                    if len(name.encode()) > 254:
                        raise ValueError('Namespaced Hi-C query name exceeds BAM limit')
                    mapped.query_name = name
                    target.write(mapped)



def bam_pairs(path):
    import pysam
    count = 0
    with pysam.AlignmentFile(str(path), 'rb') as f:
        reads = iter(f)
        for a in reads:
            b = next(reads, None)
            if b is None or a.query_name != b.query_name or not a.is_paired or not b.is_paired:
                raise ValueError('Hi-C BAM must retain adjacent primary read pairs')
            if a.is_read1 == b.is_read1 or a.is_read2 == b.is_read2:
                raise ValueError('Invalid mate flags in Hi-C BAM')
            count += 1
    if not count:
        raise ValueError('No Hi-C pairs retained after filtering')
    return count
