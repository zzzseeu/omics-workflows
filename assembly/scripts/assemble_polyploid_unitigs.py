#!/usr/bin/env python3
"""Assemble HiFi/UL into an untouched lossless primary-unitig graph for GPhase."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from fasta_utils import gfa_to_fasta
from io_utils import assembly_scratch, checksum, execute, load, publish, save


def assemble(args):
    metadata = load(args.metadata)
    options = metadata['settings']['assembly']
    if any(load(path)['status'] != 'PASS' for path in args.qcs):
        raise ValueError('Read QC failed')
    if load(args.audit)['status'] != 'PASS':
        raise ValueError('Long-read ID audit failed')
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temp:
        stage = Path(temp) / 'assembly'
        stage.mkdir()
        command = ['hifiasm', '-o', stage/'asm', '-t', args.threads, '-f', options['bloom_filter_bits'],
                   '--n-hap', metadata['settings']['ploidy'], '-l', 0, '--ul', ','.join(args.ul), *args.hifi]
        execute(command)
        graphs = [p for p in (stage/'asm.bp.p_utg.gfa', stage/'asm.p_utg.gfa') if p.is_file()]
        if len(graphs) != 1:
            raise ValueError(f'Expected one native primary-unitig graph, found: {graphs}')
        shutil.copyfile(graphs[0], stage/'unitigs.gfa')
        count = gfa_to_fasta(stage/'unitigs.gfa', stage/'unitigs.fa', True)
        identity = {'metadata': metadata, 'inputs': [load(p) for p in args.manifests]}
        aid = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(stage/'manifest.json', {'sample_id': metadata['sample_id'], 'assembly_id': aid,
             'mode': 'polyploid_hifi_ul_unitigs', 'representation': 'unitigs',
             'phase_level': 'UNPHASED_UNITIG_GRAPH', 'ploidy': metadata['settings']['ploidy'],
             'command': list(map(str, command)), 'source_gfa': graphs[0].name,
             'unitig_count': count, 'gfa_sha256': checksum(stage/'unitigs.gfa'),
             'fasta_sha256': checksum(stage/'unitigs.fa'),
             'unitig_policy': 'Native p_utg IDs/order/sequences unchanged; no contig purge, polishing or scaffolding',
             'metadata': metadata})
        publish(stage, destination, temp)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('metadata','audit','out-dir'):p.add_argument('--'+key, required=True)
    for key in ('hifi','ul','qcs','manifests'):p.add_argument('--'+key,nargs='+',required=True)
    p.add_argument('--threads',type=int,required=True)
    assemble(p.parse_args())


if __name__=='__main__':main()
