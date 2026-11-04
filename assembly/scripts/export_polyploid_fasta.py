#!/usr/bin/env python3
"""Validate GPhase graph-derived FASTA/AGP, retain unassigned unitigs and export all copies."""
import argparse
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
import re
import shutil
from agp_utils import parse_agp, reconstruct
from fasta_utils import fasta_dict, fasta_statistics, write_fasta
from io_utils import assembly_scratch, checksum, execute, load, publish, save


def export(args):
    meta=load(args.metadata);assembly=load(args.assembly);base=Path(args.gphase)/'gphase_output/scaffold_hap'
    raw=fasta_dict(args.raw_unitigs)
    if checksum(args.raw_unitigs)!=assembly['fasta_sha256']:raise ValueError('Raw unitigs changed')
    # Validate both the inspectable unitig representation and graph-merged contigs.
    components={}
    objects={}
    native={}
    for representation in ('unitig','contig'):
        components[representation]=fasta_dict(base/f'gphase_final.{representation}.fasta')
        objects[representation]=parse_agp(base/f'gphase_final.{representation}.scaffold.agp',
                                         components[representation],require_complete=False)
        native[representation]=fasta_dict(base/f'gphase_final.{representation}.scaffold.fasta')
        if set(objects[representation])!=set(native[representation]):
            raise ValueError('GPhase AGP/FASTA object names disagree')
        for name,rows in objects[representation].items():
            if reconstruct(rows,components[representation])!=native[representation][name]:
                raise ValueError(f'GPhase native FASTA differs from AGP: {name}')
    used_raw=set()
    for name,seq in components['unitig'].items():
        original=re.sub(r'_dup\d+$','',name)
        if original not in raw or seq!=raw[original]:
            raise ValueError(f'GPhase unitig copy differs from its raw sequence: {name}')
        used_raw.add(original)
    chroms=meta['settings']['chromosomes_per_haplotype'];ploidy=meta['settings']['ploidy']
    expected={(c,h) for c in range(1,chroms+1) for h in range(1,ploidy+1)}
    groups=Counter();names={};unclassified=[]
    for name in native['contig']:
        match=re.match(r'^chr(\d+)g(\d+)(?:_|$)',name)
        if not match:unclassified.append(name);continue
        pair=tuple(map(int,match.groups()))
        if pair not in expected:raise ValueError(f'Unexpected chromosome/haplotype label: {name}')
        groups[pair]+=1
        names[name]=f'Chr{pair[0]:02d}_hap{pair[1]}_part{groups[pair]}'
    complete=set(groups)==expected and all(n==1 for n in groups.values()) and not unclassified
    # No selection of a largest scaffold to force the expected chromosome count.
    dest=Path(args.out_dir);dest.parent.mkdir(parents=True,exist_ok=True)
    with assembly_scratch(dest) as temp:
        stage=Path(temp)/'delivery';stage.mkdir()
        labels=['chromosomes',*[f'hap{h}.chromosomes' for h in range(1,ploidy+1)],'unplaced']
        with ExitStack() as stack:
            files={label:stack.enter_context((stage/f'{label}.fa').open('w')) for label in labels}
            for name,seq in native['contig'].items():
                if name not in names:
                    write_fasta(files['unplaced'],name,seq);continue
                hap=int(re.match(r'^chr\d+g(\d+)',name).group(1))
                output_name=names[name].removesuffix('_part1') if groups[(int(re.match(r'^chr(\d+)',name).group(1)),hap)]==1 else names[name]
                names[name]=output_name
                write_fasta(files['chromosomes'],output_name,seq)
                write_fasta(files[f'hap{hap}.chromosomes'],output_name,seq)
            for name,seq in raw.items():
                if name not in used_raw:write_fasta(files['unplaced'],'unassigned_unitig_'+name,seq)
        shutil.copyfile(base/'gphase_final.contig.fasta',stage/'components.fa')
        shutil.copyfile(base/'gphase_final.unitig.fasta',stage/'unitig_components.fa')
        shutil.copyfile(base/'gphase_final.unitig.scaffold.agp',stage/'unitig_assembly.agp')
        shutil.copyfile(base/'gphase_final.ctg2utg.txt',stage/'contig_to_unitig.txt')
        with (stage/'assembly.agp').open('w') as out:
            for name,rows in objects['contig'].items():
                for row in rows:out.write('\t'.join([names.get(name,name),*row[1:]])+'\n')
        with (stage/'chromosome_names.tsv').open('w') as out:
            out.write('gphase_object\tdelivered_id\n')
            for name,new in names.items():out.write(f'{name}\t{new}\n')
        for name in [*labels,'components','unitig_components']:
            path=stage/f'{name}.fa'
            if path.stat().st_size:execute(['samtools','faidx',path])
            else:Path(str(path)+'.fai').write_text('')
        quality='POLYPLOID_CHROMOSOME_CANDIDATE' if complete else 'INCOMPLETE_SCAFFOLDING'
        result={'sample_id':meta['sample_id'],'assembly_id':assembly['assembly_id'],'ploidy':ploidy,
             'expected_chromosomes_per_haplotype':chroms,'quality_level':quality,
             'human_review':'NOT_REVIEWED','phase_level':'GPHASE_CHROMOSOME_PHASED_UNVALIDATED',
             'T2T_status':'NOT_ASSESSED','chromosome_group_counts':{f'Chr{c:02d}_hap{h}':groups[(c,h)] for c,h in sorted(expected)},
             'chromosome_counts':{f'hap{h}':sum(n for (c,k),n in groups.items() if k==h) for h in range(1,ploidy+1)},
             'unassigned_raw_unitigs':len(set(raw)-used_raw),'unclassified_scaffolds':unclassified,
             'native_fasta_agp_concordance':'PASS_BOTH_UNITIG_AND_CONTIG_REPRESENTATIONS',
             'unitig_copy_sequence_validation':'PASS_IDENTICAL_TO_RAW_INPUT',
             'source_coverage_validation':'GRAPH_COORDINATES_BOUNDS_NONOVERLAP; RAW_GRAPH_OVERLAP_TRIMMING_ALLOWED',
             'sets':{label:fasta_statistics(stage/f'{label}.fa') for label in labels},
             'limitations':['Four copy labels are chromosome-local; no parental or C/D/K/L subgenome identity assigned',
                            'Copy number and phase accuracy require independent k-mer and sequence validation',
                            'Graph overlap trimming and duplicated collapsed unitigs follow GPhase; no claim of exactly-once raw graph bases',
                            'Unassigned raw unitigs are retained, including possible graph redundancy',
                            'Scaffolds can contain Ns; telomere and T2T completeness are not certified',
                            'BUSCO, k-mer QV, contamination and manual Hi-C curation remain required']}
        save(stage/'manifest.json',result)
        save(stage/'review.request.json',{'status':'NOT_REVIEWED','assembly_id':assembly['assembly_id'],
             'components_sha256':checksum(stage/'components.fa'),'automatic_agp_sha256':checksum(stage/'assembly.agp'),
             'instructions':'Review remapped Hi-C contacts and graph/AGP; reviewed polyploid AGP import is not implemented'})
        with (stage/'SHA256SUMS').open('w') as out:
            for path in sorted(stage.iterdir()):
                if path.is_file() and path.name!='SHA256SUMS':out.write(f'{checksum(path)}  {path.name}\n')
        publish(stage,dest,temp)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('metadata','assembly','raw-unitigs','gphase','out-dir'):p.add_argument('--'+key,required=True)
    export(p.parse_args())


if __name__=='__main__':main()
