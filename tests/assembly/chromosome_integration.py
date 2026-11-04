#!/usr/bin/env python3
"""Real-tool synthetic joint assembly plus separate known-structure scaffold validation."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'assembly/scripts'))
from io_utils import save, load, checksum
from fasta_utils import write_fasta


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,PATH=str(Path(sys.prefix)/'bin')+os.pathsep+os.defpath,
             OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    assert not shutil.which('csub',path=env['PATH'])
    checks=['scheduler_free_PATH']
    def run(name,cmd,cwd=None,success=True):
        with (out/(name+'.log')).open('w') as log:
            r=subprocess.run([str(x) for x in cmd],cwd=cwd,env=env,stdout=log,
                             stderr=subprocess.STDOUT,timeout=3600)
        if (r.returncode==0)!=success:
            raise RuntimeError(f'{name}: exit {r.returncode}; see {out/(name+".log")}')
        checks.append(name)
        print('PASS',name,flush=True)
    run('unit_contracts',[sys.executable,'-m','unittest','discover','-s',ROOT/'tests/assembly','-p','test_*.py'])
    rng=random.Random(20261103)
    data=out/'joint fixture'
    data.mkdir()
    parent=''.join(rng.choices('ACGT',k=600000))
    other=list(parent)
    # Alternating homozygous/heterozygous blocks produce bubbles. Two globally
    # divergent, isolated toy contigs give hifiasm no phasing edges to solve.
    for i in range(0,len(other),75):
        if i % 100000 < 20000:
            other[i]=rng.choice([x for x in 'ACGT' if x!=other[i]])
    genomes=[parent,''.join(other)]
    with gzip.open(data/'hifi.fq.gz','wt') as f, gzip.open(data/'ul.fq.gz','wt') as ul:
        for hap,genome in enumerate(genomes,1):
            for i in range(1600):
                start=rng.randrange(len(genome)-15000+1)
                seq=genome[start:start+15000]
                f.write(f'@h{hap}_{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')
            for i in range(60):
                start=rng.randrange(len(genome)-100000+1)
                seq=genome[start:start+100000]
                ul.write(f'@u{hap}_{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    trans=str.maketrans('ACGT','TGCA')
    with gzip.open(data/'R1.fq.gz','wt') as a,gzip.open(data/'R2.fq.gz','wt') as b:
        for i in range(40000):
            genome=rng.choice(genomes)
            for mate,handle in enumerate((a,b),1):
                start=rng.randrange(len(genome)-150)
                seq=genome[start:start+150]
                if mate==2:seq=seq.translate(trans)[::-1]
                handle.write(f'@p{i}/{mate}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    (data/'samples.tsv').write_text('sample_id\tspecies\tprofile\nS\tSynthetic diploid\t\n')
    (data/'readsets.tsv').write_text('sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\tfile2\n'
        'S\thifi\tH\thifi\tfastq\thifi.fq.gz\t\nS\tul\tU\tont\tfastq\tul.fq.gz\t\n'
        'S\thic\tC\thic\tfastq\tR1.fq.gz\tR2.fq.gz\n')
    config={'schema_version':1,'samples_table':'samples.tsv','readsets_table':'readsets.tsv',
            'defaults':{'ploidy':2,'chromosomes_per_haplotype':1,
             'assembly':{'mode':'hifi_ul_hic','representation':'phased','threads':8,'mem_mb':16000,'bloom_filter_bits':0},
             'scaffolding':{'enabled':True,'threads':8,'mem_mb':16000,'processes':2}}}
    (data/'config.yaml').write_text(yaml.safe_dump(config))
    command=[ROOT/'run', 'assembly','--env-prefix',sys.prefix,'--config',data/'config.yaml',
             '--run-dir',out/'joint_run','--cores',8,'--mem-mb',64000,'--target','contigs']
    run('joint_dag',[*command,'--dry-run'])
    run('chromosome_dag',[*command[:-1],'report','--dry-run'])
    run('real_hifiasm_hifi_ul_hic',command)
    assembly=out/'joint_run/samples/S/04_assembly/hifiasm'
    manifest=load(assembly/'manifest.json')
    assert {x['set_id'] for x in manifest['sets']}=={'hap1','hap2'}
    assert all(x['num_contigs']>0 for x in manifest['sets'])
    assert '--primary' not in manifest['command'] and '--ul' in manifest['command']
    run('joint_resume',command)
    assert 'Nothing to be done' in (out/'joint_resume.log').read_text()
    # Known downstream truth: 2 haplotypes, 8 contigs per chromosome, two chromosome groups.
    fixture=out/'known_structure'
    fixture.mkdir()
    seqs={f'hap{hap}_ctg{i}':''.join(rng.choices('ACGT',k=1500000))
          for hap in (1,2) for i in range(1,9)}
    with (fixture/'combined.fa').open('w') as f,(fixture/'mapping.tsv').open('w') as mapping:
        mapping.write('contig_id\thaplotype\toriginal_id\n')
        for name,seq in seqs.items():
            write_fasta(f,name,seq)
            mapping.write(f'{name}\t{name[:4]}\t{name}\n')
    for hap in (1,2):
        with (fixture/f'hap{hap}.gfa').open('w') as gfa:
            for name,seq in seqs.items():
                if name.startswith(f'hap{hap}_'):
                    gfa.write(f'S\t{name}\t{seq}\tLN:i:{len(seq)}\trd:i:40\n')
    with gzip.open(fixture/'R1.fq.gz','wt') as a,gzip.open(fixture/'R2.fq.gz','wt') as b:
        for i in range(120000):
            hap=rng.choice((1,2))
            for mate,handle in enumerate((a,b),1):
                genome=seqs[f'hap{hap}_ctg{rng.randrange(1,9)}']
                start=rng.randrange(len(genome)-150)
                seq=genome[start:start+150]
                if mate==2:seq=seq.translate(trans)[::-1]
                handle.write(f'@k{i}/{mate}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    metadata=load(out/'joint_run/config/samples/S.metadata.json')
    metadata['settings']['scaffolding'].update(correct_nrounds=1,plot_bin_kb=500)
    save(fixture/'metadata.json',metadata)
    save(fixture/'assembly.json',{**manifest,'combined_sha256':checksum(fixture/'combined.fa')})
    from io_utils import signature
    spec={'sample_id':'S','readset_id':'known','library_id':'C','file1':str(fixture/'R1.fq.gz'),
          'file2':str(fixture/'R2.fq.gz'),'expected_md5':None,'expected_md5_2':None,
          'source_signature':signature(fixture/'R1.fq.gz'),'source_signature_2':signature(fixture/'R2.fq.gz')}
    save(fixture/'spec.json',spec)
    def task(name):
        return [sys.executable, ROOT / f'assembly/scripts/{name}.py']
    run('paired_full_parse',[*task('validate_hic_pairs'),'--spec',fixture/'spec.json','--out',fixture/'paired.json'])
    (fixture/'bwa').mkdir()
    run('bwa_index',['bwa','index','-p',fixture/'bwa/combined',fixture/'combined.fa'])
    run('real_bwa_samblaster_pair_filter',[*task('map_filter_hic'),'--metadata',fixture/'metadata.json',
        '--paired',fixture/'paired.json','--prefix',fixture/'bwa/combined','--bam',fixture/'library.bam',
        '--out',fixture/'library.json','--threads',8])
    run('query_adjacent_bam_merge',[*task('merge_validate_hic_bams'),'--inputs',fixture/'library.bam','--metrics',fixture/'library.json',
        '--bam',fixture/'merged.bam','--out',fixture/'linkage.json'])
    run('real_haphic_with_correction',[*task('scaffold_haplotypes'),'--metadata',fixture/'metadata.json','--fasta',fixture/'combined.fa',
        '--bam',fixture/'merged.bam','--gfas',fixture/'hap1.gfa',fixture/'hap2.gfa','--linkage',fixture/'linkage.json',
        '--out-dir',fixture/'haphic','--threads',8])
    delivery=[*task('export_chromosome_fasta'),'--metadata',fixture/'metadata.json','--assembly',fixture/'assembly.json',
        '--fasta',fixture/'combined.fa','--mapping',fixture/'mapping.tsv','--haphic',fixture/'haphic','--out-dir',fixture/'delivery']
    run('agp_validated_delivery',delivery)
    delivered=load(fixture/'delivery/manifest.json')
    assert delivered['chromosome_counts']=={'hap1':1,'hap2':1}, delivered
    assert delivered['source_coverage_validation']=='PASS_EXACTLY_ONCE_ALL_BASES'
    assert delivered['human_review']=='NOT_REVIEWED'
    run('real_contact_maps',[*task('plot_hic_contacts'),'--metadata',fixture/'metadata.json','--agp',fixture/'delivery/assembly.agp',
        '--bam',fixture/'merged.bam','--out-dir',fixture/'plots','--threads',8])
    report=fixture/'report'; report.mkdir()
    run('chromosome_assessment', [*task('collect_chromosome_metrics'), '--delivery',fixture/'delivery/manifest.json',
        '--plots',fixture/'plots','--assessment',report/'basic.json','--table',report/'basic.tsv'])
    run('chromosome_report',[*task('write_chromosome_report'),'--metadata',fixture/'metadata.json','--assembly',fixture/'assembly.json',
        '--delivery',fixture/'delivery/manifest.json','--plots',fixture/'plots','--assessment',report/'basic.json',
        '--manifest',report/'manifest.json','--markdown',report/'report.md','--html',report/'report.html'])
    # Stale approval must not turn a candidate into a reviewed reference.
    save(fixture/'stale.json',{'status':'APPROVED','reviewer':'synthetic-test','reviewed_at':'2026-11-03',
         'input_fasta_sha256':'wrong','automatic_agp_sha256':'wrong','review_agp_sha256':'wrong','chromosome_objects':['group1']})
    metadata['settings']['scaffolding'].update(review_agp=str(fixture/'haphic/04.build/scaffolds.raw.agp'),
                                             review_manifest=str(fixture/'stale.json'))
    save(fixture/'metadata.json',metadata)
    run('stale_review_blocked',delivery,success=False)
    save(out/'validation.json',{'status':'PASS','checks':checks,
         'limitations':'Synthetic functionality validation; not biological quality or T2T certification'})
    print('CHROMOSOME_EXTENSION_VALIDATION=PASS',flush=True)


if __name__=='__main__':main()
