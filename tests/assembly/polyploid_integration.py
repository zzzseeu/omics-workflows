#!/usr/bin/env python3
"""Real four-copy unitig assembly, copy-number pipeline and full DAG planning."""
import argparse
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import yaml

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'assembly/scripts'))
from io_utils import checksum, load, save


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output.resolve()
    if out==ROOT or out.is_relative_to(ROOT):raise ValueError('Outputs must be outside source')
    out.mkdir(parents=True,exist_ok=False);data=out/'fixture';data.mkdir()
    env=dict(os.environ,PATH=str(Path(sys.prefix)/'bin')+os.pathsep+os.defpath,
             OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    checks=[]
    def run(name,command):
        with (out/f'{name}.log').open('w') as log:
            subprocess.run(list(map(str,command)),env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
        checks.append(name);print('PASS',name,flush=True)
    rng=random.Random(20261103);ancestral=''.join(rng.choices('ACGT',k=600000));genomes=[]
    for hap in range(4):
        seq=list(ancestral)
        for i in range(hap,600000,85):
            if i%100000<25000:seq[i]=rng.choice([x for x in 'ACGT' if x!=seq[i]])
        genomes.append(''.join(seq))
    with gzip.open(data/'hifi.fq.gz','wt') as h,gzip.open(data/'ul.fq.gz','wt') as u:
        for hap,genome in enumerate(genomes):
            for i in range(1600):
                start=rng.randrange(len(genome)-15000);seq=genome[start:start+15000]
                h.write(f'@h{hap}_{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')
            for i in range(60):
                start=rng.randrange(len(genome)-100000);seq=genome[start:start+100000]
                u.write(f'@u{hap}_{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    with gzip.open(data/'R1.fq.gz','wt') as a,gzip.open(data/'R2.fq.gz','wt') as b:
        for i in range(4000):
            for mate,f in enumerate((a,b),1):
                genome=rng.choice(genomes);start=rng.randrange(len(genome)-150);seq=genome[start:start+150]
                f.write(f'@c{i}/{mate}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    (data/'samples.tsv').write_text('sample_id\tspecies\tprofile\nP\tSynthetic tetraploid\t\n')
    (data/'readsets.tsv').write_text('sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\tfile2\n'
        'P\th\tH\thifi\tfastq\thifi.fq.gz\t\nP\tu\tU\tont\tfastq\tul.fq.gz\t\n'
        'P\tc\tC\thic\tfastq\tR1.fq.gz\tR2.fq.gz\n')
    cfg={'schema_version':1,'samples_table':'samples.tsv','readsets_table':'readsets.tsv',
         'defaults':{'ploidy':4,'chromosomes_per_haplotype':1,
          'assembly':{'mode':'hifi_ul_hic','representation':'unitigs','threads':16,'mem_mb':16000,'bloom_filter_bits':0,'purge_level':0},
          'scaffolding':{'enabled':True,'backend':'gphase','threads':16,'mem_mb':32000,'processes':4}}}
    (data/'config.yaml').write_text(yaml.safe_dump(cfg))
    command=[ROOT/'run','assembly','--env-prefix',sys.prefix,'--config',data/'config.yaml',
             '--run-dir',out/'run','--cores',16,'--mem-mb',64000]
    run('complete_polyploid_dag',[*command,'--target','chromosomes','--dry-run'])
    run('native_four_copy_unitigs',[*command,'--target','contigs'])
    native=out/'run/samples/P/04_assembly/hifiasm';manifest=load(native/'manifest.json')
    assert manifest['ploidy']==4 and manifest['phase_level']=='UNPHASED_UNITIG_GRAPH'
    assert '--n-hap' in manifest['command'] and '--h1' not in manifest['command'] and '--primary' not in manifest['command']
    assert checksum(native/'unitigs.gfa')==checksum(native/manifest['source_gfa'])
    run('unitig_resume',[*command,'--target','contigs'])
    assert 'Nothing to be done' in (out/'unitig_resume.log').read_text()
    snakemake=load(out/'run/run_status.json')['command'][:-1]
    run('native_copy_numbers_and_hic_mapping',[*snakemake,
        'samples/P/05_linkage/copy_number/popcnv/06.genes.round.cn','samples/P/05_linkage/HiC.filtered.bam'])
    assert load(out/'run/samples/P/05_linkage/copy_number/manifest.json')['status']=='PASS'
    assert load(out/'run/samples/P/05_linkage/manifest.json')['retained_pairs']>0
    save(out/'validation.json',{'status':'PASS','checks':checks,'scope':'Synthetic four-copy unitig and popCNV execution; not biological phase/T2T validation'})
    print('POLYPLOID_UNITIG_CNV_VALIDATION=PASS',flush=True)


if __name__=='__main__':main()
