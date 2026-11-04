#!/usr/bin/env python3
"""Run post-assembly assessments and produce an evidence-bound T2T review packet."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from fasta_utils import fasta_dict
from io_utils import assembly_scratch, checksum, execute, load, publish, save
from polish_chromosomes import prepare_wgs
from postassembly_utils import inspect_sequences, map_reads, native_tidk, rebase_symlinks, write_checksums


def run_merqury(genome,copy_fastas,clean,work,threads,mem_mb):
    work.mkdir();database=work/'reads.meryl'
    execute(['meryl','k=21',f'threads={threads}',f'memory={max(1,mem_mb//2000)}','count','output',database,*clean])
    # Merqury's shell scripts do not propagate a command-line CPU budget.
    # A per-attempt shim caps each meryl process; two can coexist in a pipe.
    tools=work/'bin';tools.mkdir();meryl=Path(shutil.which('meryl')).resolve()
    shim=tools/'meryl'
    shim.write_text('#!/usr/bin/env python3\nimport os,sys\nos.execv('+repr(str(meryl))+', ['+
                    repr(str(meryl))+', '+repr(f'threads={max(1,threads//2)}')+', '+
                    repr(f'memory={max(1,mem_mb//4000)}')+', *sys.argv[1:]])\n');shim.chmod(0o755)
    env=dict(os.environ,MERQURY=str(Path(shutil.which('merqury.sh')).resolve().parent),
             PATH=str(tools)+os.pathsep+os.environ['PATH'],OMP_NUM_THREADS=str(max(1,threads//2)))
    records={}
    for label,path in [('genome',genome),*copy_fastas.items()]:
        d=work/label;d.mkdir()
        execute(['merqury.sh',database,path,'assessment'],cwd=d,env=env)
        qv=d/'assessment.qv';completeness=d/'assessment.completeness.stats'
        if not qv.is_file() or not completeness.is_file():
            raise ValueError('Merqury did not produce QV/completeness; inspect native logs')
        values=qv.read_text().splitlines()
        if len(values)!=1:raise ValueError('Unexpected Merqury QV record count')
        fields=values[0].split('\t');value=float(fields[3])
        if math.isnan(value):raise ValueError('Merqury QV is undefined')
        records[label]={'qv':'Infinity' if math.isinf(value) else value,
                        'completeness_record':completeness.read_text().strip(),
                        'qv_file_sha256':checksum(qv)}
    return records


def run_inspector(genome,hifi,work,threads):
    work.mkdir()
    map_reads(genome,hifi,work/'read_to_contig.bam',threads)
    execute(['inspector.py','-c',genome,'-r',*hifi,'-o',str(work)+'/',
             '-d','hifi','-t',threads,'--skip_read_mapping','--noplot',
             '--min_contig_length','1','--min_contig_length_assemblyerror','10000'])
    summary=work/'summary_statistics'
    if not summary.is_file():raise ValueError('Inspector summary missing')
    text=summary.read_text()
    required=['Read to Contig alignment:','Total small-scale assembly error\t','Structural error\t']
    if any(key not in text for key in required):
        raise ValueError('Inspector returned an incomplete assessment; inspect native stderr')
    if not (work/'Inspector.log').is_file() or 'TIME: QV calculation:' not in (work/'Inspector.log').read_text():
        raise ValueError('Inspector did not reach its final QV calculation')
    valid=work/'base_error_workspace/validbase'
    if not valid.is_file() or len(valid.read_text().splitlines())!=len(fasta_dict(genome)):
        raise ValueError('Inspector base-error worker did not finish for every chromosome')
    return {'status':'COMPLETED','summary':text,'summary_sha256':checksum(summary),
            'interpretation':'Read-based error evidence; no independent polyploid phasing certification'}


def assess(args):
    meta=load(args.metadata);options=meta['settings']['assessment'];polished=load(args.polishing)
    base=Path(args.polishing).parent;genome=(base/'genome.fa').resolve()
    if checksum(genome)!=polished['sets']['genome']['sha256']:
        raise ValueError('Polished FASTA changed after validation')
    copies={label:(base/(label+'.fa')).resolve() for label in polished['sets'] if label!='genome'}
    combined=fasta_dict(genome)
    for label,path in copies.items():
        if checksum(path)!=polished['sets'][label]['sha256']:raise ValueError('Polished copy FASTA changed')
        observed=fasta_dict(path)
        expected={n:s for n,s in combined.items() if polished['groups'][n]==label}
        if observed!=expected:raise ValueError('Polished copy FASTA disagrees with combined sequence')
    dest=Path(args.out_dir).resolve();dest.parent.mkdir(parents=True,exist_ok=True)
    with assembly_scratch(dest) as scratch:
        stage=Path(scratch)/'result';stage.mkdir()
        native_tidk(genome,stage/'tidk',meta['settings']['finishing']['telomere_motif'])
        chromosomes=inspect_sequences(fasta_dict(genome),options,meta['settings']['finishing']['telomere_motif'])
        evaluations={'independent_phasing':'NOT_RUN','human_structure_review':'NOT_REVIEWED',
                     'contamination':'NOT_RUN','complex_repeat_copy_validation':'NOT_RUN'}
        if options['kmer']:
            wgs=[r for r in meta['readsets'].values() if r['technology']=='wgs']
            if not wgs:evaluations['merqury']={'status':'SKIPPED_MISSING_WGS'}
            else:
                prior=base/'wgs/manifest.json'
                if prior.is_file():
                    clean=[]
                    for record in load(prior)['libraries']:
                        for name,expected in record['cleaned'].items():
                            p=base/'wgs'/name
                            if checksum(p)!=expected:raise ValueError('Polishing WGS evidence changed')
                            clean.append(p)
                else:clean=prepare_wgs(meta,stage/'wgs',args.threads)
                records=run_merqury(genome,copies,clean,stage/'merqury',args.threads,options['mem_mb'])
                evaluations['merqury']={'status':'COMPLETED','sets':records,
                    'reads_used_for_polishing':polished['status']=='COMPLETED',
                    'phasing':'No parental hap-mers; QV is not independent phase validation'}
        else:evaluations['merqury']={'status':'DISABLED'}
        evaluations['inspector']=run_inspector(genome,args.hifi,stage/'inspector',args.threads) if options['read_mapping'] else {'status':'DISABLED'}
        buscopass=None
        if options['busco']:
            buscos={}
            for hap,path in copies.items():
                execute(['busco','-i',path,'-m','genome','--offline','-l',options['busco_lineage'],
                         '-o',hap,'--out_path',stage/'busco','-c',args.threads,
                         '--download_path',stage/'busco_cache'])
                results=list((stage/'busco'/hap).glob('short_summary.specific*.json'))
                if not results:results=list((stage/'busco'/hap).rglob('short_summary.json'))
                if len(results)!=1:raise ValueError('BUSCO summary missing or ambiguous')
                data=load(results[0]);buscos[hap]={'summary':data,'summary_sha256':checksum(results[0])}
                if 'Complete percentage' not in data.get('results',{}):
                    raise ValueError('BUSCO summary lacks Complete percentage')
            buscopass=all(b['summary']['results']['Complete percentage']>=options['min_busco_complete'] for b in buscos.values())
            evaluations['busco']={'status':'COMPLETED','copies':buscos,
                'lineage_dataset_cfg_sha256':checksum(Path(options['busco_lineage'])/'dataset.cfg'),
                'interpretation':'Run per copy; duplicated BUSCOs do not automatically imply redundant sequence'}
        else:evaluations['busco']={'status':'DISABLED'}
        screened=sum(r['gap_free_two_telomeres'] for r in chromosomes)
        expected=meta['settings']['ploidy']*meta['settings']['chromosomes_per_haplotype']
        qv=evaluations['merqury'];qvpass=None
        if qv['status']=='COMPLETED':
            qvpass=all(float(r['qv'])>=options['min_qv'] for r in qv['sets'].values())
        result={'sample_id':meta['sample_id'],'status':'COMPLETED','genome_sha256':checksum(genome),
                'genome':str(genome),'sets':polished['sets'],'groups':polished['groups'],
                'source_polishing_sha256':checksum(args.polishing),'polishing_status':polished['status'],
                'expected_chromosomes':expected,'observed_chromosomes':len(chromosomes),
                'gap_free_two_telomere_chromosomes':screened,'chromosomes':chromosomes,
                'evaluations':evaluations,'qv_threshold':options['min_qv'],'qv_screen_pass':qvpass,
                'busco_screen_pass':buscopass,'busco_complete_threshold':options['min_busco_complete'],
                'T2T_status':('INCOMPLETE_QUALITY' if qvpass is False or buscopass is False else
                    'CANDIDATE_NEEDS_REVIEW' if screened==expected and len(chromosomes)==expected else 'INCOMPLETE'),
                'tools':meta['tools'],'code':meta['code'],
                'limitations':['Two terminal arrays and no Ns are screening criteria, not proof of correct repeat copy number',
                               'Final Hi-C review, independent phasing, contamination and complex-repeat validation remain mandatory',
                               'No automatic certification or status-only approval is provided']}
        save(stage/'assessment.json',result)
        save(stage/'review.request.json',{'status':'PENDING','sample_id':meta['sample_id'],
             'genome_sha256':result['genome_sha256'],'assessment_sha256':checksum(stage/'assessment.json'),
             'required_evidence':['final_assembly_HiC_review','independent_phase_validation','contamination_screen',
                                  'centromere_rDNA_copy_validation','review_of_every_accepted_finishing_edit'],
             'reviewer':None,'reviewed_at':None,'instructions':'Attach evidence bound to this final FASTA hash; changing status alone does not certify T2T'})
        with (stage/'chromosomes.tsv').open('w') as out:
            out.write('sequence_id\tlength\tN_bases\tambiguous_bases\tleft_telomere\tright_telomere\tgap_free_two_telomeres\n')
            for row in chromosomes:
                out.write('\t'.join(map(str,[row['sequence_id'],row['length'],row['N_bases'],row['ambiguous_bases'],
                    bool(row['telomeres']['left']),bool(row['telomeres']['right']),row['gap_free_two_telomeres']]))+'\n')
        rebase_symlinks(stage,dest);write_checksums(stage);publish(stage,dest,scratch)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('metadata','polishing','out-dir'):p.add_argument('--'+name,required=True)
    p.add_argument('--hifi',nargs='+',required=True);p.add_argument('--threads',type=int,required=True)
    assess(p.parse_args())


if __name__=='__main__':main()
