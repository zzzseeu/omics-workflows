#!/usr/bin/env python3
"""Validate real GPhase coordinates plus known four-copy remapping/plots/reports."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'assembly/scripts'))
from agp_utils import reconstruct
from export_polyploid_fasta import export
from fasta_utils import fasta_dict, write_fasta
from io_utils import checksum, load, save


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--known-fixture',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    checks=[]
    def script(name,*options):return [sys.executable,ROOT/f'assembly/scripts/{name}.py',*options]
    def run(name,command):
        with (out/f'{name}.log').open('w') as log:
            subprocess.run(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
        checks.append(name);print('PASS',name,flush=True)
    # This is the upstream's real triploid reference; it exercises graph trimming
    # and collapsed-copy aliases. The public route remains restricted to ploidy 4.
    reference=args.reference.resolve();raw=reference/'A09.h.asm.bp.p_utg.gfa.fa'
    save(out/'reference_meta.json',{'sample_id':'UPSTREAM_REFERENCE','settings':{'ploidy':3,'chromosomes_per_haplotype':7}})
    save(out/'reference_assembly.json',{'assembly_id':'upstream_reference','fasta_sha256':checksum(raw)})
    run('real_gphase_delivery_coordinates',script('export_polyploid_fasta','--metadata',out/'reference_meta.json',
        '--assembly',out/'reference_assembly.json','--raw-unitigs',raw,'--gphase',reference,
        '--out-dir',out/'reference_delivery'))
    fixture=args.known_fixture.resolve();native=out/'four_copy/gphase_output/scaffold_hap';native.mkdir(parents=True)
    components=fasta_dict(fixture/'combined.fa');rows=[];chromosomes={}
    for hap in range(1,5):
        source_hap=1 if hap<=2 else 2
        indices=range(1,5) if hap in (1,3) else range(5,9)
        name=f'chr1g{hap}_scaffold_1';obj=[];position=1
        for part,i in enumerate(indices):
            seq=components[f'hap{source_hap}_ctg{i}']
            if part:
                obj.append([name,str(position),str(position+99),str(len(obj)+1),'U','100','scaffold','yes','proximity_ligation'])
                position+=100
            obj.append([name,str(position),str(position+len(seq)-1),str(len(obj)+1),'W',f'hap{source_hap}_ctg{i}','1',str(len(seq)),'+'])
            position+=len(seq)
        rows.extend(obj);chromosomes[name]=reconstruct(obj,components)
    for kind in ('unitig','contig'):
        shutil.copyfile(fixture/'combined.fa',native/f'gphase_final.{kind}.fasta')
        (native/f'gphase_final.{kind}.scaffold.agp').write_text(''.join('\t'.join(row)+'\n' for row in rows))
        with (native/f'gphase_final.{kind}.scaffold.fasta').open('w') as f:
            for name,seq in chromosomes.items():write_fasta(f,name,seq)
    (native/'gphase_final.ctg2utg.txt').write_text('known structure coordinate fixture\n')
    meta=load(fixture/'metadata.json');meta['species']='Synthetic four-copy coordinate fixture'
    meta['settings'].update(ploidy=4,chromosomes_per_haplotype=1)
    meta['settings']['scaffolding'].update(review_agp=None,review_manifest=None,threads=16)
    save(out/'metadata.json',meta)
    save(out/'assembly.json',{'assembly_id':'known_four_copy','fasta_sha256':checksum(fixture/'combined.fa')})
    run('four_copy_export',script('export_polyploid_fasta','--metadata',out/'metadata.json','--assembly',out/'assembly.json',
        '--raw-unitigs',fixture/'combined.fa','--gphase',out/'four_copy','--out-dir',out/'delivery'))
    assert load(out/'delivery/manifest.json')['chromosome_counts']=={f'hap{i}':1 for i in range(1,5)}
    run('remap_and_contact_plots',script('plot_polyploid_contacts','--metadata',out/'metadata.json',
        '--fasta',out/'delivery/components.fa','--agp',out/'delivery/assembly.agp',
        '--paired',fixture/'paired.json','--threads',16,'--out-dir',out/'plots'))
    for name in ('contact_map.pdf','separate_plots.pdf'):
        assert (out/'plots'/name).read_bytes().startswith(b'%PDF')
    run('polyploid_assessment',script('collect_polyploid_metrics','--delivery',out/'delivery/manifest.json',
        '--plots',out/'plots','--json',out/'assessment/basic.json','--table',out/'assessment/basic.tsv'))
    run('polyploid_report',script('write_polyploid_report','--metadata',out/'metadata.json',
        '--delivery',out/'delivery/manifest.json','--assessment',out/'assessment/basic.json',
        '--manifest',out/'report/manifest.json','--markdown',out/'report/report.md','--html',out/'report/report.html'))
    save(out/'validation.json',{'status':'PASS','checks':checks,
         'limits':'Real upstream triploid coordinate validation and synthetic four-copy output/remapping contracts; not rice phase or T2T certification'})
    print('POLYPLOID_DELIVERY_PLOT_REPORT=PASS',flush=True)


if __name__=='__main__':main()
