if not POLYPLOID:
    rule haphic_scaffolding:
        input:
            fasta='samples/{sample}/04_assembly/hifiasm/combined.fa',
            gfas=expand('samples/{{sample}}/04_assembly/hifiasm/{hap}.gfa', hap=['hap1','hap2']),
            bam='samples/{sample}/05_linkage/HiC.filtered.bam',
            linkage='samples/{sample}/05_linkage/manifest.json',
            metadata='config/samples/{sample}.scaffolding.json',
            code=script_path("scaffold_haplotypes"),
            shared=script_helpers("scaffold_haplotypes")
        output:
            manifest='samples/{sample}/06_scaffolding/haphic/manifest.json',
            fasta='samples/{sample}/06_scaffolding/haphic/04.build/scaffolds.fa',
            raw_agp='samples/{sample}/06_scaffolding/haphic/04.build/scaffolds.raw.agp',
            agp='samples/{sample}/06_scaffolding/haphic/04.build/scaffolds.agp'
        params:
            out_dir='samples/{sample}/06_scaffolding/haphic'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources:
            mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/06_scaffolding/haphic.log'
        benchmark: 'samples/{sample}/benchmarks/06_scaffolding/haphic.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --fasta {input.fasta:q} '
            '--bam {input.bam:q} --gfas {input.gfas:q} --linkage {input.linkage:q} '
            '--out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'


    rule chromosome_delivery:
        input:
            metadata='config/samples/{sample}.delivery.json',
            assembly='samples/{sample}/04_assembly/hifiasm/manifest.json',
            fasta='samples/{sample}/04_assembly/hifiasm/combined.fa',
            mapping='samples/{sample}/04_assembly/hifiasm/contig_haplotypes.tsv',
            haphic='samples/{sample}/06_scaffolding/haphic/manifest.json',
            raw_agp='samples/{sample}/06_scaffolding/haphic/04.build/scaffolds.raw.agp',
            native='samples/{sample}/06_scaffolding/haphic/04.build/scaffolds.fa',
            reviews=lambda w: [p for k,p in config['samples'][w.sample]['settings']['scaffolding'].items()
                               if k in ('review_agp','review_manifest') and p],
            code=script_path("export_chromosome_fasta"),
            shared=script_helpers("export_chromosome_fasta")
        output:
            expand('samples/{{sample}}/06_scaffolding/delivery/{name}', name=[
                'manifest.json','chromosomes.fa','chromosomes.fa.fai','hap1.chromosomes.fa','hap2.chromosomes.fa',
                'hap1.chromosomes.fa.fai','hap2.chromosomes.fa.fai','unplaced.fa','unplaced.fa.fai',
                'chromosomes.chrom.sizes','hap1.chromosomes.chrom.sizes','hap2.chromosomes.chrom.sizes',
                'unplaced.chrom.sizes','assembly.agp','chromosomes.agp','chromosome_names.tsv',
                'SHA256SUMS','review.request.json'])
        params:
            out_dir='samples/{sample}/06_scaffolding/delivery',
            haphic_dir='samples/{sample}/06_scaffolding/haphic'
        threads: 1
        resources:
            mem_mb=16000
        log: 'samples/{sample}/logs/06_scaffolding/delivery.log'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --assembly {input.assembly:q} '
            '--fasta {input.fasta:q} --mapping {input.mapping:q} --haphic {params.haphic_dir:q} '
            '--out-dir {params.out_dir:q} > {log:q} 2>&1'


    rule chromosome_contact_maps:
        input:
            metadata='config/samples/{sample}.plot.json',
            agp='samples/{sample}/06_scaffolding/delivery/assembly.agp',
            bam='samples/{sample}/05_linkage/HiC.filtered.bam',
            code=script_path("plot_hic_contacts"),
            shared=script_helpers("plot_hic_contacts")
        output:
            pdf='samples/{sample}/06_scaffolding/contact_maps/contact_map.pdf',
            separate='samples/{sample}/06_scaffolding/contact_maps/separate_plots.pdf'
        params:
            out_dir='samples/{sample}/06_scaffolding/contact_maps'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources:
            mem_mb=32000
        log: 'samples/{sample}/logs/06_scaffolding/contact_maps.log'
        benchmark: 'samples/{sample}/benchmarks/06_scaffolding/contact_maps.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --agp {input.agp:q} '
            '--bam {input.bam:q} --out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'

else:
    rule gphase_scaffolding:
        input:
            fasta='samples/{sample}/04_assembly/hifiasm/unitigs.fa',
            gfa='samples/{sample}/04_assembly/hifiasm/unitigs.gfa',
            index='samples/{sample}/04_assembly/hifiasm/unitigs.fa.fai',
            assembly='samples/{sample}/04_assembly/hifiasm/manifest.json',
            bam='samples/{sample}/05_linkage/HiC.filtered.bam',
            linkage='samples/{sample}/05_linkage/manifest.json',
            copies='samples/{sample}/05_linkage/copy_number/popcnv/06.genes.round.cn',
            copy_manifest='samples/{sample}/05_linkage/copy_number/manifest.json',
            metadata='config/samples/{sample}.poly_scaffolding.json',
            code=script_path('scaffold_polyploid_gphase'), shared=script_helpers('scaffold_polyploid_gphase')
        output:
            manifest='samples/{sample}/06_scaffolding/gphase/manifest.json',
            files=expand('samples/{{sample}}/06_scaffolding/gphase/gphase_output/scaffold_hap/{name}', name=GPHASE_FILES)
        params: out_dir='samples/{sample}/06_scaffolding/gphase'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/06_scaffolding/gphase.log'
        benchmark: 'samples/{sample}/benchmarks/06_scaffolding/gphase.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --assembly {input.assembly:q} '
            '--fasta {input.fasta:q} --gfa {input.gfa:q} --copies {input.copies:q} --bam {input.bam:q} '
            '--out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'

    rule polyploid_chromosome_delivery:
        input:
            assembly='samples/{sample}/04_assembly/hifiasm/manifest.json',
            raw='samples/{sample}/04_assembly/hifiasm/unitigs.fa',
            gphase='samples/{sample}/06_scaffolding/gphase/manifest.json',
            files=expand('samples/{{sample}}/06_scaffolding/gphase/gphase_output/scaffold_hap/{name}',name=GPHASE_FILES),
            metadata='config/samples/{sample}.poly_delivery.json',
            code=script_path('export_polyploid_fasta'), shared=script_helpers('export_polyploid_fasta')
        output: expand('samples/{{sample}}/06_scaffolding/delivery/{name}',name=POLY_DELIVERY_FILES)
        params: out_dir='samples/{sample}/06_scaffolding/delivery', gphase_dir='samples/{sample}/06_scaffolding/gphase'
        threads: 1
        resources: mem_mb=32000
        log: 'samples/{sample}/logs/06_scaffolding/delivery.log'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --assembly {input.assembly:q} '
            '--raw-unitigs {input.raw:q} --gphase {params.gphase_dir:q} --out-dir {params.out_dir:q} > {log:q} 2>&1'

    rule polyploid_contact_maps:
        input:
            metadata='config/samples/{sample}.poly_plot.json',
            fasta='samples/{sample}/06_scaffolding/delivery/components.fa',
            agp='samples/{sample}/06_scaffolding/delivery/assembly.agp',
            paired=lambda w: paired_files(w.sample),
            code=script_path('plot_polyploid_contacts'), shared=script_helpers('plot_polyploid_contacts')
        output:
            pdf='samples/{sample}/06_scaffolding/contact_maps/contact_map.pdf',
            separate='samples/{sample}/06_scaffolding/contact_maps/separate_plots.pdf'
        params: out_dir='samples/{sample}/06_scaffolding/contact_maps'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/06_scaffolding/contact_maps.log'
        benchmark: 'samples/{sample}/benchmarks/06_scaffolding/contact_maps.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --fasta {input.fasta:q} --agp {input.agp:q} '
            '--paired {input.paired:q} --threads {threads} --out-dir {params.out_dir:q} > {log:q} 2>&1'
