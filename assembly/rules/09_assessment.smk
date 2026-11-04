if not POLYPLOID:
    rule assess_contigs:
        input:
            manifest='samples/{sample}/04_assembly/hifiasm/manifest.json',
            primary='samples/{sample}/04_assembly/hifiasm/primary.fa',
            alternate='samples/{sample}/04_assembly/hifiasm/alternate.fa',
            qcs=lambda w: sample_files(w.sample, 'samples/{sample}/02_read_qc/{readset}.json'),
            code=script_path("collect_contig_metrics"),
            shared=script_helpers("collect_contig_metrics")
        output:
            json='samples/{sample}/09_assessment/basic.json',
            table='samples/{sample}/09_assessment/basic.tsv'
        threads: 1
        params:
            assembly_dir='samples/{sample}/04_assembly/hifiasm'
        resources:
            mem_mb=512
        log: 'samples/{sample}/logs/09_assessment/basic.log'
        benchmark: 'samples/{sample}/benchmarks/09_assessment/basic.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --assembly-dir {params.assembly_dir:q} '
            '--qcs {input.qcs:q} --json {output.json:q} --table {output.table:q} '
            '> {log:q} 2>&1'


    rule assess_chromosomes:
        input:
            delivery='samples/{sample}/06_scaffolding/delivery/manifest.json',
            artifacts=lambda w: [p for p in CHROMOSOME_ARTIFACTS if f'samples/{w.sample}/' in p],
            plot='samples/{sample}/06_scaffolding/contact_maps/contact_map.pdf',
            separate='samples/{sample}/06_scaffolding/contact_maps/separate_plots.pdf',
            code=script_path('collect_chromosome_metrics'),
            shared=script_helpers('collect_chromosome_metrics')
        output:
            assessment='samples/{sample}/09_assessment/basic.json',
            table='samples/{sample}/09_assessment/basic.tsv'
        wildcard_constraints:
            sample=JOINT_PATTERN
        params:
            plots='samples/{sample}/06_scaffolding/contact_maps'
        threads: 1
        resources:
            mem_mb=1024
        log: 'samples/{sample}/logs/09_assessment/chromosomes.log'
        shell:
            '{PYTHON:q} {input.code:q} --delivery {input.delivery:q} --plots {params.plots:q} '
            '--assessment {output.assessment:q} --table {output.table:q} > {log:q} 2>&1'

else:
    rule polyploid_assessment:
        input:
            delivery='samples/{sample}/06_scaffolding/delivery/manifest.json',
            artifacts=lambda w: [p for p in CHROMOSOME_ARTIFACTS if f'samples/{w.sample}/' in p],
            pdf='samples/{sample}/06_scaffolding/contact_maps/contact_map.pdf',
            separate='samples/{sample}/06_scaffolding/contact_maps/separate_plots.pdf',
            code=script_path('collect_polyploid_metrics'), shared=script_helpers('collect_polyploid_metrics')
        output: json='samples/{sample}/09_assessment/basic.json', table='samples/{sample}/09_assessment/basic.tsv'
        params: plots='samples/{sample}/06_scaffolding/contact_maps'
        threads: 1
        resources: mem_mb=1024
        log: 'samples/{sample}/logs/09_assessment/basic.log'
        shell:
            '{PYTHON:q} {input.code:q} --delivery {input.delivery:q} --plots {params.plots:q} '
            '--json {output.json:q} --table {output.table:q} > {log:q} 2>&1'


rule assess_t2t:
    input:
        polishing='samples/{sample}/08_polishing/manifest.json',
        fasta='samples/{sample}/08_polishing/genome.fa',
        hifi=lambda w: technology_files(w.sample,'hifi'),
        wgs=wgs_source_files,
        metadata='config/samples/{sample}.t2t_assessment.json',
        code=script_path('assess_t2t'), shared=script_helpers('assess_t2t')
    output:
        assessment='samples/{sample}/09_assessment/t2t/assessment.json',
        table='samples/{sample}/09_assessment/t2t/chromosomes.tsv',
        review='samples/{sample}/09_assessment/t2t/review.request.json'
    wildcard_constraints: sample=T2T_PATTERN
    params: out_dir='samples/{sample}/09_assessment/t2t'
    threads: lambda w: config['samples'][w.sample]['settings']['assessment']['threads']
    resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['assessment']['mem_mb']
    log: 'samples/{sample}/logs/09_assessment/t2t.log'
    benchmark: 'samples/{sample}/benchmarks/09_assessment/t2t.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --polishing {input.polishing:q} '
        '--hifi {input.hifi:q} --out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'
