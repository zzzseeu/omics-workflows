if not POLYPLOID:
    rule sample_contig_report:
        input:
            assembly_manifest='samples/{sample}/04_assembly/hifiasm/manifest.json',
            primary='samples/{sample}/04_assembly/hifiasm/primary.fa',
            alternate='samples/{sample}/04_assembly/hifiasm/alternate.fa',
            primary_index='samples/{sample}/04_assembly/hifiasm/primary.fa.fai',
            alternate_index='samples/{sample}/04_assembly/hifiasm/alternate.fa.fai',
            primary_gfa='samples/{sample}/04_assembly/hifiasm/primary.gfa',
            alternate_gfa='samples/{sample}/04_assembly/hifiasm/alternate.gfa',
            assessment='samples/{sample}/09_assessment/basic.json',
            metadata='config/samples/{sample}.metadata.json',
            code=script_path("write_sample_report"),
            shared=script_helpers("write_sample_report")
        output:
            manifest='samples/{sample}/10_report/manifest.json',
            md='samples/{sample}/10_report/report.md',
            html='samples/{sample}/10_report/report.html'
        threads: 1
        params:
            assembly_dir='samples/{sample}/04_assembly/hifiasm'
        resources:
            mem_mb=256
        log: 'samples/{sample}/logs/10_report/sample.log'
        shell:
            '{PYTHON:q} {input.code:q} --assembly-dir {params.assembly_dir:q} '
            '--assessment {input.assessment:q} --metadata {input.metadata:q} '
            '--manifest {output.manifest:q} --markdown {output.md:q} --html {output.html:q} '
            '> {log:q} 2>&1'


    rule chromosome_stage_report:
        input:
            delivery='samples/{sample}/06_scaffolding/delivery/manifest.json',
            artifacts=lambda w: [p for p in CHROMOSOME_ARTIFACTS if f'samples/{w.sample}/' in p],
            plot='samples/{sample}/06_scaffolding/contact_maps/contact_map.pdf',
            separate='samples/{sample}/06_scaffolding/contact_maps/separate_plots.pdf',
            assembly='samples/{sample}/04_assembly/hifiasm/manifest.json',
            metadata='config/samples/{sample}.metadata.json',
            assessment='samples/{sample}/09_assessment/basic.json',
            code=script_path("write_chromosome_report"),
            shared=script_helpers("write_chromosome_report")
        output:
            manifest='samples/{sample}/10_report/manifest.json',
            md='samples/{sample}/10_report/report.md',
            html='samples/{sample}/10_report/report.html'
        wildcard_constraints:
            sample=JOINT_PATTERN
        params:
            plots='samples/{sample}/06_scaffolding/contact_maps'
        threads: 1
        resources:
            mem_mb=1024
        log: 'samples/{sample}/logs/10_report/chromosomes.log'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --assembly {input.assembly:q} '
            '--delivery {input.delivery:q} --plots {params.plots:q} --assessment {input.assessment:q} '
            '--manifest {output.manifest:q} --markdown {output.md:q} '
            '--html {output.html:q} > {log:q} 2>&1'


    rule batch_contig_report:
        input:
            manifests=expand('samples/{sample}/10_report/manifest.json', sample=SAMPLES),
            plan='config/plan.json',
            code=script_path("write_batch_report"),
            shared=script_helpers("write_batch_report")
        output:
            json='10_reports/batch.json',
            md='10_reports/batch.md',
            html='10_reports/batch.html'
        threads: 1
        resources:
            mem_mb=256
        log: 'logs/10_reports/batch_report.log'
        shell:
            '{PYTHON:q} {input.code:q} --manifests {input.manifests:q} '
            '--plan {input.plan:q} --json {output.json:q} '
            '--markdown {output.md:q} --html {output.html:q} > {log:q} 2>&1'

else:
    rule polyploid_sample_report:
        input:
            metadata='config/samples/{sample}.metadata.json',
            delivery='samples/{sample}/06_scaffolding/delivery/manifest.json',
            assessment='samples/{sample}/09_assessment/basic.json',
            artifacts=lambda w: [p for p in CHROMOSOME_ARTIFACTS if f'samples/{w.sample}/' in p],
            code=script_path('write_polyploid_report'), shared=script_helpers('write_polyploid_report')
        output:
            manifest='samples/{sample}/10_report/manifest.json',
            md='samples/{sample}/10_report/report.md', html='samples/{sample}/10_report/report.html'
        threads: 1
        resources: mem_mb=1024
        log: 'samples/{sample}/logs/10_report/sample.log'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --delivery {input.delivery:q} '
            '--assessment {input.assessment:q} --manifest {output.manifest:q} '
            '--markdown {output.md:q} --html {output.html:q} > {log:q} 2>&1'

    rule batch_contig_report:
        input:
            manifests=expand('samples/{sample}/10_report/manifest.json', sample=SAMPLES),
            plan='config/plan.json',
            code=script_path("write_batch_report"),
            shared=script_helpers("write_batch_report")
        output:
            json='10_reports/batch.json',
            md='10_reports/batch.md',
            html='10_reports/batch.html'
        threads: 1
        resources:
            mem_mb=256
        log: 'logs/10_reports/batch_report.log'
        shell:
            '{PYTHON:q} {input.code:q} --manifests {input.manifests:q} '
            '--plan {input.plan:q} --json {output.json:q} '
            '--markdown {output.md:q} --html {output.html:q} > {log:q} 2>&1'


rule t2t_report:
    input:
        assessment='samples/{sample}/09_assessment/t2t/assessment.json',
        finishing='samples/{sample}/07_finishing/manifest.json',
        polishing='samples/{sample}/08_polishing/manifest.json',
        code=script_path('write_t2t_report'),shared=script_helpers('write_t2t_report')
    output:
        manifest='samples/{sample}/10_report/t2t/manifest.json',
        md='samples/{sample}/10_report/t2t/report.md',html='samples/{sample}/10_report/t2t/report.html'
    wildcard_constraints: sample=T2T_PATTERN
    threads: 1
    resources: mem_mb=512
    log: 'samples/{sample}/logs/10_report/t2t.log'
    shell:
        '{PYTHON:q} {input.code:q} --assessment {input.assessment:q} --finishing {input.finishing:q} '
        '--polishing {input.polishing:q} --manifest {output.manifest:q} '
        '--markdown {output.md:q} --html {output.html:q} > {log:q} 2>&1'
