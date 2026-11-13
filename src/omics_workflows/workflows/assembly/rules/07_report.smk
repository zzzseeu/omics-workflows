rule sample_report:
    input:
        assembly=final_assembly,
        assessment='samples/{sample}/06_assessment',
        metadata='config/samples/{sample}.report.json',
        code=script_path('07_write_assembly_report'),
        helpers=script_helpers('07_write_assembly_report')
    output: directory('samples/{sample}/07_report')
    threads: 1
    resources: mem_mb=512
    log: 'samples/{sample}/logs/07_report.log'
    shell:
        '{PYTHON:q} {input.code:q} --assembly-manifest {input.assembly:q}/manifest.json '
        '--assessment {input.assessment:q}/assessment.json --out-dir {output:q} > {log:q} 2>&1'


rule batch_report:
    input:
        reports=expand('samples/{sample}/07_report', sample=SAMPLES),
        code=script_path('07_write_assembly_report'),
        helpers=script_helpers('07_write_assembly_report')
    output:
        html='07_report/summary.html',
        table='07_report/summary.tsv'
    params:
        manifests=expand('samples/{sample}/07_report/manifest.json', sample=SAMPLES)
    threads: 1
    resources: mem_mb=512
    log: 'logs/07_batch_report.log'
    shell:
        '{PYTHON:q} {input.code:q} --reports {params.manifests:q} --out-dir 07_report > {log:q} 2>&1'
