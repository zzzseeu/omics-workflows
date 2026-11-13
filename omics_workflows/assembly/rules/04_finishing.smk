rule finish_chromosomes:
    input:
        assembly=stage_input,
        sources=external_sources,
        reads='samples/{sample}/02_read_qc',
        metadata='config/samples/{sample}.finishing.json',
        code=script_path('04_finish_chromosomes'),
        helpers=script_helpers('04_finish_chromosomes')
    output: directory('samples/{sample}/04_finishing')
    params:
        source=lambda w: config['samples'][w.sample]['settings']['finishing']['source_delivery'] or
            f'samples/{w.sample}/03_assembly/delivery/manifest.json'
    threads: lambda w: sample_setting(w, 'finishing', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'finishing', 'mem_mb')
    log: 'samples/{sample}/logs/04_finishing.log'
    shell:
        '{PYTHON:q} {input.code:q} --assembly-manifest {params.source:q} '
        '--reads-manifest {input.reads:q}/manifest.json --metadata {input.metadata:q} '
        '--threads {threads} --out-dir {output:q} > {log:q} 2>&1'
