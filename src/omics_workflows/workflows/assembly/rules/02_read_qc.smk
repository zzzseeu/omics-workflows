rule quality_reads:
    input:
        prepared='samples/{sample}/01_inputs',
        metadata='config/samples/{sample}.read_qc.json',
        code=script_path('02_quality_reads'),
        helpers=script_helpers('02_quality_reads')
    output: directory('samples/{sample}/02_read_qc')
    threads: lambda w: sample_setting(w, 'read_qc', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'read_qc', 'mem_mb')
    log: 'samples/{sample}/logs/02_read_qc.log'
    shell:
        '{PYTHON:q} {input.code:q} --inputs-manifest {input.prepared:q}/manifest.json '
        '--metadata {input.metadata:q} --threads {threads} --out-dir {output:q} > {log:q} 2>&1'
