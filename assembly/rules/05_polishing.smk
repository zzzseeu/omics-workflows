rule polish_chromosomes:
    input:
        assembly='samples/{sample}/04_finishing',
        reads='samples/{sample}/02_read_qc',
        metadata='config/samples/{sample}.polishing.json',
        code=script_path('05_polish_chromosomes'),
        helpers=script_helpers('05_polish_chromosomes')
    output: directory('samples/{sample}/05_polishing')
    threads: lambda w: sample_setting(w, 'polishing', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'polishing', 'mem_mb')
    log: 'samples/{sample}/logs/05_polishing.log'
    shell:
        '{PYTHON:q} {input.code:q} --assembly-manifest {input.assembly:q}/manifest.json '
        '--reads-manifest {input.reads:q}/manifest.json --metadata {input.metadata:q} '
        '--threads {threads} --out-dir {output:q} > {log:q} 2>&1'
