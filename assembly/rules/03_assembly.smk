rule assemble_contigs:
    input:
        reads='samples/{sample}/02_read_qc',
        metadata='config/samples/{sample}.contigs.json',
        code=script_path('03_assemble_genome'),
        helpers=script_helpers('03_assemble_genome')
    output: directory('samples/{sample}/03_assembly/contigs')
    threads: lambda w: sample_setting(w, 'assembly', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'assembly', 'mem_mb')
    log: 'samples/{sample}/logs/03_contigs.log'
    shell:
        '{PYTHON:q} {input.code:q} --stage contigs --reads-manifest {input.reads:q}/manifest.json '
        '--metadata {input.metadata:q} --threads {threads} --out-dir {output:q} > {log:q} 2>&1'


rule scaffold_chromosomes:
    input:
        reads='samples/{sample}/02_read_qc',
        contigs='samples/{sample}/03_assembly/contigs',
        metadata='config/samples/{sample}.scaffolding.json',
        code=script_path('03_assemble_genome'),
        helpers=script_helpers('03_assemble_genome')
    output: directory('samples/{sample}/03_assembly/delivery')
    threads: lambda w: sample_setting(w, 'scaffolding', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'scaffolding', 'mem_mb')
    log: 'samples/{sample}/logs/03_scaffolding.log'
    shell:
        '{PYTHON:q} {input.code:q} --stage chromosomes --reads-manifest {input.reads:q}/manifest.json '
        '--contigs-manifest {input.contigs:q}/manifest.json --metadata {input.metadata:q} '
        '--threads {threads} --out-dir {output:q} > {log:q} 2>&1'
