rule prepare_hifi_readset:
    input:
        source=lambda w: config['samples'][w.sample]['readsets'][w.readset]['file1'],
        spec='config/readsets/{sample}/{readset}.json',
        reused=lambda w: [config['samples'][w.sample]['readsets'][w.readset]['prepared_reuse']['identity'][k]
                          for k in ('prepared_fastq','manifest','qc')]
                         if config['samples'][w.sample]['readsets'][w.readset].get('prepared_reuse') else [],
        code=lambda w: script_path('reuse_prepared_reads' if config['samples'][w.sample]['readsets'][w.readset].get('prepared_reuse') else 'prepare_long_reads'),
        shared=lambda w: script_helpers('reuse_prepared_reads' if config['samples'][w.sample]['readsets'][w.readset].get('prepared_reuse') else 'prepare_long_reads')
    output:
        fastq='samples/{sample}/01_inputs/{readset}.fastq.gz',
        manifest='samples/{sample}/01_inputs/{readset}.json'
    wildcard_constraints:
        readset=LONG_READSET_PATTERN
    threads: min(2, config['budget']['cores'])
    resources:
        mem_mb=256
    log: 'samples/{sample}/logs/01_inputs/{readset}.log'
    benchmark: 'samples/{sample}/benchmarks/01_inputs/{readset}.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --spec {input.spec:q} '
        '--fastq {output.fastq:q} --manifest {output.manifest:q} '
        '--threads {threads} > {log:q} 2>&1'


rule validate_hic_pair:
    input:
        source1=lambda w: config['samples'][w.sample]['readsets'][w.readset]['file1'],
        source2=lambda w: config['samples'][w.sample]['readsets'][w.readset]['file2'],
        spec='config/readsets/{sample}/{readset}.json',
        code=script_path("validate_hic_pairs"),
        shared=script_helpers("validate_hic_pairs")
    output:
        json='samples/{sample}/01_inputs/{readset}.paired.json'
    threads: 1
    resources:
        mem_mb=512
    log: 'samples/{sample}/logs/01_inputs/{readset}.paired.log'
    benchmark: 'samples/{sample}/benchmarks/01_inputs/{readset}.paired.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --spec {input.spec:q} --out {output.json:q} > {log:q} 2>&1'


if POLYPLOID:
    rule audit_polyploid_read_ids:
        input:
            hifi=lambda w: technology_files(w.sample,'hifi'),
            ul=lambda w: technology_files(w.sample,'ont'),
            code=script_path('audit_long_read_ids'), shared=script_helpers('audit_long_read_ids')
        output: 'samples/{sample}/01_inputs/long_read_ids.json'
        threads: 1
        resources: mem_mb=2048
        log: 'samples/{sample}/logs/01_inputs/long_read_ids.log'
        shell:
            '{PYTHON:q} {input.code:q} --hifi {input.hifi:q} --ul {input.ul:q} --out {output:q} > {log:q} 2>&1'
