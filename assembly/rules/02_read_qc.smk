rule hifi_read_qc:
    input:
        fastq='samples/{sample}/01_inputs/{readset}.fastq.gz',
        manifest='samples/{sample}/01_inputs/{readset}.json',
        prior_qc=lambda w: [config['samples'][w.sample]['readsets'][w.readset]['prepared_reuse']['identity']['qc']]
                           if config['samples'][w.sample]['readsets'][w.readset].get('prepared_reuse') else [],
        code=script_path("check_read_quality"),
        shared=script_helpers("check_read_quality")
    output:
        stats='samples/{sample}/02_read_qc/{readset}.tsv',
        json='samples/{sample}/02_read_qc/{readset}.json'
    threads: min(2, config['budget']['cores'])
    resources:
        mem_mb=1024
    params:
        seqkit=config['tools']['seqkit']['path']
    log: 'samples/{sample}/logs/02_read_qc/{readset}.log'
    benchmark: 'samples/{sample}/benchmarks/02_read_qc/{readset}.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --fastq {input.fastq:q} --prepared-manifest {input.manifest:q} '
        '--stats {output.stats:q} --json {output.json:q} --seqkit {params.seqkit:q} --threads {threads} > {log:q} 2>&1'
