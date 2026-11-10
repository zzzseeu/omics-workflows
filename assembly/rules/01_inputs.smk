rule prepare_inputs:
    input:
        metadata='config/samples/{sample}.inputs.json',
        sources=raw_files,
        code=script_path('01_prepare_inputs'),
        helpers=script_helpers('01_prepare_inputs')
    output: directory('samples/{sample}/01_inputs')
    threads: min(2, config['budget']['cores'])
    resources: mem_mb=min(2048, config['budget']['mem_mb'])
    log: 'samples/{sample}/logs/01_inputs.log'
    shell:
        '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} '
        '--threads {threads} --out-dir {output:q} > {log:q} 2>&1'
