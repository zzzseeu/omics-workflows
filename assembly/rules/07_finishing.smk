rule finish_chromosomes:
    input:
        delivery=finishing_delivery,
        sources=finishing_source_files,
        hifi=lambda w: technology_files(w.sample,'hifi'),
        ul=lambda w: technology_files(w.sample,'ont'),
        metadata='config/samples/{sample}.finishing.json',
        code=script_path('finish_chromosomes'), shared=script_helpers('finish_chromosomes')
    output:
        manifest='samples/{sample}/07_finishing/manifest.json',
        fasta='samples/{sample}/07_finishing/genome.fa',
        fai='samples/{sample}/07_finishing/genome.fa.fai',
        edits='samples/{sample}/07_finishing/edits.tsv'
    wildcard_constraints: sample=T2T_PATTERN
    params: out_dir='samples/{sample}/07_finishing'
    threads: lambda w: config['samples'][w.sample]['settings']['finishing']['threads']
    resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['finishing']['mem_mb']
    log: 'samples/{sample}/logs/07_finishing/finishing.log'
    benchmark: 'samples/{sample}/benchmarks/07_finishing/finishing.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --delivery {input.delivery:q} '
        '--hifi {input.hifi:q} --ul {input.ul:q} --out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'
