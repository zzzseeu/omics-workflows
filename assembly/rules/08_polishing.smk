rule polish_chromosomes:
    input:
        finishing='samples/{sample}/07_finishing/manifest.json',
        fasta='samples/{sample}/07_finishing/genome.fa',
        hifi=lambda w: technology_files(w.sample,'hifi'),
        wgs=wgs_source_files,
        metadata='config/samples/{sample}.polishing.json',
        code=script_path('polish_chromosomes'), shared=script_helpers('polish_chromosomes')
    output:
        manifest='samples/{sample}/08_polishing/manifest.json',
        fasta='samples/{sample}/08_polishing/genome.fa',
        fai='samples/{sample}/08_polishing/genome.fa.fai'
    wildcard_constraints: sample=T2T_PATTERN
    params: out_dir='samples/{sample}/08_polishing'
    threads: lambda w: config['samples'][w.sample]['settings']['polishing']['threads']
    resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['polishing']['mem_mb']
    log: 'samples/{sample}/logs/08_polishing/polishing.log'
    benchmark: 'samples/{sample}/benchmarks/08_polishing/polishing.tsv'
    shell:
        '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --finishing {input.finishing:q} '
        '--hifi {input.hifi:q} --out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'
