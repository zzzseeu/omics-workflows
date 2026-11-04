if not POLYPLOID:
    def hic_references(sample):
        mode = config['samples'][sample]['settings']['scaffolding'].get('hic_reference', 'combined')
        return ['hap1', 'hap2'] if mode == 'per_haplotype' else ['combined']

    rule hic_bwa_index:
        input:
            fasta='samples/{sample}/04_assembly/hifiasm/{reference}.fa'
        output:
            expand('samples/{{sample}}/05_linkage/bwa_index/{{reference}}.{suffix}', suffix=['amb','ann','bwt','pac','sa'])
        wildcard_constraints:
            reference='combined|hap1|hap2'
        params:
            prefix='samples/{sample}/05_linkage/bwa_index/{reference}',
            bwa=lambda w: config['tools']['bwa']['path']
        threads: 1
        resources:
            mem_mb=12000
        log: 'samples/{sample}/logs/05_linkage/bwa_index.{reference}.log'
        shell:
            '{params.bwa:q} index -p {params.prefix:q} {input.fasta:q} > {log:q} 2>&1'


    rule hic_library_mapping:
        input:
            indexes=lambda w: [
                f'samples/{w.sample}/05_linkage/bwa_index/{ref}.{suffix}'
                for ref in hic_references(w.sample)
                for suffix in ['amb','ann','bwt','pac','sa']],
            paired='samples/{sample}/01_inputs/{readset}.paired.json',
            metadata='config/samples/{sample}.linkage.json',
            code=script_path("map_filter_hic"),
            shared=script_helpers("map_filter_hic")
        output:
            bam='samples/{sample}/05_linkage/libraries/{readset}.filtered.bam',
            metrics='samples/{sample}/05_linkage/libraries/{readset}.json'
        params:
            prefix=lambda w: [f'samples/{w.sample}/05_linkage/bwa_index/{ref}'
                for ref in hic_references(w.sample)]
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources:
            mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/05_linkage/{readset}.log'
        benchmark: 'samples/{sample}/benchmarks/05_linkage/{readset}.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --paired {input.paired:q} '
            '--prefix {params.prefix:q} --bam {output.bam:q} --out {output.metrics:q} '
            '--threads {threads} > {log:q} 2>&1'


    rule merge_hic_libraries:
        input:
            bams=lambda w: [f'samples/{w.sample}/05_linkage/libraries/{r}.filtered.bam' for r in sample_readsets(w.sample, 'hic')],
            metrics=lambda w: [f'samples/{w.sample}/05_linkage/libraries/{r}.json' for r in sample_readsets(w.sample, 'hic')],
            code=script_path("merge_validate_hic_bams"),
            shared=script_helpers("merge_validate_hic_bams")
        output:
            bam='samples/{sample}/05_linkage/HiC.filtered.bam',
            metrics='samples/{sample}/05_linkage/manifest.json'
        threads: 1
        resources:
            mem_mb=1024
        log: 'samples/{sample}/logs/05_linkage/merge.log'
        benchmark: 'samples/{sample}/benchmarks/05_linkage/merge.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --inputs {input.bams:q} --metrics {input.metrics:q} '
            '--bam {output.bam:q} --out {output.metrics:q} > {log:q} 2>&1'

else:
    rule polyploid_bwa_index:
        input: fasta='samples/{sample}/04_assembly/hifiasm/unitigs.fa'
        output: expand('samples/{{sample}}/05_linkage/bwa_index/unitigs.{suffix}',suffix=['amb','ann','bwt','pac','sa'])
        params: prefix='samples/{sample}/05_linkage/bwa_index/unitigs', tool=config['tools']['bwa']['path']
        threads: 1
        resources: mem_mb=16000
        log: 'samples/{sample}/logs/05_linkage/bwa_index.log'
        shell: '{params.tool:q} index -p {params.prefix:q} {input.fasta:q} > {log:q} 2>&1'

    rule polyploid_hic_mapping:
        input:
            indexes=expand('samples/{{sample}}/05_linkage/bwa_index/unitigs.{suffix}',suffix=['amb','ann','bwt','pac','sa']),
            paired='samples/{sample}/01_inputs/{readset}.paired.json',
            code=script_path('map_polyploid_hic'), shared=script_helpers('map_polyploid_hic')
        output:
            bam='samples/{sample}/05_linkage/libraries/{readset}.filtered.bam',
            metrics='samples/{sample}/05_linkage/libraries/{readset}.json'
        params: prefix='samples/{sample}/05_linkage/bwa_index/unitigs'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/05_linkage/{readset}.log'
        benchmark: 'samples/{sample}/benchmarks/05_linkage/{readset}.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --paired {input.paired:q} --prefix {params.prefix:q} '
            '--bam {output.bam:q} --out {output.metrics:q} --threads {threads} > {log:q} 2>&1'

    rule polyploid_merge_hic:
        input:
            bams=lambda w: [f'samples/{w.sample}/05_linkage/libraries/{r}.filtered.bam' for r in sample_readsets(w.sample,'hic')],
            metrics=lambda w: [f'samples/{w.sample}/05_linkage/libraries/{r}.json' for r in sample_readsets(w.sample,'hic')],
            code=script_path('merge_validate_hic_bams'), shared=script_helpers('merge_validate_hic_bams')
        output: bam='samples/{sample}/05_linkage/HiC.filtered.bam', metrics='samples/{sample}/05_linkage/manifest.json'
        threads: 1
        resources: mem_mb=2048
        log: 'samples/{sample}/logs/05_linkage/merge.log'
        shell:
            '{PYTHON:q} {input.code:q} --inputs {input.bams:q} --metrics {input.metrics:q} '
            '--bam {output.bam:q} --out {output.metrics:q} > {log:q} 2>&1'

    rule polyploid_unitig_copy_number:
        input:
            fasta='samples/{sample}/04_assembly/hifiasm/unitigs.fa',
            index='samples/{sample}/04_assembly/hifiasm/unitigs.fa.fai',
            hifi=lambda w: technology_files(w.sample,'hifi'),
            code=script_path('estimate_unitig_copies'), shared=script_helpers('estimate_unitig_copies')
        output:
            manifest='samples/{sample}/05_linkage/copy_number/manifest.json',
            copies='samples/{sample}/05_linkage/copy_number/popcnv/06.genes.round.cn'
        params: out_dir='samples/{sample}/05_linkage/copy_number'
        threads: lambda w: config['samples'][w.sample]['settings']['scaffolding']['threads']
        resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['scaffolding']['mem_mb']
        log: 'samples/{sample}/logs/05_linkage/copy_number.log'
        benchmark: 'samples/{sample}/benchmarks/05_linkage/copy_number.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --fasta {input.fasta:q} --hifi {input.hifi:q} '
            '--out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'
