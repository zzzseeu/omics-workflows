if not POLYPLOID:
    rule hifiasm_contigs:
        input:
            fastqs=lambda w: sample_files(w.sample, 'samples/{sample}/01_inputs/{readset}.fastq.gz'),
            manifests=lambda w: sample_files(w.sample, 'samples/{sample}/01_inputs/{readset}.json'),
            qcs=lambda w: sample_files(w.sample, 'samples/{sample}/02_read_qc/{readset}.json'),
            spec='config/samples/{sample}.assembly.json',
            code=script_path("assemble_hifi_contigs"),
            shared=script_helpers("assemble_hifi_contigs")
        output:
            manifest='samples/{sample}/04_assembly/hifiasm/manifest.json',
            primary='samples/{sample}/04_assembly/hifiasm/primary.fa',
            alternate='samples/{sample}/04_assembly/hifiasm/alternate.fa',
            primary_gfa='samples/{sample}/04_assembly/hifiasm/primary.gfa',
            alternate_gfa='samples/{sample}/04_assembly/hifiasm/alternate.gfa'
        params:
            out_dir='samples/{sample}/04_assembly/hifiasm'
        threads: lambda w: config['samples'][w.sample]['settings']['assembly']['threads']
        resources:
            mem_mb=lambda w: config['samples'][w.sample]['settings']['assembly']['mem_mb']
        log: 'samples/{sample}/logs/04_assembly/hifiasm.log'
        benchmark: 'samples/{sample}/benchmarks/04_assembly/hifiasm.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --spec {input.spec:q} '
            '--fastqs {input.fastqs:q} --manifests {input.manifests:q} '
            '--qcs {input.qcs:q} --out-dir {params.out_dir:q} '
            '--threads {threads} > {log:q} 2>&1'


    rule joint_hifiasm:
        input:
            hifi=lambda w: technology_files(w.sample, 'hifi'),
            ul=lambda w: technology_files(w.sample, 'ont'),
            paired=lambda w: paired_files(w.sample),
            manifests=lambda w: sample_files(w.sample, 'samples/{sample}/01_inputs/{readset}.json'),
            qcs=lambda w: sample_files(w.sample, 'samples/{sample}/02_read_qc/{readset}.json'),
            metadata='config/samples/{sample}.joint.json',
            code=script_path("assemble_phased_contigs"),
            shared=script_helpers("assemble_phased_contigs")
        output:
            manifest='samples/{sample}/04_assembly/hifiasm/manifest.json',
            hap1='samples/{sample}/04_assembly/hifiasm/hap1.fa',
            hap2='samples/{sample}/04_assembly/hifiasm/hap2.fa',
            gfa1='samples/{sample}/04_assembly/hifiasm/hap1.gfa',
            gfa2='samples/{sample}/04_assembly/hifiasm/hap2.gfa',
            combined='samples/{sample}/04_assembly/hifiasm/combined.fa',
            mapping='samples/{sample}/04_assembly/hifiasm/contig_haplotypes.tsv'
        wildcard_constraints:
            sample=JOINT_PATTERN
        params:
            out_dir='samples/{sample}/04_assembly/hifiasm'
        threads: lambda w: config['samples'][w.sample]['settings']['assembly']['threads']
        resources:
            mem_mb=lambda w: config['samples'][w.sample]['settings']['assembly']['mem_mb']
        log: 'samples/{sample}/logs/04_assembly/hifiasm.log'
        benchmark: 'samples/{sample}/benchmarks/04_assembly/hifiasm.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} '
            '--hifi {input.hifi:q} --ul {input.ul:q} --paired {input.paired:q} '
            '--qcs {input.qcs:q} --manifests {input.manifests:q} '
            '--out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'


    rule index_contig_set:
        input:
            fasta='samples/{sample}/04_assembly/hifiasm/{set_id}.fa'
        output:
            index='samples/{sample}/04_assembly/hifiasm/{set_id}.fa.fai'
        wildcard_constraints:
            set_id='primary|alternate|hap1|hap2|combined'
        params:
            samtools=config['tools']['samtools']['path']
        threads: 1
        resources:
            mem_mb=256
        log: 'samples/{sample}/logs/04_assembly/index/{set_id}.log'
        shell:
            '(if [ -s {input.fasta:q} ]; then {params.samtools:q} faidx {input.fasta:q}; '
            'else : > {output.index:q}; fi) > {log:q} 2>&1'

else:
    rule polyploid_unitigs:
        input:
            hifi=lambda w: technology_files(w.sample, 'hifi'),
            ul=lambda w: technology_files(w.sample, 'ont'),
            qcs=lambda w: sample_files(w.sample, 'samples/{sample}/02_read_qc/{readset}.json'),
            manifests=lambda w: sample_files(w.sample, 'samples/{sample}/01_inputs/{readset}.json'),
            audit='samples/{sample}/01_inputs/long_read_ids.json',
            metadata='config/samples/{sample}.poly_assembly.json',
            code=script_path('assemble_polyploid_unitigs'), shared=script_helpers('assemble_polyploid_unitigs')
        output:
            manifest='samples/{sample}/04_assembly/hifiasm/manifest.json',
            fasta='samples/{sample}/04_assembly/hifiasm/unitigs.fa',
            gfa='samples/{sample}/04_assembly/hifiasm/unitigs.gfa'
        params: out_dir='samples/{sample}/04_assembly/hifiasm'
        threads: lambda w: config['samples'][w.sample]['settings']['assembly']['threads']
        resources: mem_mb=lambda w: config['samples'][w.sample]['settings']['assembly']['mem_mb']
        log: 'samples/{sample}/logs/04_assembly/hifiasm.log'
        benchmark: 'samples/{sample}/benchmarks/04_assembly/hifiasm.tsv'
        shell:
            '{PYTHON:q} {input.code:q} --metadata {input.metadata:q} --hifi {input.hifi:q} --ul {input.ul:q} '
            '--qcs {input.qcs:q} --manifests {input.manifests:q} --audit {input.audit:q} '
            '--out-dir {params.out_dir:q} --threads {threads} > {log:q} 2>&1'

    rule index_polyploid_unitigs:
        input: 'samples/{sample}/04_assembly/hifiasm/unitigs.fa'
        output: 'samples/{sample}/04_assembly/hifiasm/unitigs.fa.fai'
        params: tool=config['tools']['samtools']['path']
        threads: 1
        resources: mem_mb=1024
        log: 'samples/{sample}/logs/04_assembly/index.log'
        shell: '{params.tool:q} faidx {input:q} > {log:q} 2>&1'
