rule assess_assembly:
    input:
        assembly=final_assembly,
        reads='samples/{sample}/02_read_qc',
        metadata='config/samples/{sample}.assessment.json',
        code=script_path('06_assess_assembly'),
        helpers=script_helpers('06_assess_assembly'),
        lineage=lambda w: [sample_setting(w, 'assessment', 'busco_lineage')] if
            sample_setting(w, 'assessment', 'busco_lineage') else []
    output: directory('samples/{sample}/06_assessment')
    params:
        endpoint=lambda w: '' if w.sample in CHROMOSOME_SAMPLES else '--contigs-only',
        lineage=lambda w: sample_setting(w, 'assessment', 'busco_lineage') or 'MISSING_LINEAGE',
        motif=lambda w: sample_setting(w, 'finishing', 'telomere_motif') or 'MISSING_MOTIF',
        kmer=lambda w: sample_setting(w, 'assessment', 'kmer_size'),
        mapq=lambda w: sample_setting(w, 'assessment', 'min_mapq'),
        window=lambda w: sample_setting(w, 'assessment', 'coverage_window'),
        minimum=lambda w: sample_setting(w, 'assessment', 'quast_min_contig'),
        telomere_window=lambda w: sample_setting(w, 'assessment', 'terminal_window'),
        telomere_offset=lambda w: sample_setting(w, 'assessment', 'terminal_max_distance'),
        telomere_min=lambda w: sample_setting(w, 'assessment', 'telomere_min_length'),
        retain=lambda w: '--keep-intermediates' if sample_setting(w, 'assessment', 'keep_intermediates') else ''
    threads: lambda w: sample_setting(w, 'assessment', 'threads')
    resources: mem_mb=lambda w: sample_setting(w, 'assessment', 'mem_mb')
    log: 'samples/{sample}/logs/06_assessment.log'
    shell:
        '{PYTHON:q} {input.code:q} --assembly-manifest {input.assembly:q}/manifest.json '
        '--reads-manifest {input.reads:q}/manifest.json --metadata {input.metadata:q} '
        '--busco-lineage {params.lineage:q} --telomere-motif {params.motif:q} '
        '--threads {threads} --mem-mb {resources.mem_mb} --kmer-size {params.kmer} '
        '--min-mapq {params.mapq} --coverage-window {params.window} --quast-min-contig {params.minimum} '
        '--terminal-window {params.telomere_window} --terminal-max-distance {params.telomere_offset} '
        '--telomere-min-length {params.telomere_min} {params.endpoint} {params.retain} '
        '--out-dir {output:q} > {log:q} 2>&1'
