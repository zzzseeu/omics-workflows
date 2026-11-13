"""Physical-library Hi-C mapping, filtering, pair validation and contact maps."""

from pathlib import Path
from types import SimpleNamespace

from omics_workflows.assembly.lib.io import assembly_scratch as native_scratch
from omics_workflows.assembly.lib.io import (
    checksum,
    execute,
    load,
    pipeline,
    save,
    signature,
    source_check,
)


def combine_haplotype_bams(inputs, output):
    """按参考名称把两份单倍型BAM重建到联合头，并为读段ID添加参考来源前缀。

    同时转换自身及配对参考编号；同一物理分子在不同参考上的比对属于多份观察，
    不能作为互相独立的分相证据，前缀用于避免混淆这些观察。
    """
    import pysam

    headers = []
    for path in inputs:
        with pysam.AlignmentFile(str(path), "rb") as bam:
            headers.append(bam.header.to_dict())
    header = dict(headers[0])
    header["SQ"] = [sq for h in headers for sq in h["SQ"]]
    names = [sq["SN"] for sq in header["SQ"]]
    if len(names) != len(set(names)):
        raise ValueError("Per-haplotype BAM reference names must be disjoint")
    for i, h in enumerate(headers, 1):
        if any(not sq["SN"].startswith(f"hap{i}_") for sq in h["SQ"]):
            raise ValueError(f"Expected hap{i}_ contig names in per-haplotype BAM")
    header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(str(output), "wb", header=header) as target:
        for i, path in enumerate(inputs, 1):
            with pysam.AlignmentFile(str(path), "rb") as source:
                for read in source:
                    # Parsing SAM against the union header translates RNEXT as
                    # well as RNAME; concatenating these BAMs would corrupt IDs.
                    mapped = pysam.AlignedSegment.fromstring(read.to_string(), header)
                    name = f"hap{i}:{read.query_name}"
                    if len(name.encode()) > 254:
                        raise ValueError("Namespaced Hi-C query name exceeds BAM limit")
                    mapped.query_name = name
                    target.write(mapped)


def bam_pairs(path):
    """逐对检查BAM中相邻主要读段的名称与R1/R2标志，返回保留对数；空结果或配对破损时报错。"""
    import pysam

    count = 0
    with pysam.AlignmentFile(str(path), "rb") as f:
        reads = iter(f)
        for a in reads:
            b = next(reads, None)
            if b is None or a.query_name != b.query_name or not a.is_paired or not b.is_paired:
                raise ValueError("Hi-C BAM must retain adjacent primary read pairs")
            if a.is_read1 == b.is_read1 or a.is_read2 == b.is_read2:
                raise ValueError("Invalid mate flags in Hi-C BAM")
            count += 1
    if not count:
        raise ValueError("No Hi-C pairs retained after filtering")
    return count


def map_hic(args):
    """按一个物理Hi-C文库执行BWA比对、samblaster标记重复和成对MAPQ/NM过滤，写BAM及保留率。"""
    metadata, paired = load(args.metadata), load(args.paired)
    source = paired["source_spec"]
    source_check(source)
    options = metadata["settings"]["scaffolding"]
    mode = options.get("hic_reference", "combined")
    prefixes = args.prefix if isinstance(args.prefix, list) else [args.prefix]
    if mode not in ("combined", "per_haplotype") or len(prefixes) != (2 if mode == "per_haplotype" else 1):
        raise ValueError("Hi-C reference mode and BWA index prefixes disagree")
    out = Path(args.bam)
    out.parent.mkdir(parents=True, exist_ok=True)
    with native_scratch(out, prefix=".hic-") as temp:
        references, filtered_inputs = [], []
        for index, prefix in enumerate(prefixes, 1):
            label = f"hap{index}" if mode == "per_haplotype" else "combined"
            raw = Path(temp) / f"{label}.mapped.bam"
            filtered = Path(temp) / f"{label}.filtered.bam"
            pipeline(
                [
                    [
                        "bwa",
                        "mem",
                        "-5SP",
                        "-t",
                        max(1, args.threads - 3),
                        "-R",
                        f'@RG\\tID:{source["readset_id"]}\\tSM:{metadata["sample_id"]}'
                        f'\\tLB:{source["library_id"]}',
                        prefix,
                        source["file1"],
                        source["file2"],
                    ],
                    ["samblaster"],
                    ["samtools", "view", "-F", "3340", "-b", "-"],
                ],
                raw,
            )
            pipeline(
                [
                    [
                        metadata["tools"]["filter_bam"]["path"],
                        raw,
                        options["mapq"],
                        "--nm",
                        options["nm_exclusive"],
                        "--remove-singletons",
                        "--threads",
                        max(1, args.threads - 2),
                    ],
                    ["samtools", "view", "-b", "-"],
                ],
                filtered,
            )
            execute(["samtools", "quickcheck", "-v", filtered])
            count = bam_pairs(filtered)
            fraction = count / paired["num_pairs"]
            result = {"reference": label, "retained_pairs": count, "retained_fraction": fraction}
            references.append(result)
            save(
                Path(temp) / "retention.json", {"input_pairs": paired["num_pairs"], "references": references}
            )
            if fraction < options.get("min_retained_fraction", 0.001):
                raise ValueError(
                    f'Hi-C retained only {count}/{paired["num_pairs"]} pairs '
                    f"({fraction:.6%}) against {label}; below min_retained_fraction. "
                    "Inspect ambiguous MAPQ=0 alignments and reference duplication before scaffolding."
                )
            filtered_inputs.append(filtered)
        filtered = filtered_inputs[0]
        if mode == "per_haplotype":
            filtered = Path(temp) / "union.filtered.bam"
            combine_haplotype_bams(filtered_inputs, filtered)
        execute(["samtools", "quickcheck", "-v", filtered])
        pairs = bam_pairs(filtered)
        if pairs != sum(r["retained_pairs"] for r in references):
            raise ValueError("Hi-C reference observation counts were not conserved")
        source_check(source)
        filtered.replace(out)
        save(
            args.out,
            {
                "status": "PASS",
                "readset_id": source["readset_id"],
                "input_pairs": paired["num_pairs"],
                "retained_pairs": pairs,
                "retained_fraction": pairs / (paired["num_pairs"] * len(prefixes)),
                "reference_pair_opportunities": paired["num_pairs"] * len(prefixes),
                "bam_sha256": checksum(out),
                "hic_reference": mode,
                "references": references,
                "count_unit": "reference_pair_observations" if mode == "per_haplotype" else "physical_pairs",
                "shared_molecules_across_references": mode == "per_haplotype",
                "phase_validation": "NOT_ASSESSED",
                "mapping": "bwa mem -5SP; samblaster per physical library; -F3340; pairwise MAPQ/NM",
                "thresholds": {"mapq_min": options["mapq"], "NM_exclusive": options["nm_exclusive"]},
                "sort_order": "query-adjacent, not coordinate-sorted",
            },
        )


def merge_hic(args):
    # Concatenation preserves adjacent pairs. Name sorting across libraries is unnecessary.
    """合并不同物理文库的已过滤BAM，保留相邻配对关系并记录总保留对数，不跨文库重新去重。"""
    out = Path(args.bam)
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_suffix(".partial.bam")
    try:
        import pysam

        headers = []
        for path in args.inputs:
            with pysam.AlignmentFile(str(path), "rb") as bam:
                headers.append(list(zip(bam.references, bam.lengths)))
        if any(h != headers[0] for h in headers[1:]):
            raise ValueError("Hi-C library BAM reference dictionaries differ; refusing unsafe concatenation")
        execute(["samtools", "cat", "-o", temporary, *args.inputs])
        execute(["samtools", "quickcheck", "-v", temporary])
        expected = sum(load(p)["retained_pairs"] for p in args.metrics)
        actual = bam_pairs(temporary)
        if actual != expected:
            raise ValueError("Merged Hi-C pair count differs from library totals")
        temporary.replace(out)
        libraries = [load(p) for p in args.metrics]
        modes = {r.get("hic_reference", "combined") for r in libraries}
        if len(modes) != 1:
            raise ValueError("Hi-C libraries use inconsistent reference modes")
        mode = modes.pop()
        save(
            args.out,
            {
                "status": "PASS",
                "retained_pairs": actual,
                "bam_sha256": checksum(out),
                "hic_reference": mode,
                "phase_validation": "NOT_ASSESSED",
                "count_unit": "reference_pair_observations" if mode == "per_haplotype" else "physical_pairs",
                "libraries": libraries,
            },
        )
    finally:
        temporary.unlink(missing_ok=True)


def pair_records(reads, directory):
    """把已核验的Hi-C清洁读段转换为物理文库比对规格文件，返回这些规格文件路径列表。"""
    paths = []
    for rid, record in reads["readsets"].items():
        if record["technology"] != "hic":
            continue
        spec = {
            "sample_id": reads["sample_id"],
            "readset_id": rid,
            "library_id": record["library_id"],
            "file1": record["paths"][0],
            "file2": record["paths"][1],
            "source_signature": signature(record["paths"][0]),
            "source_signature_2": signature(record["paths"][1]),
        }
        path = directory / f"{rid}.paired.json"
        save(path, {"source_spec": spec, "num_pairs": record["statistics"]["num_pairs"]})
        paths.append(path)
    if not paths:
        raise ValueError("Chromosome scaffolding requires clean Hi-C pairs")
    return paths


def map_contacts(references, paired, metadata, directory, threads, raw=False):
    """对指定参考及物理文库生成合并Hi-C BAM，返回BAM路径和比对统计。

    raw模式为GPhase保留MAPQ0观察；普通模式采用约定的成对过滤，不混用两种统计口径。
    """
    directory.mkdir(parents=True, exist_ok=True)
    prefixes = []
    for index, reference in enumerate(references):
        prefix = directory / f"reference{index}"
        execute(["bwa", "index", "-p", prefix, Path(reference).resolve()])
        prefixes.append(prefix)
    bams, metrics = [], []
    for number, path in enumerate(paired):
        bam = directory / f"library{number}.bam"
        stat = directory / f"library{number}.json"
        record = load(path)
        if raw:
            spec = record["source_spec"]
            group = f'@RG\\tID:{spec["readset_id"]}\\tSM:{spec["sample_id"]}\\tLB:{spec["library_id"]}'
            pipeline(
                [
                    [
                        "bwa",
                        "mem",
                        "-5SP",
                        "-t",
                        max(1, threads - 3),
                        "-R",
                        group,
                        prefixes[0],
                        spec["file1"],
                        spec["file2"],
                    ],
                    ["samblaster"],
                    ["samtools", "view", "-F", "3340", "-b", "-"],
                ],
                bam,
            )
            save(stat, {"retained_pairs": bam_pairs(bam), "input_pairs": record["num_pairs"]})
        else:
            map_hic(
                SimpleNamespace(
                    metadata=metadata, paired=path, prefix=prefixes, bam=bam, out=stat, threads=threads
                )
            )
        bams.append(bam)
        metrics.append(stat)
    merged = directory / "hic.bam"
    merge_hic(SimpleNamespace(inputs=bams, metrics=metrics, bam=merged, out=directory / "mapping.json"))
    return merged, load(directory / "mapping.json")


def plot_contacts(reference, agp, reads, metadata, directory, threads):
    """在最终组装坐标下重做Hi-C比对和接触图，返回比对统计；保留图后清理BAM和索引。"""
    directory.mkdir(parents=True, exist_ok=True)
    paired = pair_records(reads, directory)
    bam, statistics = map_contacts([reference], paired, metadata, directory / "mapping", threads)
    options = load(metadata)["settings"]["scaffolding"]
    execute(
        [
            "haphic",
            "plot",
            Path(agp).resolve(),
            bam.resolve(),
            "--bin_size",
            options["plot_bin_kb"],
            "--min_len",
            "0",
            "--separate_plots",
            "--normalization",
            "log10",
            "--threads",
            threads,
        ],
        cwd=directory,
    )
    for name in ("contact_map.pdf", "separate_plots.pdf"):
        if not (directory / name).is_file():
            raise ValueError(f"Missing Hi-C contact plot: {name}")
    return statistics
