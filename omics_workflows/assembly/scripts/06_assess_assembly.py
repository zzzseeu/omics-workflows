#!/usr/bin/env python3
"""Seven final assembly checks with fixed contracts and compact scientific outputs."""
import re
import shutil

# Direct script execution also works outside the repository.
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.resources import workflow_defaults
from omics_workflows.assembly.lib.assessment import (
    busco_metrics,
    lai_metrics,
    mapping_metrics,
    merqury_metrics,
    quast_metrics,
    telomere_metrics,
)
from omics_workflows.assembly.lib.hic import plot_contacts
from omics_workflows.assembly.lib.io import assembly_scratch, checksum, load, publish, save
from omics_workflows.assembly.lib.postassembly import (
    artifact,
    assembly_bundle,
    assembly_identity,
    read_bundle,
    write_sequences,
    write_table,
)
from omics_workflows.assembly.lib.sequences import fasta_dict

CATEGORIES = ("continuity", "genes", "kmers", "read_support", "ltr", "chromosome_structure", "telomeres")


def sequence_sets(source, genome, groups, work):
    """按显式分组写单套/单倍型及含未挂载序列的全集FASTA，返回集合到路径的映射。"""
    paths = {}
    for label in sorted(set(groups.values())):
        paths[label] = work / f"{label}.fa"
        write_sequences(
            paths[label], {name: sequence for name, sequence in genome.items() if groups[name] == label}
        )
    competitive = dict(genome)
    for name, sequence in fasta_dict(Path(source["base"]) / "unplaced.fa").items():
        if name in competitive:
            raise ValueError(f"Chromosome and unplaced sequence IDs collide: {name}")
        competitive[name] = sequence
    combined = work / "assembly_all.fa"
    write_sequences(combined, competitive)
    paths["assembly_all"] = combined
    return paths


def copy_metric_file(source, target):
    """创建目标父目录并复制必要的原生统计文件，不改写其内容。"""
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def structural_metrics(source, genome, reads, metadata, work, output, threads):
    """以最终FASTA建立等价AGP重绘Hi-C图，并返回染色体数、挂载比例和未挂载长度。

    不沿用完善/纠错前的旧坐标；接触图仍保留未人工审阅状态。
    """
    agp = work / "final.agp"
    with agp.open("w", encoding="utf-8") as handle:
        for name, sequence in genome.items():
            handle.write(f"{name}\t1\t{len(sequence)}\t1\tW\t{name}\t1\t{len(sequence)}\t+\n")
    statistics = plot_contacts(work / "assembly_all.fa", agp, reads, metadata, work / "hic", threads)
    output.mkdir()
    for name in ("contact_map.pdf", "separate_plots.pdf"):
        copy_metric_file(work / "hic" / name, output / name)
    nuclear = sum(len(sequence) - sequence.count("N") for sequence in genome.values())
    unplaced = fasta_dict(Path(source["base"]) / "unplaced.fa")
    total = nuclear + sum(len(sequence) - sequence.count("N") for sequence in unplaced.values())
    return {
        "status": "COMPLETED",
        "chromosome_counts": dict(Counter(source["groups"].values())),
        "anchored_non_N_fraction": nuclear / total if total else None,
        "unplaced_length": sum(map(len, unplaced.values())),
        "hic_mapping": statistics,
        "human_review": "NOT_REVIEWED",
        "coordinates": "final FASTA identity AGP",
    }


def validate_telomeres(rows, genome):
    """核验tidk逐末端记录的序列、左右端、状态及唯一性，返回以序列和端点为键的字典。"""
    ends = {}
    for row in rows:
        if row.get("id") not in genome or row.get("side") not in ("start", "end"):
            raise ValueError("tidk output has an unexpected chromosome or end")
        if row["status"] not in ("present", "absent", "wrong_strand", "not_terminal"):
            raise ValueError("Unsupported tidk end status")
        key = (row["id"], row["side"])
        if key in ends:
            raise ValueError("Duplicate tidk chromosome-end record")
        ends[key] = row
    if len(ends) != len(genome) * 2:
        raise ValueError("tidk did not evaluate both ends of every chromosome")
    return ends


def run_checks(args, source, genome, groups, reads, metadata, stage, work):
    """执行七类固定评估，保留必要原生结果并写metrics.tsv，返回各检查的状态及指标字典。"""
    paths = sequence_sets(source, genome, groups, work)
    metrics = []
    evaluations = {}
    copies = {label: path for label, path in paths.items() if label != "assembly_all"}
    all_sets = {"genome": Path(source["genome_path"]), **paths}
    quast = {}
    for label, fasta in all_sets.items():
        folder = work / "quast" / label
        quast[label] = quast_metrics(fasta, label, folder, args.threads, args.quast_min_contig)
        copy_metric_file(folder / "report.tsv", stage / "statistics/quast" / f"{label}.tsv")
        for metric in ("N50", "L50", "Total length", "Largest contig", "GC (%)", "# N's per 100 kbp"):
            for assembly_name, value in quast[label].get(metric, {}).items():
                metrics.append(
                    {
                        "sequence_set": assembly_name,
                        "metric": metric,
                        "value": value,
                        "unit": "QUAST_native",
                        "category": "continuity",
                        "status": "COMPLETED",
                    }
                )
    evaluations["continuity"] = {
        "status": "COMPLETED",
        "sets": quast,
        "min_contig_length": args.quast_min_contig,
    }
    buscos = {}
    for label, fasta in copies.items():
        values, summary, tables = busco_metrics(
            fasta, label, work / "busco", args.busco_lineage, args.threads
        )
        buscos[label] = values
        copy_metric_file(summary, stage / "statistics/busco" / f"{label}.json")
        if len(tables) != 1:
            raise ValueError("BUSCO full table missing or ambiguous")
        copy_metric_file(tables[0], stage / "statistics/busco" / f"{label}.tsv")
        metrics.extend(
            {
                "sequence_set": label,
                "metric": symbol,
                "value": value,
                "unit": "%",
                "category": "genes",
                "status": "COMPLETED",
            }
            for symbol, value in values.items()
        )
    evaluations["genes"] = {
        "status": "COMPLETED",
        "sets": buscos,
        "lineage": str(Path(args.busco_lineage).resolve()),
        "lineage_sha256": checksum(Path(args.busco_lineage) / "dataset.cfg"),
    }
    wgs = [
        path
        for record in reads["readsets"].values()
        if record["technology"] == "wgs"
        for path in record["paths"]
    ]
    if wgs:
        values = merqury_metrics(
            Path(source["genome_path"]),
            paths,
            wgs,
            work / "merqury",
            args.threads,
            args.mem_mb,
            args.kmer_size,
        )
        evaluations["kmers"] = {
            "status": "COMPLETED",
            "sets": values,
            "k": args.kmer_size,
            "scope": "haplotype completeness uses whole-sample read k-mers",
        }
        for label, value in values.items():
            metrics.extend(
                [
                    {
                        "sequence_set": label,
                        "metric": "QV",
                        "value": value["qv"],
                        "unit": "QV",
                        "category": "kmers",
                        "status": "COMPLETED",
                    },
                    {
                        "sequence_set": label,
                        "metric": "kmer_completeness",
                        "value": value["completeness_percent"],
                        "unit": "%",
                        "category": "kmers",
                        "status": "COMPLETED",
                    },
                ]
            )
            for path in (work / "merqury" / label).glob("assessment*"):
                if path.is_file() and (path.suffix in (".qv", ".stats", ".hist", ".pdf", ".png")):
                    copy_metric_file(path, stage / "statistics/merqury" / label / path.name)
    else:
        evaluations["kmers"] = {"status": "MISSING_DATA", "reason": "Ordinary WGS not supplied"}
    hifi = [
        path
        for record in reads["readsets"].values()
        if record["technology"] == "hifi"
        for path in record["paths"]
    ]
    evaluations["read_support"] = mapping_metrics(
        paths["assembly_all"],
        hifi,
        genome,
        work / "mapping",
        args.threads,
        args.min_mapq,
        args.coverage_window,
    )
    for name in ("chromosome_coverage.tsv", "coverage.regions.bed.gz", "coverage.mosdepth.global.dist.txt"):
        copy_metric_file(work / "mapping" / name, stage / "statistics/mapping" / name)
    for metric in ("mapping_percent", "split_read_percent", "covered_fraction"):
        metrics.append(
            {
                "sequence_set": "genome",
                "metric": metric,
                "value": evaluations["read_support"][metric],
                "unit": "%" if metric.endswith("percent") else "fraction",
                "category": "read_support",
                "status": "COMPLETED",
            }
        )
    lais = {}
    for label, fasta in copies.items():
        values, table = lai_metrics(
            fasta,
            {name: sequence for name, sequence in genome.items() if groups[name] == label},
            work / "lai" / label,
            args.threads,
        )
        lais[label] = values
        if table:
            copy_metric_file(table, stage / "statistics/lai" / f"{label}.tsv")
        metrics.append(
            {
                "sequence_set": label,
                "metric": "LAI",
                "value": values.get("LAI", ""),
                "unit": "LAI",
                "category": "ltr",
                "status": values["status"],
            }
        )
    evaluations["ltr"] = {
        "status": (
            "COMPLETED"
            if all(value["status"] == "COMPLETED" for value in lais.values())
            else "NOT_APPLICABLE"
        ),
        "sets": lais,
    }
    evaluations["chromosome_structure"] = structural_metrics(
        source, genome, reads, metadata, work, stage / "contact_maps", args.threads
    )
    rows, summary, table = telomere_metrics(
        Path(source["genome_path"]),
        work / "telomeres",
        args.telomere_motif,
        args.terminal_window,
        args.telomere_min_length,
        args.terminal_max_distance,
    )
    ends = validate_telomeres(rows, genome)
    copy_metric_file(table, stage / "telomeres.tsv")
    per_copy = {}
    for label in copies:
        names = [name for name, group in groups.items() if group == label]
        per_copy[label] = {
            "chromosomes": len(names),
            "present_ends": sum(
                ends[(name, side)]["status"] == "present" for name in names for side in ("start", "end")
            ),
            "both_ends": sum(
                all(ends[(name, side)]["status"] == "present" for side in ("start", "end")) for name in names
            ),
        }
    evaluations["telomeres"] = {"status": "COMPLETED", "native_summary": summary, "sets": per_copy}
    write_table(
        stage / "chromosomes.tsv",
        [
            "sequence_id",
            "sequence_set",
            "length",
            "N_bases",
            "gap_count",
            "ambiguous_bases",
            "left_telomere",
            "right_telomere",
        ],
        [
            {
                "sequence_id": name,
                "sequence_set": groups[name],
                "length": len(sequence),
                "N_bases": sequence.count("N"),
                "gap_count": len(re.findall("N+", sequence)),
                "ambiguous_bases": sum(base not in "ACGT" for base in sequence),
                "left_telomere": ends[(name, "start")]["status"],
                "right_telomere": ends[(name, "end")]["status"],
            }
            for name, sequence in genome.items()
        ],
    )
    for label, values in per_copy.items():
        metrics.extend(
            {
                "sequence_set": label,
                "metric": key,
                "value": value,
                "unit": "count",
                "category": "telomeres",
                "status": "COMPLETED",
            }
            for key, value in values.items()
        )
    structure = evaluations["chromosome_structure"]
    metrics.extend(
        {
            "sequence_set": "genome",
            "metric": key,
            "value": structure[key],
            "unit": "fraction" if key.endswith("fraction") else "bp",
            "category": "chromosome_structure",
            "status": "COMPLETED",
        }
        for key in ("anchored_non_N_fraction", "unplaced_length")
    )
    for category, result in evaluations.items():
        if result["status"] != "COMPLETED":
            metrics.append(
                {
                    "sequence_set": "genome",
                    "metric": "availability",
                    "value": "",
                    "unit": "",
                    "category": category,
                    "status": result["status"],
                }
            )
    write_table(
        stage / "metrics.tsv", ["category", "sequence_set", "metric", "value", "unit", "status"], metrics
    )
    return evaluations


def assess(args):
    """核验最终染色体交付、读段和输出份数，执行评估并发布同版本质量结果及逐染色体统计。

    将缺数据、不适用和人工审阅分别记录，工具失败抛出异常，不发布新的完成清单。
    """
    source, genome, groups = assembly_bundle(args.assembly_manifest)
    reads = read_bundle(args.reads_manifest, source["sample_id"])
    if not (Path(args.busco_lineage) / "dataset.cfg").is_file():
        raise ValueError("BUSCO requires an offline lineage containing dataset.cfg")
    # 新结果显式声明份数；旧清单只在缺少该字段时按已有分组兼容。
    counts = Counter(groups.values())
    expected = source["expected_chromosomes_per_haplotype"]
    declared_copies = source.get("output_copies", len(counts))
    if len(counts) != declared_copies or any(count != expected for count in counts.values()):
        raise ValueError("Chromosome sequence sets do not match declared output copies")
    metadata = (
        load(args.metadata)
        if args.metadata
        else {
            "sample_id": source["sample_id"],
            "tools": {},
            "settings": yaml.safe_load(workflow_defaults("assembly")),
        }
    )
    metadata["settings"]["scaffolding"]["hic_reference"] = "combined"
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        work = Path(temporary) / "work"
        work.mkdir()
        stage = Path(temporary) / "result"
        stage.mkdir()
        metadata_path = work / "metadata.json"
        save(metadata_path, metadata)
        evaluations = run_checks(args, source, genome, groups, reads, metadata_path, stage, work)
        if set(evaluations) != set(CATEGORIES):
            raise ValueError("Final assessment must report exactly seven categories")
        result = {
            **assembly_identity(source),
            "execution_status": "COMPLETED",
            "status": (
                "PARTIAL"
                if any(value["status"] == "MISSING_DATA" for value in evaluations.values())
                else "COMPLETED"
            ),
            "genome_sha256": checksum(source["genome_path"]),
            "files": [artifact(path, stage) for path in sorted(stage.rglob("*")) if path.is_file()],
            "evaluations": evaluations,
            "quality_status": "NEEDS_REVIEW",
            "T2T_status": "NOT_CERTIFIED",
            "phasing": "NOT_APPLICABLE" if source["strategy"] == "single" else "NOT_INDEPENDENTLY_VALIDATED",
            "assembly_manifest": str(Path(args.assembly_manifest).resolve()),
            "assembly_manifest_sha256": checksum(args.assembly_manifest),
            "reads_manifest": str(Path(args.reads_manifest).resolve()),
            "reads_manifest_sha256": checksum(args.reads_manifest),
            "parameters": {key: value for key, value in vars(args).items() if key != "metadata"},
            "tools": metadata.get("tools", {}),
            "code": metadata.get("code", {}),
            "limitations": [
                (
                    "Input reads may also have been used in assembly/polishing; these are "
                    "not independent validation"
                ),
                "Low scores, missing data and failed tool execution are distinct states",
                "Hi-C review, phase evidence and complex-repeat validation are reported separately",
            ],
        }
        save(stage / "assessment.json", result)
        if args.keep_intermediates:
            shutil.move(work, stage / "intermediates")
        publish(stage, destination, temporary)


def assess_contigs(args):
    """为HiFi-only终点汇总primary/alternate的连续性统计，明确标记未运行染色体七类评估。"""
    source = load(args.assembly_manifest)
    base = Path(args.assembly_manifest).parent
    from omics_workflows.assembly.lib.sequences import fasta_statistics

    rows = []
    statistics = {}
    for label in ("primary", "alternate"):
        path = base / (label + ".fa")
        if path.is_file():
            statistics[label] = fasta_statistics(path)
            for metric, value in statistics[label].items():
                if metric != "sha256":
                    rows.append(
                        {
                            "category": "continuity",
                            "sequence_set": label,
                            "metric": metric,
                            "value": value,
                            "unit": "FASTA_native",
                            "status": "COMPLETED",
                        }
                    )
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        write_table(
            stage / "metrics.tsv", ["category", "sequence_set", "metric", "value", "unit", "status"], rows
        )
        save(
            stage / "assessment.json",
            {
                "sample_id": source["sample_id"],
                "endpoint": "contigs",
                "status": "CONTIG_STATISTICS_ONLY",
                "execution_status": "COMPLETED",
                "sets": statistics,
                "assembly_manifest_sha256": checksum(args.assembly_manifest),
                "limitations": [
                    "Chromosome-scale seven-category final assessment not run for contig endpoint"
                ],
            },
        )
        publish(stage, destination, temporary)


def main():
    """解析步骤06参数，按contigs-only或染色体路线分派检查，并校验资源及分析阈值。"""
    parser = argument_parser(description="步骤06：评估最终组装质量，或仅汇总contig统计；不修改组装序列。")
    parser.add_argument(
        "--assembly-manifest",
        required=True,
        help="类型：JSON文件路径；必填。最终组装结果清单，评估与交付序列须为同一版本。",
    )
    parser.add_argument(
        "--reads-manifest",
        required=True,
        help="类型：JSON 文件路径；必填。步骤02质控后的读段清单，含样本身份、文件位置和校验值。",
    )
    parser.add_argument(
        "--busco-lineage",
        required=True,
        help="类型：目录路径；必填。适合物种的离线BUSCO lineage，须含dataset.cfg；--contigs-only模式不使用。",
    )
    parser.add_argument(
        "--telomere-motif",
        required=True,
        help="类型：字符串，由2至20个大写ACGT组成；必填。如水稻常用TTTAGGG，应按材料背景指定；--contigs-only模式不使用。",
    )
    parser.add_argument(
        "--metadata", help="类型：JSON 文件路径；可选。沿用流程生成的分析设置；未提供时读取模块默认设置。"
    )
    parser.add_argument(
        "--out-dir", required=True, help="类型：目录路径；必填。保存本步骤的核心结果、统计表和结果清单。"
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="类型：整数；默认8。CPU线程数；染色体最终评估至少4。"
    )
    parser.add_argument(
        "--mem-mb",
        type=int,
        default=64000,
        help="类型：整数，单位MB；默认64000，染色体最终评估至少1024。提供给评估工具的内存预算。",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=21,
        help="类型：整数，范围15至31；默认21。meryl/Merqury评估使用的k-mer长度。",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=1, help="类型：整数，范围0至60；默认1。HiFi覆盖统计的最低比对MAPQ。"
    )
    parser.add_argument(
        "--coverage-window",
        type=int,
        default=10000,
        help="类型：正整数，单位bp；默认10000。mosdepth覆盖深度统计窗口长度。",
    )
    parser.add_argument(
        "--quast-min-contig",
        type=int,
        default=1,
        help="类型：正整数，单位bp；默认1。QUAST统计纳入的最短序列长度。",
    )
    parser.add_argument(
        "--terminal-window",
        type=int,
        default=10000,
        help="类型：正整数，单位bp；默认10000。tidk检查的染色体末端窗口长度。",
    )
    parser.add_argument(
        "--terminal-max-distance",
        type=int,
        default=1000,
        help="类型：非负整数，单位bp；默认1000，须小于--terminal-window。端粒重复阵列距序列端点的最大允许距离。",
    )
    parser.add_argument(
        "--telomere-min-length",
        type=int,
        default=200,
        help="类型：正整数，单位bp；默认200。tidk检出的端粒重复阵列最短长度。",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。保留中间文件用于排错，增加磁盘占用。",
    )
    parser.add_argument(
        "--contigs-only",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。仅统计contig连续性，不执行染色体七类最终评估；仍需提供接口中的必填参数。",
    )
    args = parser.parse_args()
    if args.contigs_only:
        assess_contigs(args)
        return
    valid_sizes = min(args.coverage_window, args.quast_min_contig, args.telomere_min_length) > 0
    valid_telomere = 0 <= args.terminal_max_distance < args.terminal_window and bool(
        re.fullmatch("[ACGT]{2,20}", args.telomere_motif)
    )
    valid_resources = args.threads >= 4 and args.mem_mb >= 1024
    if (
        not valid_resources
        or not 15 <= args.kmer_size <= 31
        or not 0 <= args.min_mapq <= 60
        or not valid_sizes
        or not valid_telomere
    ):
        parser.error("Invalid resources, k-mer, mapping or telomere parameters")
    assess(args)


if __name__ == "__main__":
    main()
