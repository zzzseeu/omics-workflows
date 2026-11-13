#!/usr/bin/env python3
"""NextPolish2 with ordinary WGS k-mers and explicit missing-data records."""
import gzip
import shutil

# Direct script execution also works outside the repository.
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.resources import workflow_defaults
from omics_workflows.assembly.lib.io import assembly_scratch, checksum, execute, load, publish, save
from omics_workflows.assembly.lib.postassembly import (
    assembly_bundle,
    assembly_identity,
    map_reads,
    passthrough,
    read_bundle,
    retain_files,
    write_bundle,
)
from omics_workflows.assembly.lib.sequences import fasta_dict


def merge_wgs(paths, destination):
    """按给定顺序合并普通WGS为有效gzip流，兼容原始文本和gzip输入，不重复做质量过滤。"""
    with destination.open("wb") as target:
        for path in paths:
            with Path(path).open("rb") as handle:
                compressed = handle.read(2) == b"\x1f\x8b"
                handle.seek(0)
                if compressed:
                    shutil.copyfileobj(handle, target)
                else:
                    with gzip.GzipFile(fileobj=target, mode="wb", compresslevel=3) as member:
                        shutil.copyfileobj(handle, member)


def audit_polished(before, after, work, threads):
    """比对纠错前后序列，核对名称、集合、方向和覆盖一致性，返回通过检查的新序列字典。"""
    old, new = fasta_dict(before), fasta_dict(after)
    if set(old) != set(new):
        raise ValueError("Polishing changed chromosome IDs or copy count")
    for name, sequence in new.items():
        if abs(len(sequence) - len(old[name])) > max(1000, len(old[name]) * 0.05):
            raise ValueError("Unexpected length change after polishing: " + name)
    paf = work / "polishing_changes.paf"
    with paf.open("w", encoding="utf-8") as out:
        execute(["minimap2", "-c", "--cs=long", "-x", "asm5", "-t", threads, before, after], stdout=out)
    intervals = defaultdict(list)
    for line in paf.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        if f[0] != f[5] or f[4] != "+":
            raise ValueError("Polished sequence changes copy/strand assignment; review required")
        intervals[f[0]].append((int(f[2]), int(f[3])))
    for name, seq in new.items():
        end = covered = 0
        for start, stop in sorted(intervals[name]):
            covered += max(0, stop - max(start, end))
            end = max(end, stop)
        if covered < len(seq) * 0.95:
            raise ValueError("Insufficient same-copy polishing concordance: " + name)
    return new


def polish(args):
    """从清洁WGS建立yak计数并用HiFi和NextPolish2纠错，发布核验后的序列及前后对应。

    按设置区分关闭、主动跳过和缺WGS，不能用Hi-C替代普通WGS，也不静默改变拷贝归属。
    """
    defaults = yaml.safe_load(workflow_defaults("assembly"))
    meta = load(args.metadata) if args.metadata else {"settings": defaults, "tools": {}}
    options = meta["settings"]["polishing"]
    source, _, groups = assembly_bundle(args.finishing)
    base = Path(source["base"])
    if args.min_contig_length is not None:
        options["min_contig_length"] = args.min_contig_length
    if args.on_missing_wgs is not None:
        options["on_missing_wgs"] = args.on_missing_wgs
    wgs = []
    if args.reads_manifest:
        reads = read_bundle(args.reads_manifest, source["sample_id"])
        args.hifi = [
            p
            for record in reads["readsets"].values()
            if record["technology"] == "hifi"
            for p in record["paths"]
        ]
        wgs = [
            p
            for record in reads["readsets"].values()
            if record["technology"] == "wgs"
            for p in record["paths"]
        ]
    elif args.wgs_r1 or args.wgs_r2:
        if len(args.wgs_r1) != len(args.wgs_r2):
            raise ValueError("WGS R1/R2 lists must have matching library counts")
        from omics_workflows.assembly.lib.postassembly import fastq_statistics

        for first, second in zip(args.wgs_r1, args.wgs_r2):
            fastq_statistics([first, second])
            wgs.extend([first, second])
    enabled = options["enabled"] and not args.skip if args.metadata else not args.skip
    if enabled and (args.threads < 4 or not 0 <= args.min_mapq <= 60 or options["min_contig_length"] < 1000):
        raise ValueError("Polishing requires >=4 threads, MAPQ 0..60 and min_contig_length>=1000")
    if enabled and not wgs and options["on_missing_wgs"] == "error":
        raise ValueError("NextPolish2 requires ordinary WGS reads; Hi-C cannot substitute")
    if any(k < 15 or k > 31 for k in args.kmer_sizes) or len(set(args.kmer_sizes)) != len(args.kmer_sizes):
        raise ValueError("yak k-mer sizes must be unique integers from 15 through 31")
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        if not enabled or not wgs:
            passthrough(args.finishing, stage, "DISABLED" if not enabled else "SKIPPED_MISSING_WGS")
        else:
            if not args.hifi:
                raise ValueError("NextPolish2 requires HiFi reads")
            genome = Path(source["genome_path"])
            merged = stage / "wgs.fastq.gz"
            merge_wgs(wgs, merged)
            databases = []
            for kmer in args.kmer_sizes:
                path = stage / f"k{kmer}.yak"
                execute(
                    ["yak", "count", "-k", kmer, "-b", "30", "-t", args.threads, "-o", path, merged, merged]
                )
                databases.append(path)
            bam = stage / "hifi.bam"
            map_reads(genome, args.hifi, bam, args.threads)
            unique = stage / "hifi.unique.bam"
            execute(["samtools", "view", "-b", "-q", args.min_mapq, "-F", "2308", "-o", unique, bam])
            execute(["samtools", "index", unique])
            with (stage / "polished.fa").open("w", encoding="utf-8") as target:
                execute(
                    [
                        "nextPolish2",
                        "-t",
                        args.threads,
                        "-m",
                        "ref",
                        "-L",
                        options["min_contig_length"],
                        "--uppercase",
                        unique,
                        genome,
                        *databases,
                    ],
                    stdout=target,
                )
            sequence = audit_polished(genome, stage / "polished.fa", stage, args.threads)
            (stage / "unplaced.fa").symlink_to((base / "unplaced.fa").resolve())
            stats = write_bundle(stage, sequence, groups)
            save(
                stage / "manifest.json",
                {
                    **assembly_identity(source),
                    "status": "COMPLETED",
                    "groups": groups,
                    "sets": stats,
                    "parameters": {**options, "kmer_sizes": args.kmer_sizes, "min_mapq": args.min_mapq},
                    "source_manifest": str(Path(args.finishing).resolve()),
                    "source_manifest_sha256": checksum(args.finishing),
                    "tools": meta.get("tools", {}),
                    "phase_validation": (
                        "NOT_APPLICABLE" if source.get("strategy") == "single" else "NOT_ASSESSED"
                    ),
                },
            )
            if not args.keep_intermediates:
                keep = {"genome.fa", "genome.fa.fai", "unplaced.fa", "manifest.json", "polishing_changes.paf"}
                keep.update(group + suffix for group in set(groups.values()) for suffix in (".fa", ".fa.fai"))
                retain_files(stage, keep)
        publish(stage, destination, temporary)


def main():
    """解析步骤05参数并执行纠错或显式跳过，保留原分析开关和缺数据处理语义。"""
    parser = argument_parser(description="步骤05：利用清洁HiFi和普通WGS，通过NextPolish2进行碱基纠错。")
    parser.add_argument(
        "--metadata", help="类型：JSON 文件路径；可选。沿用流程生成的分析设置；未提供时读取模块默认设置。"
    )
    parser.add_argument(
        "--assembly-manifest",
        "--finishing",
        dest="finishing",
        required=True,
        help="类型：JSON文件路径；必填。待纠错的组装结果清单；--finishing是同义参数。",
    )
    parser.add_argument(
        "--reads-manifest",
        help="类型：JSON文件路径；可选。质控后读段清单；提供后自动读取HiFi与普通WGS，替代直接指定的读段列表。",
    )
    parser.add_argument(
        "--hifi",
        nargs="*",
        default=[],
        help="类型：FASTQ/FASTQ.gz文件路径列表，以空格分隔；默认空列表。未提供--reads-manifest时指定清洁HiFi读段。",
    )
    parser.add_argument(
        "--wgs-r1",
        nargs="*",
        default=[],
        help="类型：FASTQ/FASTQ.gz文件路径列表，以空格分隔；默认空列表。普通WGS的R1，按文库顺序与--wgs-r2一一对应；不接受Hi-C代替。",
    )
    parser.add_argument(
        "--wgs-r2",
        nargs="*",
        default=[],
        help="类型：FASTQ/FASTQ.gz文件路径列表，以空格分隔；默认空列表。普通WGS的R2，数量和文库顺序须与--wgs-r1一致。",
    )
    parser.add_argument(
        "--out-dir", required=True, help="类型：目录路径；必填。保存本步骤的核心结果、统计表和结果清单。"
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="类型：整数；默认8。CPU线程数；启用纠错时至少4。"
    )
    parser.add_argument(
        "--kmer-sizes",
        nargs="+",
        type=int,
        default=[21, 31],
        help="类型：整数列表，以空格分隔；默认21 31。yak计数使用的k-mer长度，各值须互不重复且在15至31之间。",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=30, help="类型：整数，范围0至60；默认30。用于纠错的HiFi比对最低MAPQ。"
    )
    parser.add_argument(
        "--min-contig-length",
        type=int,
        help="类型：整数，单位bp，至少1000。NextPolish2处理的最短序列长度；沿用metadata，无metadata时默认10000。",
    )
    parser.add_argument(
        "--on-missing-wgs",
        choices=["skip", "error"],
        help="类型：字符串枚举。skip在缺少普通WGS时跳过纠错，error报错；沿用metadata，无metadata时默认skip。",
    )
    parser.add_argument(
        "--skip",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。跳过纠错并链接原组装；有metadata时还受polishing.enabled控制，无metadata时默认尝试纠错。",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。保留中间文件用于排错，增加磁盘占用。",
    )
    polish(parser.parse_args())


if __name__ == "__main__":
    main()
