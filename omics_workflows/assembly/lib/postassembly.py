"""Assembly manifests, local sequence evidence and post-assembly contracts."""

import csv
import itertools
import re
import shutil
import sys
from collections import Counter, OrderedDict
from pathlib import Path

import pysam

from omics_workflows.assembly.lib.io import checksum, execute, load, pipeline, save, signature
from omics_workflows.assembly.lib.sequences import fasta_dict, fasta_statistics, write_fasta


def reverse_complement(sequence):
    """返回由大写ACGTN组成序列的反向互补字符串。"""
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_sequences(path, sequences):
    """把序列ID到碱基字符串的映射写为FASTA，保留映射迭代顺序。"""
    with Path(path).open("w", encoding="utf-8") as handle:
        for name, seq in sequences.items():
            write_fasta(handle, name, seq)


def terminal_repeats(sequence, motif, window=10000, distance=1000, repeats=10):
    """在真实序列末端附近查找连续端粒阵列，返回左右端状态。

    搜索重复单元的循环移位及两种方向；内部命中不算末端，间断阵列保留未解决状态，
    即使tidk在附近检出重复也不直接接受为可完善位点。
    """
    motifs = {motif[i:] + motif[:i] for i in range(len(motif))}
    rc = reverse_complement(motif)
    motifs |= {rc[i:] + rc[:i] for i in range(len(rc))}
    result = {}
    for side, region in [("left", sequence[:window]), ("right", sequence[-window:])]:
        hits = []
        for m in motifs:
            for match in re.finditer("(?:" + m + "){" + str(repeats) + ",}", region):
                offset = match.start() if side == "left" else len(region) - match.end()
                if offset <= distance:
                    hits.append({"distance": offset, "array_bases": match.end() - match.start(), "motif": m})
        result[side] = max(hits, key=lambda h: h["array_bases"]) if hits else None
    return result


def inspect_sequences(sequences, options, motif):
    """按序列汇总长度、缺口和左右端粒状态，返回记录列表，供完善前后比较。"""
    result = []
    for name, seq in sequences.items():
        ends = terminal_repeats(
            seq,
            motif,
            options["terminal_window"],
            options["terminal_max_distance"],
            options["terminal_min_repeats"],
        )
        gaps = [
            {"start": m.start(), "end": m.end(), "length": m.end() - m.start()}
            for m in re.finditer("N+", seq)
        ]
        ambiguity = sum(base not in "ACGT" for base in seq)
        result.append(
            {
                "sequence_id": name,
                "length": len(seq),
                "gaps": gaps,
                "N_bases": seq.count("N"),
                "ambiguous_bases": ambiguity,
                "telomeres": ends,
                "gap_free_two_telomeres": not ambiguity and bool(ends["left"]) and bool(ends["right"]),
            }
        )
    return result


def map_reads(reference, reads, output, threads, preset="map-hifi", _min_mapq=None):
    """把长读段回贴到参考，生成排序并建立索引的BAM；预设由HiFi或ONT技术决定。

    比对和排序并发共享线程预算，避免两个进程各自占满声明资源。
    """
    reference, output = Path(reference).resolve(), Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    commands = [
        [
            "minimap2",
            "-a",
            "-x",
            preset,
            "--secondary=yes",
            "-N",
            "10",
            "-t",
            max(1, threads - 2),
            reference,
            *[str(Path(p).resolve()) for p in reads],
        ],
        ["samtools", "sort", "-@", "1", "-m", "512M", "-o", output, "-"],
    ]
    # minimap2 treats multiple FASTQ arguments as paired fragments. Long-read
    # libraries are therefore streamed as ONE concatenated input, never as mates.
    if len(reads) > 1:
        stream = (
            "import sys,gzip,shutil\n"
            "for p in sys.argv[1:]:\n"
            ' with open(p,"rb") as h: magic=h.read(2)\n'
            " opener=gzip.open if magic==bytes([31,139]) else open\n"
            ' with opener(p,"rb") as h: shutil.copyfileobj(h,sys.stdout.buffer,8388608)\n'
        )
        commands[0] = commands[0][: -len(reads)] + ["-"]
        commands[0][commands[0].index("-t") + 1] = max(1, threads - 3)
        commands.insert(0, [sys.executable, "-c", stream, *[str(Path(p).resolve()) for p in reads]])
    sink = output.with_suffix(".pipeline.stdout")
    pipeline(commands, sink)
    sink.unlink()
    execute(["samtools", "quickcheck", "-v", output])
    execute(["samtools", "index", output])


def supported_alignments(bam, chrom, start, end, min_mapq, min_identity):
    """按区间迭代唯一、主要且达到MAPQ和一致性要求的比对，排除不可靠或重复支持。"""
    seen = set()
    for record in bam.fetch(chrom, max(0, start), end):
        if (
            record.is_unmapped
            or record.is_secondary
            or record.is_supplementary
            or record.mapping_quality < min_mapq
        ):
            continue
        aligned = record.query_alignment_length or 0
        if not aligned or not record.has_tag("NM") or 1 - record.get_tag("NM") / aligned < min_identity:
            continue
        if record.query_name in seen:
            raise ValueError("Duplicate primary long-read ID: " + record.query_name)
        seen.add(record.query_name)
        yield record


def write_bundle(stage, sequences, groups):
    """写入全集和按单套/单倍型分组的FASTA及索引，返回各集合统计；groups定义每条序列归属。"""
    write_sequences(stage / "genome.fa", sequences)
    execute(["samtools", "faidx", stage / "genome.fa"])
    stats = {"genome": fasta_statistics(stage / "genome.fa")}
    for hap in sorted(set(groups.values())):
        write_sequences(
            stage / (hap + ".fa"), OrderedDict((n, s) for n, s in sequences.items() if groups[n] == hap)
        )
        execute(["samtools", "faidx", stage / (hap + ".fa")])
        stats[hap] = fasta_statistics(stage / (hap + ".fa"))
    if (stage / "unplaced.fa").is_file():
        stats["unplaced"] = fasta_statistics(stage / "unplaced.fa")
    return stats


def write_checksums(stage):
    """为暂存目录中的普通文件写SHA256SUMS，相对路径记录；跳过符号链接和校验表自身。"""
    with (stage / "SHA256SUMS").open("w", encoding="utf-8") as out:
        for f in sorted(stage.rglob("*")):
            if f.is_file() and not f.is_symlink() and f.name != "SHA256SUMS":
                out.write(f"{checksum(f)}  {f.relative_to(stage)}\n")


def rebase_symlinks(stage, destination):
    """在发布目录前调整内部链接，使移动后的证据仍可访问；外部原始输入链接保持原目标。"""
    stage = Path(stage).resolve()
    links = {}
    for path in stage.rglob("*"):
        if not path.is_symlink():
            continue
        target = path.resolve()
        if target.is_relative_to(stage):
            target = Path(destination).resolve() / target.relative_to(stage)
            path.unlink()
            path.symlink_to(target)
        links[str(path.relative_to(stage))] = str(target)
    if links:
        save(stage / "symlinks.json", links)


def native_tidk(genome, outdir, motif):
    """按给定重复单元调用tidk ends生成原生端粒结果；不把端粒检出等同于T2T认证。"""
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    execute(
        [
            "tidk",
            "search",
            "--string",
            motif,
            "--window",
            "1000",
            "--output",
            "repeats",
            "--dir",
            outdir,
            Path(genome).resolve(),
        ]
    )
    if not list(outdir.glob("repeats*")):
        raise ValueError("tidk completed without repeat output")


def artifact(path, base):
    """返回文件相对路径、SHA-256和生产时签名组成的记录，用于后续步骤核验来源。"""
    path = Path(path)
    return {"file": str(path.relative_to(base)), "sha256": checksum(path), "signature": signature(path)}


def checked_path(base, record):
    """根据基准目录和清单记录返回有效文件路径。

    先核对生产时签名，变化后才重算校验值，避免重复扫描大体积读段；缺失或损坏时报错。
    """
    path = (Path(base) / record["file"]).absolute()
    if not path.is_file():
        raise ValueError(f"Missing declared artifact: {path}")
    if record.get("signature") != signature(path):
        if checksum(path) != record["sha256"]:
            raise ValueError(f"Artifact checksum mismatch: {path}")
    return path


def read_bundle(path, sample=None):
    """读取并核验PASS状态的读段清单，检查样本身份和文件，返回补充实际路径的读段记录。"""
    result = load(path)
    if sample and result["sample_id"] != sample:
        raise ValueError("Read manifest belongs to a different sample")
    if result.get("status") != "PASS":
        raise ValueError("Read preparation or quality control did not pass")
    base = Path(path).parent
    for record in result["readsets"].values():
        record["paths"] = [str(checked_path(base, item)) for item in record["files"]]
    return result


def fastq_statistics(paths):
    """完整解析单端或双端FASTQ，返回读段数、碱基数、平均长度和读段N50。

    双端逐条核对数量、ID和顺序，只保留必要计数，不把格式解析替代测序质量过滤。
    """
    count = bases = 0
    lengths = Counter()
    readers = [pysam.FastxFile(str(path)) for path in paths]
    try:
        for records in itertools.zip_longest(*readers):
            if any(record is None for record in records):
                raise ValueError("Paired FASTQ record counts differ")
            ids = [record.name.removesuffix("/1").removesuffix("/2") for record in records]
            if len(set(ids)) != 1:
                raise ValueError(f"Paired FASTQ IDs/order differ: {ids}")
            for record in records:
                length = len(record.sequence)
                if not length or record.quality is None or len(record.quality) != length:
                    raise ValueError(f"Malformed FASTQ sequence/quality: {record.name}")
                count += 1
                bases += length
                lengths[length] += 1
    finally:
        for reader in readers:
            reader.close()
    if not count:
        raise ValueError("Empty FASTQ after preparation/filtering")
    accumulated = n50 = 0
    for length, number in sorted(lengths.items(), reverse=True):
        accumulated += length * number
        if accumulated >= bases / 2:
            n50 = length
            break
    return {
        "num_reads": count,
        "num_bases": bases,
        "mean_length": bases / count,
        "read_n50": n50,
        "num_pairs": count // 2 if len(paths) == 2 else None,
    }


def write_table(path, fields, rows):
    """按fields列顺序把字典记录写成TSV，输出表头并忽略未声明的额外键。"""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def assembly_bundle(path):
    """读取染色体或后处理交付，核验FASTA、索引、清单及分组，返回清单、序列和归属字典。

    拒绝不完整挂载结果；依据显式策略和分组读取，不从文件名推断单倍型或证明相位。
    """
    result = load(path)
    base = Path(path).parent
    if result.get("quality_level") == "INCOMPLETE_SCAFFOLDING":
        raise ValueError("Incomplete chromosome scaffolding requires review")
    genome_label = "genome" if "genome" in result["sets"] else "chromosomes"
    genome = base / f"{genome_label}.fa"
    if checksum(genome) != result["sets"][genome_label]["sha256"]:
        raise ValueError("Assembly FASTA changed after delivery")
    sequence = fasta_dict(genome)
    groups = result.get("groups", {})
    if not groups:
        for label in result["sets"]:
            if label.endswith(".chromosomes"):
                group = label.split(".")[0]
                for name, value in fasta_dict(base / f"{label}.fa").items():
                    if name in groups or sequence.get(name) != value:
                        raise ValueError("Overlapping or inconsistent sequence sets")
                    groups[name] = group
    if set(groups) != set(sequence):
        raise ValueError("Assembly must assign every chromosome to an explicit sequence set")
    for label in set(groups.values()):
        key = label if genome_label == "genome" else label + ".chromosomes"
        copy = base / f"{key}.fa"
        if checksum(copy) != result["sets"][key]["sha256"]:
            raise ValueError(f"Sequence set checksum mismatch: {label}")
        expected = {name: value for name, value in sequence.items() if groups[name] == label}
        if fasta_dict(copy) != expected:
            raise ValueError(f"Sequence set differs from combined chromosomes: {label}")
    unplaced = base / "unplaced.fa"
    if not unplaced.is_file():
        raise ValueError("Assembly lacks the declared unplaced FASTA")
    if "unplaced" in result["sets"] and checksum(unplaced) != result["sets"]["unplaced"]["sha256"]:
        raise ValueError("Unplaced FASTA checksum mismatch")
    result["groups"] = groups
    result["genome_path"] = str(genome.resolve())
    result["base"] = str(base.resolve())
    return result, sequence, groups


def assembly_identity(source):
    """从组装清单提取样本、真实倍性、策略、输出份数和每份染色体数，返回身份字典。"""
    return {
        key: source[key]
        for key in ("sample_id", "ploidy", "strategy", "output_copies", "expected_chromosomes_per_haplotype")
        if key in source
    }


def passthrough(source_path, destination, status):
    """把未修改的组装链接到下一编号步骤，写入指定状态并返回清单，避免物理复制整份基因组。"""
    source, _, groups = assembly_bundle(source_path)
    base = Path(source["base"])
    genome_label = "genome" if "genome" in source["sets"] else "chromosomes"
    (destination / "genome.fa").symlink_to((base / f"{genome_label}.fa").resolve())
    (destination / "genome.fa.fai").symlink_to((base / f"{genome_label}.fa.fai").resolve())
    for label in set(groups.values()):
        key = label if genome_label == "genome" else label + ".chromosomes"
        for suffix in (".fa", ".fa.fai"):
            (destination / (label + suffix)).symlink_to((base / (key + suffix)).resolve())
    (destination / "unplaced.fa").symlink_to((base / "unplaced.fa").resolve())
    sets = {
        label: fasta_statistics(destination / f"{label}.fa")
        for label in ["genome", *sorted(set(groups.values())), "unplaced"]
    }
    result = {
        **assembly_identity(source),
        "status": status,
        "groups": groups,
        "sets": sets,
        "source_manifest": str(Path(source_path).resolve()),
        "source_manifest_sha256": checksum(source_path),
    }
    save(destination / "manifest.json", result)
    return result


def retain_files(directory, names):
    """仅在流程拥有的暂存目录中保留names列出的文件，清理其余工具工作目录和中间文件。"""
    for path in Path(directory).iterdir():
        if path.name not in names:
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink()
