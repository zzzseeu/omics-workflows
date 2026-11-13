"""Native assembly quality tools and strict result parsers."""

import csv
import gzip
import math
import os
import re
import shutil
from pathlib import Path

import pysam

from omics_workflows.assembly.lib.io import checksum, execute, load
from omics_workflows.assembly.lib.postassembly import map_reads, write_table


def run_merqury(genome, copy_fastas, clean, work, threads, mem_mb, kmer=21):
    """用普通WGS建立meryl计数库，评估染色体全集及各序列集合，返回原生结果文件索引。

    通过局部包装限制Merqury子进程的CPU和内存，避免shell脚本超出步骤预算。
    """
    work.mkdir()
    database = work / "reads.meryl"
    execute(
        [
            "meryl",
            f"k={kmer}",
            f"threads={threads}",
            f"memory={max(1,mem_mb//2000)}",
            "count",
            "output",
            database,
            *clean,
        ]
    )
    # Merqury's shell scripts do not propagate a command-line CPU budget.
    # A per-attempt shim caps each meryl process; two can coexist in a pipe.
    tools = work / "bin"
    tools.mkdir()
    meryl = Path(shutil.which("meryl")).resolve()
    shim = tools / "meryl"
    shim.write_text(
        "#!/usr/bin/env python3\nimport os,sys\nos.execv("
        + repr(str(meryl))
        + ", ["
        + repr(str(meryl))
        + ", "
        + repr(f"threads={max(1,threads//2)}")
        + ", "
        + repr(f"memory={max(1,mem_mb//4000)}")
        + ", *sys.argv[1:]])\n",
        encoding="utf-8",
    )
    shim.chmod(0o755)
    env = dict(
        os.environ,
        MERQURY=str(Path(shutil.which("merqury.sh")).resolve().parent),
        PATH=str(tools) + os.pathsep + os.environ["PATH"],
        OMP_NUM_THREADS=str(max(1, threads // 2)),
    )
    records = {}
    for label, path in [("genome", genome), *copy_fastas.items()]:
        d = work / label
        d.mkdir()
        execute(["merqury.sh", database, path, "assessment"], cwd=d, env=env)
        qv = d / "assessment.qv"
        completeness = d / "assessment.completeness.stats"
        if not qv.is_file() or not completeness.is_file():
            raise ValueError("Merqury did not produce QV/completeness; inspect native logs")
        values = qv.read_text(encoding="utf-8").splitlines()
        if len(values) != 1:
            raise ValueError("Unexpected Merqury QV record count")
        fields = values[0].split("\t")
        value = float(fields[3])
        if math.isnan(value):
            raise ValueError("Merqury QV is undefined")
        records[label] = {
            "qv": "Infinity" if math.isinf(value) else value,
            "completeness_record": completeness.read_text(encoding="utf-8").strip(),
            "qv_file_sha256": checksum(qv),
        }
    return records


def number(value):
    """把原生指标转为整数或浮点数，拒绝NaN；允许无错误k-mer情况下的无限QV。"""
    result = float(value)
    if math.isnan(result):
        raise ValueError("Native metric is NaN")
    return int(result) if math.isfinite(result) and result.is_integer() else result


def quast_metrics(fasta, label, folder, threads, minimum):
    """运行QUAST并解析report.tsv，返回按指标和序列集合组织的值，原生统计文件保留在folder。"""
    execute(
        [
            "quast.py",
            "--large",
            "--min-contig",
            minimum,
            "--threads",
            threads,
            "--split-scaffolds",
            "--no-plots",
            "--no-html",
            "--no-icarus",
            "--labels",
            label,
            "--output-dir",
            folder,
            fasta,
        ]
    )
    report = folder / "report.tsv"
    with report.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0][0] != "Assembly" or not any(row[0] == "N50" for row in rows[1:]):
        raise ValueError("Incomplete QUAST report")
    return {row[0]: dict(zip(rows[0][1:], row[1:])) for row in rows[1:] if len(row) > 1}


def busco_metrics(fasta, label, folder, lineage, threads):
    """使用离线lineage运行BUSCO genome模式，返回指标字典、摘要路径及full_table路径列表。"""
    execute(
        [
            "busco",
            "-i",
            fasta,
            "-m",
            "genome",
            "--offline",
            "-l",
            lineage,
            "-o",
            label,
            "--out_path",
            folder,
            "--download_path",
            folder / "cache",
            "-c",
            threads,
        ]
    )
    run = folder / label
    summaries = list(run.glob("short_summary.specific*.json"))
    if len(summaries) != 1:
        raise ValueError("BUSCO specific summary is missing or ambiguous")
    data = load(summaries[0])
    values = data.get("results", {})
    required = [
        "Complete percentage",
        "Single copy percentage",
        "Multi copy percentage",
        "Fragmented percentage",
        "Missing percentage",
    ]
    if any(key not in values for key in required):
        raise ValueError("BUSCO lacks C/S/D/F/M percentages")
    result = {symbol: number(values[key]) for symbol, key in zip("CSDFM", required)}
    if any(not 0 <= value <= 100 for value in result.values()):
        raise ValueError("Invalid BUSCO percentage")
    return result, summaries[0], list(run.rglob("full_table.tsv"))


def merqury_metrics(genome, sets, reads, work, threads, memory, kmer):
    """运行并解析Merqury的QV、完整率等结果，返回各序列集合的指标及原生文件位置。"""
    records = run_merqury(genome, sets, reads, work, threads, memory, kmer)
    for record in records.values():
        lines = [
            line
            for line in record["completeness_record"].splitlines()
            if line.strip() and not line.startswith("#")
        ]
        if not lines:
            raise ValueError("Merqury completeness is empty")
        record["completeness_percent"] = number(lines[-1].split()[-1])
        if not 0 <= record["completeness_percent"] <= 100:
            raise ValueError("Invalid Merqury completeness")
    return records


def interval_length(intervals):
    """计算半开区间并集的总长度，拒绝非法区间，避免重叠LTR命中重复计数。"""
    end = total = 0
    for start, stop in sorted(intervals):
        if start < 0 or stop < start:
            raise ValueError("Invalid annotation interval")
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def ltr_fractions(intact, annotation, lengths):
    """解析完整LTR坐标和RepeatMasker注释，返回完整LTR及总LTR覆盖比例。

    按基因组总长度计算，重叠区间仅计一次；完整LTR坐标允许反向记录，但必须位于序列内。
    """
    intact_hits, total_hits = {}, {}
    for line in intact.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        match = re.fullmatch(r"(.+):(\d+)\.\.(\d+)", line.split()[0])
        if not match or match[1] not in lengths:
            raise ValueError("Unrecognized LTR_retriever intact-element coordinate")
        start, stop = sorted((int(match[2]), int(match[3])))
        if start < 1 or stop > lengths[match[1]]:
            raise ValueError("LTR coordinate lies outside its sequence")
        intact_hits.setdefault(match[1], []).append((start - 1, stop))
    for line in annotation.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        if len(fields) < 7 or fields[4] not in lengths:
            raise ValueError("Invalid LTR RepeatMasker annotation")
        if len(fields) < 11 or "LTR" not in fields[10] or int(fields[6]) - int(fields[5]) + 1 < 80:
            continue
        total_hits.setdefault(fields[4], []).append((int(fields[5]) - 1, int(fields[6])))
    total = sum(lengths.values())
    return {
        "intact_fraction": sum(interval_length(hits) for hits in intact_hits.values()) / total,
        "total_ltr_fraction": sum(interval_length(hits) for hits in total_hits.values()) / total,
    }


def lai_metrics(_fasta, sequences, work, threads):
    """运行LTR候选识别、LTR_retriever和LAI，返回状态/指标字典及LAI表路径。

    原生工具使用短ID，交付序列名不改动；低于论文适用范围时明确返回不适用，
    兼容genome.fa和genome.fa.mod结果前缀，不将缺失结果默认为零。
    """
    work.mkdir(parents=True)
    genome = work / "genome.fa"
    identifiers = {f"seq{i:06d}": name for i, name in enumerate(sequences, 1)}
    from omics_workflows.assembly.lib.postassembly import write_sequences

    write_sequences(genome, {short: sequences[name] for short, name in identifiers.items()})
    execute(
        [
            "gt",
            "suffixerator",
            "-db",
            genome,
            "-indexname",
            genome,
            "-tis",
            "-suf",
            "-lcp",
            "-des",
            "-ssp",
            "-sds",
            "-dna",
        ]
    )
    harvest = work / "harvest.scn"
    with harvest.open("w", encoding="utf-8") as output:
        execute(
            [
                "gt",
                "ltrharvest",
                "-index",
                genome,
                "-minlenltr",
                "100",
                "-maxlenltr",
                "7000",
                "-mintsd",
                "4",
                "-maxtsd",
                "6",
                "-motif",
                "TGCA",
                "-motifmis",
                "1",
                "-similar",
                "85",
                "-vic",
                "10",
                "-seed",
                "20",
                "-seqids",
                "yes",
            ],
            stdout=output,
        )
    execute(
        [
            "LTR_FINDER_parallel",
            "-seq",
            genome,
            "-threads",
            threads,
            "-harvest_out",
            "-size",
            "1000000",
            "-time",
            "300",
        ],
        cwd=work,
    )
    finder = work / "genome.fa.finder.combine.scn"
    if not finder.is_file():
        raise ValueError("LTR_FINDER_parallel did not produce its candidate table")
    candidates = work / "raw.scn"
    with candidates.open("w", encoding="utf-8") as target:
        for path in (harvest, finder):
            with path.open(encoding="utf-8") as handle:
                shutil.copyfileobj(handle, target)
    if not any(
        line.strip() and not line.startswith("#")
        for line in candidates.read_text(encoding="utf-8").splitlines()
    ):
        return {"status": "NOT_APPLICABLE", "reason": "No LTR candidates; LAI applicability not met"}, None
    execute(["LTR_retriever", "-genome", genome, "-inharvest", candidates, "-threads", threads], cwd=work)
    prefixes = [work / name for name in ("genome.fa", "genome.fa.mod")]
    completed = [prefix for prefix in prefixes if Path(str(prefix) + ".pass.list").is_file()]
    if len(completed) != 1:
        raise ValueError("LTR_retriever intact-element output is missing or ambiguous")
    prefix = completed[0]
    intact, annotation = Path(str(prefix) + ".pass.list"), Path(str(prefix) + ".out")
    if intact.is_file() and not any(
        line.strip() and not line.startswith("#") for line in intact.read_text(encoding="utf-8").splitlines()
    ):
        return {"status": "NOT_APPLICABLE", "reason": "No intact LTR elements after filtering"}, None
    fractions = ltr_fractions(
        intact, annotation, {short: len(sequences[name]) for short, name in identifiers.items()}
    )
    if fractions["intact_fraction"] < 0.001 or fractions["total_ltr_fraction"] < 0.05:
        return {
            "status": "NOT_APPLICABLE",
            "reason": "LTR content below LAI paper applicability",
            **fractions,
        }, None
    source = Path(str(prefix) + ".out.LAI")
    lines = source.read_text(encoding="utf-8").splitlines()
    global_rows = [line.split() for line in lines if line.startswith("whole_genome")]
    if len(global_rows) != 1:
        raise ValueError("LAI whole-genome row missing or ambiguous")
    output = work / "lai.tsv"
    with output.open("w", encoding="utf-8") as target:
        for line in lines:
            fields = line.split("\t")
            if fields and fields[0] in identifiers:
                fields[0] = identifiers[fields[0]]
            target.write("\t".join(fields) + "\n")
    return {"status": "COMPLETED", "LAI": number(global_rows[0][-1]), **fractions}, output


def mapping_metrics(reference, reads, chromosomes, folder, threads, min_mapq, window):
    """用HiFi回贴及mosdepth计算比对、分裂比对、覆盖广度和深度，返回指标字典。

    使用全部竞争序列进行比对，按指定MAPQ和窗口统计目标染色体，避免漏掉未挂载竞争序列。
    """
    folder.mkdir(parents=True)
    bam = folder / "hifi.bam"
    map_reads(reference, reads, bam, threads)
    primary = mapped = split = 0
    with pysam.AlignmentFile(str(bam)) as alignments:
        for read in alignments.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            primary += 1
            if not read.is_unmapped:
                mapped += 1
                split += read.has_tag("SA")
    if not primary:
        raise ValueError("HiFi mapping contains no primary read records")
    prefix = folder / "coverage"
    execute(
        [
            "mosdepth",
            "-t",
            max(1, threads - 1),
            "-Q",
            min_mapq,
            "--by",
            window,
            "--thresholds",
            "1",
            "--no-per-base",
            prefix,
            bam,
        ]
    )
    covered = {name: 0 for name in chromosomes}
    with gzip.open(folder / "coverage.thresholds.bed.gz", "rt", encoding="utf-8") as handle:
        header = next(handle, "").split()
        if header[:3] != ["#chrom", "start", "end"] or "1X" not in header:
            raise ValueError("Unexpected mosdepth threshold table")
        covered_column = header.index("1X")
        for line in handle:
            fields = line.split()
            if len(fields) != len(header):
                raise ValueError("Malformed mosdepth threshold row")
            if fields[0] in covered:
                count = int(fields[covered_column])
                if not 0 <= count <= int(fields[2]) - int(fields[1]):
                    raise ValueError("Invalid mosdepth covered base count")
                covered[fields[0]] += count
    depth = {}
    with (folder / "coverage.mosdepth.summary.txt").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chrom"] in chromosomes:
                depth[row["chrom"]] = number(row["mean"])
    if set(depth) != set(chromosomes):
        raise ValueError("mosdepth summary does not cover every chromosome")
    write_table(
        folder / "chromosome_coverage.tsv",
        ["sequence_id", "length", "covered_bases", "covered_fraction", "mean_depth"],
        [
            {
                "sequence_id": name,
                "length": len(sequence),
                "covered_bases": covered[name],
                "covered_fraction": covered[name] / len(sequence),
                "mean_depth": depth[name],
            }
            for name, sequence in chromosomes.items()
        ],
    )
    return {
        "status": "COMPLETED",
        "primary_reads": primary,
        "mapped_reads": mapped,
        "mapping_percent": mapped * 100 / primary,
        "split_read_percent": split * 100 / mapped if mapped else 0,
        "covered_fraction": sum(covered.values()) / sum(map(len, chromosomes.values())),
        "coverage_mapq": min_mapq,
        "reference_scope": "all chromosomes plus unplaced competitors",
    }


def telomere_metrics(genome, folder, motif, window, minimum, offset):
    """运行tidk ends并读取原生逐末端结果，返回记录、JSON摘要和TSV路径，不自行重定义状态。"""
    folder.mkdir(parents=True)
    execute(
        [
            "tidk",
            "ends",
            "--string",
            motif,
            "--window",
            window,
            "--min-length",
            minimum,
            "--max-offset",
            offset,
            "--output",
            "telomeres",
            "--dir",
            folder,
            genome,
        ]
    )
    table = folder / "telomeres.ends.tsv"
    if not table.is_file() or not (folder / "telomeres.ends.json").is_file():
        raise ValueError("tidk ends did not produce its native tables")
    with table.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("tidk ends result is empty")
    return rows, load(folder / "telomeres.ends.json"), table
