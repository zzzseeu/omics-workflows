"""FASTA, GFA and AGP parsing and chromosome delivery indexing."""

from collections import OrderedDict, defaultdict
from pathlib import Path

from omics_workflows.assembly.lib.io import checksum, execute

COMPLEMENT = str.maketrans("ACGTRYMKBDHVNSW", "TGCAYRKMVHDBNSW")


def gfa_to_fasta(gfa, fasta, required):
    """提取GFA的S记录为FASTA并返回序列数；检查重复ID，required决定空结果是否报错。"""
    count, seen = 0, set()
    with Path(fasta).open("w", encoding="utf-8") as out:
        if Path(gfa).exists():
            with Path(gfa).open(encoding="utf-8") as handle:
                for line in handle:
                    if not line.startswith("S\t"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3 or fields[2] == "*" or fields[1] in seen:
                        raise ValueError(f"Invalid sequence-containing GFA segment in {gfa}")
                    seen.add(fields[1])
                    out.write(f">{fields[1]}\n{fields[2]}\n")
                    count += 1
    if required and not count:
        raise ValueError("No primary contigs produced by hifiasm")
    return count


def fasta_statistics(path):
    """流式读取FASTA并检查ID和序列有效性，返回长度、N50/L50、GC、N碱基数及文件指纹。"""
    lengths, seen, name, length, n_count, total_gc = [], set(), None, 0, 0, 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if not length:
                        raise ValueError("Empty FASTA sequence")
                    lengths.append(length)
                name = line[1:].split()[0]
                if name in seen:
                    raise ValueError("Duplicate contig ID")
                seen.add(name)
                length = 0
            else:
                seq = line.strip().upper()
                if seq and name is None:
                    raise ValueError("FASTA sequence without a header")
                if any(base not in "ACGTRYSWKMBDHVN" for base in seq):
                    raise ValueError("Unexpected FASTA alphabet")
                length += len(seq)
                n_count += seq.count("N")
                total_gc += seq.count("G") + seq.count("C")
    if name is not None:
        if not length:
            raise ValueError("Empty FASTA sequence")
        lengths.append(length)
    total, cumulative, n50, l50 = sum(lengths), 0, 0, 0
    for i, value in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += value
        if cumulative >= total / 2:
            n50, l50 = value, i
            break
    return {
        "num_contigs": len(lengths),
        "total_bases": total,
        "N50": n50,
        "L50": l50,
        "longest_contig": max(lengths, default=0),
        "N_bases": n_count,
        "GC_percent": 100 * total_gc / total if total else None,
        "sha256": checksum(path),
    }


def prefixed_gfa(source, target, prefix):
    """给hifiasm GFA的片段及链接引用同步添加前缀，保持图连接关系，拒绝不支持的记录类型。"""
    with Path(source).open(encoding="utf-8") as src, Path(target).open("w", encoding="utf-8") as out:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in {"S", "A", "L"}:
                fields[1] = prefix + fields[1]
                if fields[0] == "L":
                    fields[3] = prefix + fields[3]
            elif fields[0] != "H":
                raise ValueError(f"Unsupported hifiasm GFA record: {fields[0]}")
            out.write("\t".join(fields) + "\n")


def fasta_dict(path):
    """读取FASTA为保持原顺序的大写序列字典；重复ID或空序列时报错。"""
    import pysam

    result = OrderedDict()
    with pysam.FastxFile(str(path)) as f:
        for r in f:
            if r.name in result or not r.sequence:
                raise ValueError("Duplicate ID or empty FASTA sequence")
            result[r.name] = r.sequence.upper()
    return result


def write_fasta(handle, name, sequence):
    """向已打开的文本句柄写入一条FASTA记录，每行60个碱基，不负责关闭句柄。"""
    handle.write(f">{name}\n")
    for i in range(0, len(sequence), 60):
        handle.write(sequence[i : i + 60] + "\n")


def parse_agp(path, sequences, require_complete=True):
    """解析并返回按目标序列分组的AGP记录，核对坐标连续性、片段边界、方向及源序列覆盖。

    require_complete为真时要求源序列完整且不重复使用；为假时允许未覆盖部分，仍拒绝重叠使用。
    """
    objects, coverage = OrderedDict(), defaultdict(list)
    closed_objects, current = set(), None
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip().split("\t")
            if len(c) != 9:
                raise ValueError("AGP requires nine tab-separated fields")
            obj, start, end, part, kind = c[:5]
            start, end, part = int(start), int(end), int(part)
            if obj != current:
                if obj in closed_objects:
                    raise ValueError("AGP object rows must be contiguous")
                if current is not None:
                    closed_objects.add(current)
                current = obj
            rows = objects.setdefault(obj, [])
            if start != (int(rows[-1][2]) + 1 if rows else 1) or end < start or part != len(rows) + 1:
                raise ValueError("AGP object coordinates or part numbers do not tile")
            if kind == "W":
                component, a, b, orientation = c[5], int(c[6]), int(c[7]), c[8]
                if component not in sequences or not 1 <= a <= b <= len(sequences[component]):
                    raise ValueError("AGP component out of bounds or absent from input")
                if b - a != end - start or orientation not in {"+", "-"}:
                    raise ValueError("AGP component length/orientation is invalid")
                coverage[component].append((a, b))
            elif kind in {"N", "U"}:
                if int(c[5]) != end - start + 1 or c[7] != "yes" or c[8] != "proximity_ligation":
                    raise ValueError("Invalid or unsupported AGP scaffold gap")
            else:
                raise ValueError(f"Unsupported AGP component type: {kind}")
            rows.append(c)
    # 单套/二倍体交付要求完整源序列覆盖；图导出可放宽完整性，但不允许区间重叠。
    for cid, sequence in sequences.items():
        end = 0
        for a, b in sorted(coverage[cid]):
            if a <= end or (require_complete and a != end + 1):
                raise ValueError(f"AGP duplicates or loses source sequence: {cid}")
            end = b
        if require_complete and end != len(sequence):
            raise ValueError(f"AGP loses source sequence: {cid}")
    return objects


def reconstruct(rows, sequences):
    """根据AGP记录和组成序列字典重建一条序列，按方向取片段并用N表示连接缺口。"""
    parts = []
    for row in rows:
        if row[4] in {"N", "U"}:
            parts.append("N" * int(row[5]))
        else:
            seq = sequences[row[5]][int(row[6]) - 1 : int(row[7])]
            parts.append(seq if row[8] == "+" else seq.translate(COMPLEMENT)[::-1])
    return "".join(parts)


def delivery_labels(ploidy, strategy=None):
    """按交付份数和策略返回染色体全集、各单套/单倍型及未挂载集合的标准文件标签。"""
    copies = (
        ["reference.chromosomes"]
        if strategy == "single"
        else [f"hap{i}.chromosomes" for i in range(1, ploidy + 1)]
    )
    return ["chromosomes", *copies, "unplaced"]


def complete_delivery(stage, result, phasing_tool, scaffolding_tools):
    """补齐不同挂载路线共用的FASTA索引、染色体AGP和交付说明，并更新result字典。

    具体源序列覆盖由各导出器检查：二倍体contig要求精确覆盖，GPhase允许图修剪和折叠拷贝展开。
    此处统一文件接口，不验证亲本相位，也不把输出格式正确解释为生物学验收通过。
    """
    stage = Path(stage)
    labels = delivery_labels(result.get("output_copies", result["ploidy"]), result.get("strategy"))
    if set(result["sets"]) != set(labels):
        raise ValueError("Delivery sequence sets do not match the configured ploidy")
    for label in [*labels, "components"]:
        fasta = stage / f"{label}.fa"
        index = Path(str(fasta) + ".fai")
        if fasta.stat().st_size:
            execute(["samtools", "faidx", fasta])
        else:
            index.write_text("", encoding="utf-8")
        if label in labels:
            with (stage / f"{label}.chrom.sizes").open("w", encoding="utf-8") as out:
                for line in index.read_text(encoding="utf-8").splitlines():
                    out.write("\t".join(line.split("\t")[:2]) + "\n")
    chromosomes = {
        line.split("\t")[0]
        for line in (stage / "chromosomes.fa.fai").read_text(encoding="utf-8").splitlines()
    }
    with (
        (stage / "assembly.agp").open(encoding="utf-8") as source,
        (stage / "chromosomes.agp").open("w", encoding="utf-8") as out,
    ):
        for line in source:
            if line.strip() and not line.startswith("#") and line.split("\t", 1)[0] in chromosomes:
                out.write(line)
    result.update(
        analysis_stage="phasing_scaffolding",
        phasing_tool=phasing_tool,
        scaffolding_tools=scaffolding_tools,
        phase_validation="NOT_ASSESSED",
        artifacts={
            "sequence_sets": {label: f"{label}.fa" for label in labels},
            "components": "components.fa",
            "agp": "assembly.agp",
            "chromosome_agp": "chromosomes.agp",
            "chromosome_table": "chromosome_names.tsv",
            "review_request": "review.request.json",
        },
    )
