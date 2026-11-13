#!/usr/bin/env python3
"""Create a small synthetic two-sample project outside the workflow repository."""
import csv
import gzip
import hashlib
import io
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.resources import package_directory, source_boundary


def fastq(path, records):
    # Sequence/read sampling and gzip timestamps are deterministic.
    """把名称/序列记录写为确定性gzip FASTQ，固定压缩时间并为合成读段写统一质量字符。"""
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii") as out:
                for name, sequence in records:
                    out.write(f"@{name}\n{sequence}\n+\n" + "I" * len(sequence) + "\n")


def main():
    """在新的仓库外目录生成可复现的两样本合成项目及配置；仅用于演示和功能验证，不代表真实材料。"""
    parser = argument_parser(description="生成可复现的小型合成数据项目，仅用于熟悉流程和软件功能验证。")
    parser.add_argument(
        "--project-dir",
        required=True,
        type=Path,
        help="类型：目录路径；必填。新建合成案例项目，须位于代码目录之外且尚不存在；已有目录不会被覆盖。",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20261030,
        help="类型：整数；默认20261030。合成序列和读段抽样的随机种子，用于复现示例数据。",
    )
    args = parser.parse_args()
    project = args.project_dir.resolve()
    with package_directory() as package:
        boundary = source_boundary(package)
    if project == boundary or project.is_relative_to(boundary):
        parser.error("Synthetic data and results must be outside the code repository")
    if project.exists():
        parser.error("Project directory already exists; choose a new one (nothing overwritten)")
    samtools = Path(sys.prefix) / "bin/samtools"
    if not samtools.is_file():
        parser.error("Run with the unified environment Python; samtools is required")
    project.mkdir(parents=True)
    data, config_dir = project / "data", project / "config"
    data.mkdir()
    with package_directory() as package:
        shutil.copytree(package / "assembly/examples/hifi_only", config_dir)
    rng = random.Random(args.seed)
    for sid in ["sample_A", "sample_B"]:
        reference = "".join(rng.choices("ACGT", k=100000))
        (data / f"{sid}.synthetic_reference.fa").write_text(
            f">{sid}_synthetic\n{reference}\n", encoding="utf-8"
        )
        records = []
        for i in range(400):
            start = rng.randrange(0, 90001)
            records.append((f"{sid}_read_{i}", reference[start : start + 10000]))
        if sid == "sample_A":
            sam = data / "sample_A.cell01.sam"
            with sam.open("w", encoding="utf-8") as out:
                out.write("@HD\tVN:1.6\tSO:unknown\n")
                for name, sequence in records[:200]:
                    out.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t" + "I" * len(sequence) + "\n")
            subprocess.run(
                [str(samtools), "view", "-b", "-o", str(data / "sample_A.cell01.bam"), str(sam)], check=True
            )
            fastq(data / "sample_A.cell02.fastq.gz", records[200:])
        else:
            fastq(data / "sample_B.cell01.fastq.gz", records)
    table = config_dir / "readsets.tsv"
    with table.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["md5_1"] = hashlib.md5((config_dir / row["file1"]).read_bytes()).hexdigest()
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*rows[0].keys()], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    config_path = config_dir / "config.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["analysis_profiles"]["hifi_contigs"] = {
        "assembly": {"threads": 2, "mem_mb": 2048, "bloom_filter_bits": 0, "purge_level": 0}
    }
    for sid in ["sample_A", "sample_B"]:
        config["samples"][sid]["ploidy"] = 1
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    (config_dir / "samples.tsv").write_text(
        "sample_id\tspecies\tprofile\n"
        "sample_A\tSynthetic sequence A\thifi_contigs\n"
        "sample_B\tSynthetic sequence B\thifi_contigs\n",
        encoding="utf-8",
    )
    (project / "DEMO_ONLY.json").write_text(
        json.dumps(
            {
                "purpose": "Execution example only; not biological or T2T validation",
                "seed": args.seed,
                "samples": ["sample_A", "sample_B"],
                "reference_bases_per_sample": 100000,
                "reads_per_sample": 400,
                "bases_per_read": 10000,
                "profile": "small synthetic data only",
                "note": "BAM bytes may differ across samtools versions and project paths",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"SYNTHETIC_DEMO_PROJECT={project}")
    print("Config: " + str(config_path))
    print("Synthetic-only profile: 2 threads / 2048 MB per assembly; NOT for real genomes")


if __name__ == "__main__":
    main()
