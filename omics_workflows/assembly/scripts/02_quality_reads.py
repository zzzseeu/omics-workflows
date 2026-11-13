#!/usr/bin/env python3
"""Apply the agreed technology-specific read QC and publish one clean-read manifest."""
import shutil

# Direct script execution also works outside the repository.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.assembly.lib.io import (
    assembly_scratch,
    checksum,
    execute,
    load,
    pipeline,
    publish,
    save,
)
from omics_workflows.assembly.lib.postassembly import (
    artifact,
    fastq_statistics,
    read_bundle,
    write_table,
)


def filter_paired(record, _work, output, reports, options, threads):
    """对一组WGS或Hi-C双端读段运行fastp，返回清洁R1/R2路径，并写HTML/JSON质量报告。"""
    files = [output / f'{record["readset_id"]}.R{mate}.fastq.gz' for mate in (1, 2)]
    execute(
        [
            "fastp",
            "-i",
            record["paths"][0],
            "-I",
            record["paths"][1],
            "-o",
            files[0],
            "-O",
            files[1],
            "-w",
            min(threads, 16),
            "--length_required",
            options["fastp_min_length"],
            "--qualified_quality_phred",
            options["fastp_quality"],
            "--json",
            reports / "fastp.json",
            "--html",
            reports / "fastp.html",
        ]
    )
    return files


def filter_long(record, work, output, reports, options, threads):
    """按技术生成NanoPlot报告并用HiFiAdapterFilt或NanoFilt过滤，返回清洁长读段路径。

    HiFi接头过滤使用工作副本，防止上游脚本解压时修改原始输入；清理非必要报告中间文件。
    """
    source = Path(record["paths"][0])
    target = output / f'{record["readset_id"]}.fastq.gz'
    execute(
        [
            "NanoPlot",
            "--fastq",
            source,
            "--outdir",
            reports,
            "--threads",
            threads,
            "--no_static",
            "--plots",
            "dot",
            "--tsv_stats",
        ]
    )
    if record["technology"] == "hifi":
        native = work / "adapter"
        native.mkdir()
        # The upstream shell script decompresses its inputs. Never give it a
        # link to an immutable producer or a user-owned raw input.
        shutil.copyfile(source, native / "reads.fastq.gz")
        execute(
            [
                "hifiadapterfilt.sh",
                "-p",
                "reads",
                "-t",
                threads,
                "-l",
                options["hifi_adapter_length"],
                "-m",
                options["hifi_adapter_identity"],
                "-o",
                "filtered",
            ],
            cwd=native,
        )
        candidates = list(native.rglob("*.filt.fastq.gz"))
        if len(candidates) != 1:
            raise ValueError("HiFiAdapterFilt did not produce exactly one filtered FASTQ")
        shutil.move(candidates[0], target)
        for stat in native.rglob("*.stats"):
            shutil.copyfile(stat, reports / "adapter.stats")
    else:
        pipeline(
            [
                ["gzip", "-dc", source],
                ["NanoFilt", "-q", options["ont_min_quality"], "-l", options["ont_min_length"]],
                ["pigz", "-p", max(1, threads - 2), "-c"],
            ],
            target,
        )
    # HTML reports embed their figures; retain HTML assets and summary only.
    for path in reports.iterdir():
        if path.is_file() and path.suffix not in (".html", ".txt", ".tsv", ".stats"):
            path.unlink()
    return [target]


def quality_sample(args):
    """读取步骤01清单并按技术进行质控，发布02_read_qc的清洁读段、报告、计数和身份清单。"""
    metadata = load(args.metadata)
    prepared = read_bundle(args.inputs_manifest, metadata["sample_id"])
    options = metadata["settings"]["read_qc"]
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        output = stage / "reads"
        output.mkdir(parents=True)
        records = {}
        for rid, record in prepared["readsets"].items():
            work = Path(temporary) / rid
            work.mkdir()
            reports = stage / "reports" / rid
            reports.mkdir(parents=True)
            operation = filter_paired if record["technology"] in ("hic", "wgs") else filter_long
            files = operation(record, work, output, reports, options, args.threads)
            statistics = fastq_statistics(files)
            records[rid] = {key: record[key] for key in ("readset_id", "library_id", "technology")}
            records[rid].update(
                files=[artifact(path, stage) for path in files],
                statistics=statistics,
                input_statistics=record["statistics"],
                filtering=(
                    "fastp"
                    if len(files) == 2
                    else "HiFiAdapterFilt" if record["technology"] == "hifi" else "NanoFilt"
                ),
            )
        write_table(
            stage / "reads.tsv",
            [
                "readset_id",
                "technology",
                "filtering",
                "num_reads",
                "num_bases",
                "mean_length",
                "read_n50",
                "num_pairs",
            ],
            [{**record, **record["statistics"]} for record in records.values()],
        )
        save(
            stage / "manifest.json",
            {
                "sample_id": metadata["sample_id"],
                "status": "PASS",
                "readsets": records,
                "parameters": options,
                "input_manifest_sha256": checksum(args.inputs_manifest),
                "tools": metadata.get("tools", {}),
            },
        )
        publish(stage, destination, temporary)


def main():
    """解析步骤02参数并对一个样本执行约定的测序数据质控路线。"""
    parser = argument_parser(description="步骤02：按测序技术进行质量检查、接头和低质量读段过滤。")
    parser.add_argument(
        "--inputs-manifest", required=True, help="类型：JSON文件路径；必填。步骤01准备后的读段清单。"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="类型：JSON 文件路径；必填。流程生成的样本参数、软件来源和分析设置。",
    )
    parser.add_argument(
        "--out-dir", required=True, help="类型：目录路径；必填。保存本步骤的核心结果、统计表和结果清单。"
    )
    parser.add_argument("--threads", type=int, default=8, help="类型：整数；默认8。分配给本步骤的CPU线程数。")
    quality_sample(parser.parse_args())


if __name__ == "__main__":
    main()
