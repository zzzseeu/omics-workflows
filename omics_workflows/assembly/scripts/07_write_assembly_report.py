#!/usr/bin/env python3
"""Publish one compact final report bound to the exact evaluated FASTA."""
import csv
import html

# Direct script execution also works outside the repository.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.assembly.lib.io import assembly_scratch, checksum, load, publish, save
from omics_workflows.assembly.lib.postassembly import assembly_bundle, write_table


def html_table(fields, rows):
    """把列名和字典记录转换为转义后的HTML表格字符串，防止样本或指标文本被当作HTML执行。"""
    header = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f'<td>{html.escape(str(row.get(field, "")))}</td>' for field in fields) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def page(title, body):
    """为已构建的报告正文添加UTF-8中文HTML页面和简单样式，返回完整页面字符串。"""
    return (
        '<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
        f"<title>{html.escape(title)}</title><style>body{{font:16px sans-serif;margin:2em;line-height:1.6}}"
        "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:.4em .7em}"
        f"th{{background:#f2f4f7}}</style><h1>{html.escape(title)}</h1>{body}</html>"
    )


def report(args):
    """核对最终组装与评估指纹，生成样本HTML报告和来源清单，不重新分析数据或自动认证T2T。"""
    assessment = load(args.assessment)
    source_path = Path(args.assembly_manifest)
    if assessment["assembly_manifest_sha256"] != checksum(source_path):
        raise ValueError("Assessment belongs to a different assembly manifest")
    if assessment.get("endpoint") == "contigs":
        source = load(source_path)
        identity = source["sample_id"]
    else:
        source, _, _ = assembly_bundle(source_path)
        identity = source["sample_id"]
        if checksum(source["genome_path"]) != assessment["genome_sha256"]:
            raise ValueError("Final assembly changed after quality assessment")
    with (Path(args.assessment).parent / "metrics.tsv").open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    body = (
        f'<p>样本：{html.escape(identity)}；组装策略：{html.escape(source.get("strategy", "single"))}；'
        f'评估状态：{html.escape(assessment["status"])}。</p>'
        "<p>指标按原生工具报告。未评估项目、相位验证和人工审阅单独记录，运行完成不表示 T2T 认证。</p>"
    )
    body += html_table(fields, rows)
    if assessment.get("endpoint") != "contigs":
        body += (
            '<p><a href="../06_assessment/chromosomes.tsv">逐染色体统计</a> · '
            '<a href="../06_assessment/telomeres.tsv">端粒结果</a> · '
            '<a href="../06_assessment/contact_maps/contact_map.pdf">Hi-C 接触图</a></p>'
        )
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        (stage / "report.html").write_text(page(f"{identity} 组装质量报告", body), encoding="utf-8")
        save(
            stage / "manifest.json",
            {
                "sample_id": identity,
                "status": assessment["status"],
                "assembly_manifest": str(source_path.resolve()),
                "assembly_manifest_sha256": checksum(source_path),
                "assessment": str(Path(args.assessment).resolve()),
                "assessment_sha256": checksum(args.assessment),
                "report_sha256": checksum(stage / "report.html"),
                "T2T_status": "NOT_CERTIFIED",
            },
        )
        publish(stage, destination, temporary)


def batch(args):
    """核验各样本报告指纹，生成批次状态TSV和HTML汇总，不把部分评估提升为完整验收。"""
    rows = []
    for path in args.reports:
        record = load(path)
        if checksum(Path(path).parent / "report.html") != record["report_sha256"]:
            raise ValueError("Sample HTML report changed before batch collection")
        rows.append({key: record[key] for key in ("sample_id", "status", "T2T_status")})
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        fields = ["sample_id", "status", "T2T_status"]
        write_table(stage / "summary.tsv", fields, rows)
        (stage / "summary.html").write_text(page("组装批次汇总", html_table(fields, rows)), encoding="utf-8")
        publish(stage, destination, temporary)


def main():
    """解析步骤07参数，按样本清单加评估结果或样本报告列表选择单样本/批次模式。"""
    parser = argument_parser(description="步骤07：生成样本报告或批次汇总，核对结果来源；不重复运行分析。")
    parser.add_argument(
        "--assembly-manifest",
        help="类型：JSON文件路径。样本报告必填；最终组装清单，与--assessment配套；批次模式无需提供。",
    )
    parser.add_argument(
        "--assessment", help="类型：JSON文件路径。样本报告必填；步骤06的assessment.json，须与组装版本一致。"
    )
    parser.add_argument(
        "--reports",
        nargs="+",
        help="类型：JSON文件路径列表，以空格分隔。批次模式必填；各样本步骤07的manifest.json；提供后按批次模式运行。",
    )
    parser.add_argument(
        "--out-dir", required=True, help="类型：目录路径；必填。保存本步骤的核心结果、统计表和结果清单。"
    )
    args = parser.parse_args()
    if args.reports:
        batch(args)
    elif args.assembly_manifest and args.assessment:
        report(args)
    else:
        parser.error("Provide --reports, or --assembly-manifest with --assessment")


if __name__ == "__main__":
    main()
