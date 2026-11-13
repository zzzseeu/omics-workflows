#!/usr/bin/env python3
"""Assemble one reference, diploid haplotypes or a polyploid graph and scaffold it."""
import hashlib
import json
import shutil

# Direct script execution also works outside the repository.
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
# Imports below follow the standalone CLI path bootstrap.
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.workflows.assembly.lib.assembly import (
    estimate_unitig_copies,
    graph_assemble,
    joint_assemble,
    primary_assemble,
    run_gphase,
)
from omics_workflows.workflows.assembly.lib.hic import map_contacts, pair_records, plot_contacts
from omics_workflows.workflows.assembly.lib.io import assembly_scratch, checksum, execute, load, publish, save
from omics_workflows.workflows.assembly.lib.postassembly import artifact, read_bundle, retain_files
from omics_workflows.workflows.assembly.lib.scaffolding import deliver, export_polyploid, scaffold


def contigs(args, metadata, reads):
    """根据single、diploid或polyploid策略调用初步组装实现，保留核心序列/图并发布统一清单。"""
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    strategy = metadata["settings"]["assembly"]["strategy"]
    hifi = [
        p for record in reads["readsets"].values() if record["technology"] == "hifi" for p in record["paths"]
    ]
    ul = [
        p for record in reads["readsets"].values() if record["technology"] == "ont" for p in record["paths"]
    ]
    with assembly_scratch(destination) as temporary:
        scratch = Path(temporary)
        stage = scratch / "result"
        qc = scratch / "qc.json"
        save(qc, {"status": "PASS"})
        paired = pair_records(reads, scratch) if metadata["settings"]["scaffolding"]["enabled"] else []
        if strategy == "single":
            spec = scratch / "spec.json"
            save(
                spec,
                {
                    "sample_id": metadata["sample_id"],
                    "assembly": metadata["settings"]["assembly"],
                    "hifiasm": metadata["tools"].get("hifiasm"),
                    "code": metadata.get("code", {}),
                },
            )
            primary_assemble(
                SimpleNamespace(
                    spec=spec, fastqs=hifi, qcs=[qc], manifests=[], out_dir=stage, threads=args.threads, ul=ul
                )
            )
            (stage / "combined.fa").symlink_to("primary.fa")
            with (stage / "contig_haplotypes.tsv").open("w", encoding="utf-8") as out:
                out.write("contig_id\thaplotype\toriginal_id\n")
                with (stage / "primary.fa").open(encoding="utf-8") as source:
                    for line in source:
                        if line.startswith(">"):
                            name = line[1:].split()[0]
                            out.write(f"{name}\treference\t{name}\n")
            record = load(stage / "manifest.json")
            record.update(combined_sha256=checksum(stage / "combined.fa"), phase_level="UNPHASED_REFERENCE")
            keep = {
                "primary.fa",
                "primary.fa.fai",
                "primary.gfa",
                "alternate.fa",
                "alternate.fa.fai",
                "alternate.gfa",
                "combined.fa",
                "contig_haplotypes.tsv",
            }
        elif strategy == "diploid":
            joint_assemble(
                SimpleNamespace(
                    metadata=args.metadata,
                    out_dir=stage,
                    threads=args.threads,
                    hifi=hifi,
                    ul=ul,
                    paired=paired,
                    qcs=[qc],
                    manifests=[],
                )
            )
            record = load(stage / "manifest.json")
            keep = {
                "hap1.fa",
                "hap2.fa",
                "hap1.fa.fai",
                "hap2.fa.fai",
                "hap1.gfa",
                "hap2.gfa",
                "combined.fa",
                "combined.fa.fai",
                "contig_haplotypes.tsv",
            }
        else:
            graph_assemble(
                SimpleNamespace(
                    metadata=args.metadata,
                    audit=qc,
                    qcs=[qc],
                    hifi=hifi,
                    ul=ul,
                    manifests=[],
                    out_dir=stage,
                    threads=args.threads,
                )
            )
            record = load(stage / "manifest.json")
            keep = {"unitigs.fa", "unitigs.fa.fai", "unitigs.gfa"}
        record.update(
            strategy=strategy,
            ploidy=metadata["settings"]["ploidy"],
            output_copies=metadata["settings"]["output_copies"],
            input_manifest_sha256=checksum(args.reads_manifest),
        )
        record["mode"] = metadata["settings"]["assembly"]["mode"]
        record["assembly_id"] = hashlib.sha256(
            json.dumps(
                {
                    "reads": record["input_manifest_sha256"],
                    "settings": metadata["settings"],
                    "tools": metadata.get("tools", {}),
                    "code": metadata.get("code", {}),
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()
        for fasta in stage.glob("*.fa"):
            if fasta.stat().st_size:
                execute(["samtools", "faidx", fasta])
            else:
                Path(str(fasta) + ".fai").write_text("", encoding="utf-8")
        retain_files(stage, keep | {"manifest.json", ".omics-managed"})
        record["files"] = [
            artifact(path, stage)
            for path in sorted(stage.iterdir())
            if path.name not in ("manifest.json", ".omics-managed")
        ]
        save(stage / "manifest.json", record)
        publish(stage, destination, temporary)


def chromosomes(args, metadata, reads):
    """把初步组装接入HapHiC或GPhase完成染色体交付，并核验输出份数、序列和统计。"""
    destination = Path(args.out_dir).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    assembly = Path(args.contigs_manifest).parent.resolve()
    strategy = metadata["settings"]["assembly"]["strategy"]
    with assembly_scratch(destination) as temporary:
        work = Path(temporary)
        stage = work / "result"
        paired = pair_records(reads, work)
        options = metadata["settings"]["scaffolding"]
        if strategy == "polyploid":
            reference = assembly / "unitigs.fa"
            bam, _ = map_contacts([reference], paired, args.metadata, work / "hic", args.threads, raw=True)
            coverage = work / "coverage"
            hifi = [
                p
                for record in reads["readsets"].values()
                if record["technology"] == "hifi"
                for p in record["paths"]
            ]
            estimate_unitig_copies(
                SimpleNamespace(fasta=reference, hifi=hifi, out_dir=coverage, threads=args.threads)
            )
            native = work / "gphase"
            run_gphase(
                SimpleNamespace(
                    metadata=args.metadata,
                    assembly=args.contigs_manifest,
                    fasta=reference,
                    gfa=assembly / "unitigs.gfa",
                    copies=coverage / "popcnv/06.genes.round.cn",
                    bam=bam,
                    out_dir=native,
                    threads=args.threads,
                )
            )
            export_polyploid(
                SimpleNamespace(
                    metadata=args.metadata,
                    assembly=args.contigs_manifest,
                    raw_unitigs=reference,
                    gphase=native,
                    out_dir=stage,
                )
            )
            shutil.copyfile(coverage / "popcnv/06.genes.round.cn", stage / "unitig_copy_number.tsv")
        else:
            reference = assembly / "combined.fa"
            references = (
                [assembly / f"hap{i}.fa" for i in (1, 2)]
                if options["hic_reference"] == "per_haplotype"
                else [reference]
            )
            bam, linkage = map_contacts(references, paired, args.metadata, work / "hic", args.threads)
            linkage_path = work / "linkage.json"
            save(linkage_path, linkage)
            native = work / "haphic"
            gfas = [assembly / f"hap{i}.gfa" for i in (1, 2)] if strategy == "diploid" else []
            scaffold(
                SimpleNamespace(
                    metadata=args.metadata,
                    fasta=reference,
                    bam=bam,
                    linkage=linkage_path,
                    out_dir=native,
                    threads=args.threads,
                    gfas=gfas,
                )
            )
            deliver(
                SimpleNamespace(
                    metadata=args.metadata,
                    assembly=args.contigs_manifest,
                    fasta=reference,
                    mapping=assembly / "contig_haplotypes.tsv",
                    haphic=native,
                    out_dir=stage,
                )
            )
        record = load(stage / "manifest.json")
        if record["quality_level"] == "INCOMPLETE_SCAFFOLDING":
            raise ValueError("Chromosome groups are incomplete; review native scaffolding output")
        # Re-map even graph-merged components so the plotted AGP and BAM agree.
        statistics = plot_contacts(
            stage / "components.fa",
            stage / "assembly.agp",
            reads,
            args.metadata,
            stage / "contact_maps",
            args.threads,
        )
        retain_files(stage / "contact_maps", {"contact_map.pdf", "separate_plots.pdf"})
        record.update(
            strategy=strategy,
            output_copies=metadata["settings"]["output_copies"],
            hic_mapping=statistics,
            input_reads_manifest_sha256=checksum(args.reads_manifest),
        )
        record["artifacts"].pop("chromosome_agp", None)
        record["artifacts"].pop("review_request", None)
        save(stage / "manifest.json", record)
        keep = {
            "manifest.json",
            "assembly.agp",
            "components.fa",
            "components.fa.fai",
            "chromosome_names.tsv",
            "contact_maps",
            "SHA256SUMS",
            ".omics-managed",
            "unitig_copy_number.tsv",
            "unitig_assembly.agp",
            "contig_to_unitig.txt",
            "unitig_components.fa",
            "unitig_components.fa.fai",
        }
        keep.update(label + suffix for label in record["sets"] for suffix in (".fa", ".fa.fai"))
        retain_files(stage, keep)
        from omics_workflows.workflows.assembly.lib.postassembly import write_checksums

        write_checksums(stage)
        publish(stage, destination, temporary)


def main():
    """解析步骤03阶段参数，核验读段及样本元数据后分派contigs或chromosomes分析。"""
    parser = argument_parser(description="步骤03：初步组装，以及按策略进行染色体分相与挂载。")
    parser.add_argument(
        "--metadata",
        required=True,
        help="类型：JSON 文件路径；必填。流程生成的样本参数、软件来源和分析设置。",
    )
    parser.add_argument(
        "--reads-manifest",
        required=True,
        help="类型：JSON 文件路径；必填。步骤02质控后的读段清单，含样本身份、文件位置和校验值。",
    )
    parser.add_argument(
        "--contigs-manifest",
        help="类型：JSON文件路径；执行chromosomes时必填。步骤03初步组装的结果清单；contigs阶段无需提供。",
    )
    parser.add_argument(
        "--stage",
        choices=["contigs", "chromosomes"],
        required=True,
        help="类型：字符串枚举；必填。contigs生成初步组装或原始unitig图；chromosomes进行染色体分相与挂载。",
    )
    parser.add_argument(
        "--out-dir", required=True, help="类型：目录路径；必填。保存本步骤的核心结果、统计表和结果清单。"
    )
    parser.add_argument("--threads", type=int, default=8, help="类型：整数；默认8。分配给本步骤的CPU线程数。")
    args = parser.parse_args()
    metadata = load(args.metadata)
    reads = read_bundle(args.reads_manifest, metadata["sample_id"])
    if args.stage == "chromosomes" and not args.contigs_manifest:
        parser.error("--contigs-manifest is required for chromosomes")
    (contigs if args.stage == "contigs" else chromosomes)(args, metadata, reads)


if __name__ == "__main__":
    main()
