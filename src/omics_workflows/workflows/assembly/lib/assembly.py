"""Single-reference, diploid and polyploid assembly route implementations."""

import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

from omics_workflows.workflows.assembly.lib.io import (
    assembly_scratch,
    checksum,
    execute,
    load,
    pipeline,
    publish,
    save,
    source_check,
)
from omics_workflows.workflows.assembly.lib.sequences import gfa_to_fasta, prefixed_gfa


def primary_assemble(args):
    """使用清洁HiFi及可选超长读段运行hifiasm主组装，写primary/alternate的GFA、FASTA和清单。

    每次重算使用全新缓存；替代序列不作为第二份完整单倍型，拒绝覆盖非流程管理目录。
    """
    spec = load(args.spec)
    options = spec["assembly"]
    qcs = [load(p) for p in args.qcs]
    if any(q["status"] != "PASS" for q in qcs):
        raise ValueError("Read QC did not pass")
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Every recomputation gets a fresh cache; old hifiasm .bin files are never silently reused.
    with tempfile.TemporaryDirectory(prefix=".hifiasm-", dir=destination.parent) as temp:
        stage = Path(temp) / "assembly"
        stage.mkdir()
        prefix = stage / "asm"
        command = [
            "hifiasm",
            "-o",
            str(prefix),
            "-t",
            str(args.threads),
            "--primary",
            "-f",
            str(options["bloom_filter_bits"]),
        ]
        if options["purge_level"] is not None:
            command += ["-l", str(options["purge_level"])]
        if getattr(args, "ul", None):
            command += ["--ul", ",".join(args.ul)]
        command += args.fastqs
        execute(command)
        print("HIFIASM_OUTPUTS", json.dumps(sorted(p.name for p in stage.iterdir())), flush=True)
        sets = []
        for label, suffix in [("primary", "p_ctg"), ("alternate", "a_ctg")]:
            candidates = [stage / f"asm.bp.{suffix}.gfa", stage / f"asm.{suffix}.gfa"]
            existing = [p for p in candidates if p.is_file()]
            if len(existing) > 1:
                raise ValueError(f"Ambiguous hifiasm {label} GFA output: {existing}")
            if label == "primary" and not existing:
                raise ValueError(f"Missing primary GFA; expected one of {candidates}")
            gfa = existing[0] if existing else candidates[0]
            normalized_gfa = stage / f"{label}.gfa"
            if gfa.exists():
                shutil.copyfile(gfa, normalized_gfa)
            else:
                normalized_gfa.write_text("", encoding="utf-8")
            fasta = stage / f"{label}.fa"
            count = gfa_to_fasta(gfa, fasta, label == "primary")
            sets.append(
                {
                    "set_id": label,
                    "representation": label,
                    "num_contigs": count,
                    "has_sequence": bool(count),
                    "fasta": fasta.name,
                    "fai": fasta.name + ".fai",
                    "normalized_gfa": normalized_gfa.name,
                    "source_gfa": gfa.name if gfa.exists() else None,
                }
            )
        inputs = [load(p) for p in args.manifests]
        identity = {"spec": spec, "input_sha256": [p["source_sha256"] for p in inputs]}
        assembly_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(
            stage / "manifest.json",
            {
                "sample_id": spec["sample_id"],
                "assembly_id": assembly_id,
                "mode": "hifi_only",
                "representation": "primary_alternate",
                "phase_level": "NOT_CERTIFIED",
                "sets": sets,
                "command": command,
                "settings": options,
                "input_manifests": inputs,
                "tools": {"hifiasm": spec["hifiasm"]},
                "code": spec["code"],
                "cache_policy": "fresh_on_recompute",
            },
        )
        (stage / ".omics-managed").write_text("omics-workflows assembly output\n", encoding="utf-8")
        if destination.exists():
            # Snakemake creates output parents before executing a file-output rule.
            if any(destination.iterdir()) and not (destination / ".omics-managed").is_file():
                raise ValueError(f"Refusing to replace unmanaged assembly directory: {destination}")
            destination.replace(Path(temp) / "previous")
        stage.replace(destination)


def joint_assemble(args):
    """使用HiFi、Hi-C及可选超长读段运行hifiasm二倍体分相，保留两份单倍型及合并序列和来源。"""
    metadata = load(args.metadata)
    options = metadata["settings"]["assembly"]
    for q in args.qcs:
        if load(q)["status"] != "PASS":
            raise ValueError("Read QC did not pass")
    hic = [load(p) for p in args.paired]
    for item in hic:
        source_check(item["source_spec"])
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temp:
        stage = Path(temp) / "assembly"
        stage.mkdir()
        command = [
            "hifiasm",
            "-o",
            stage / "asm",
            "-t",
            args.threads,
            "-f",
            options["bloom_filter_bits"],
            "--h1",
            ",".join(x["source_spec"]["file1"] for x in hic),
            "--h2",
            ",".join(x["source_spec"]["file2"] for x in hic),
        ]
        if args.ul:
            command += ["--ul", ",".join(args.ul)]
        if options["purge_level"] is not None:
            command += ["-l", options["purge_level"]]
        command += args.hifi
        execute(command)
        sets = []
        with (
            (stage / "combined.fa").open("wb") as combined,
            (stage / "contig_haplotypes.tsv").open("w", encoding="utf-8") as mapping,
        ):
            mapping.write("contig_id\thaplotype\toriginal_id\n")
            for hap in ("hap1", "hap2"):
                source = stage / f"asm.hic.{hap}.p_ctg.gfa"
                if not source.is_file():
                    raise ValueError(f"Missing phased hifiasm output: {source}")
                target = stage / f"{hap}.gfa"
                prefixed_gfa(source, target, hap + "_")
                fasta = stage / f"{hap}.fa"
                count = gfa_to_fasta(target, fasta, True)
                with fasta.open("rb") as f:
                    shutil.copyfileobj(f, combined)
                with fasta.open(encoding="utf-8") as f:
                    for line in f:
                        if line.startswith(">"):
                            cid = line[1:].strip()
                            mapping.write(f"{cid}\t{hap}\t{cid[len(hap)+1:]}\n")
                sets.append(
                    {
                        "set_id": hap,
                        "fasta": fasta.name,
                        "normalized_gfa": target.name,
                        "source_gfa": source.name,
                        "num_contigs": count,
                    }
                )
        execute(["samtools", "faidx", stage / "combined.fa"])
        for item in hic:
            source_check(item["source_spec"])
        inputs = [load(p) for p in args.manifests] + hic
        identity = {
            "options": options,
            "inputs": inputs,
            "tools": metadata["tools"],
            "code": metadata["code"],
        }
        assembly_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(
            stage / "manifest.json",
            {
                "sample_id": metadata["sample_id"],
                "assembly_id": assembly_id,
                "mode": options["mode"],
                "representation": "phased",
                "phase_level": "HIC_PHASED_UNVALIDATED",
                "sets": sets,
                "command": [str(x) for x in command],
                "metadata": metadata,
                "input_manifests": inputs,
                "combined_sha256": checksum(stage / "combined.fa"),
                "limitations": [
                    (
                        "Haplotype labels are not parental identities; chromosome-wide "
                        "phase switches remain possible"
                    )
                ],
            },
        )
        publish(stage, destination, temp)


def graph_assemble(args):
    """运行多倍体路线所需的hifiasm原始primary unitig图，保留原始节点ID、顺序和序列。"""
    metadata = load(args.metadata)
    options = metadata["settings"]["assembly"]
    if any(load(path)["status"] != "PASS" for path in args.qcs):
        raise ValueError("Read QC failed")
    if load(args.audit)["status"] != "PASS":
        raise ValueError("Long-read ID audit failed")
    destination = Path(args.out_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(destination) as temp:
        stage = Path(temp) / "assembly"
        stage.mkdir()
        command = [
            "hifiasm",
            "-o",
            stage / "asm",
            "-t",
            args.threads,
            "-f",
            options["bloom_filter_bits"],
            "--n-hap",
            metadata["settings"]["ploidy"],
            "-l",
            0,
            "--ul",
            ",".join(args.ul),
            *args.hifi,
        ]
        execute(command)
        graphs = [p for p in (stage / "asm.bp.p_utg.gfa", stage / "asm.p_utg.gfa") if p.is_file()]
        if len(graphs) != 1:
            raise ValueError(f"Expected one native primary-unitig graph, found: {graphs}")
        shutil.copyfile(graphs[0], stage / "unitigs.gfa")
        count = gfa_to_fasta(stage / "unitigs.gfa", stage / "unitigs.fa", True)
        identity = {"metadata": metadata, "inputs": [load(p) for p in args.manifests]}
        aid = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        save(
            stage / "manifest.json",
            {
                "sample_id": metadata["sample_id"],
                "assembly_id": aid,
                "mode": "polyploid_hifi_ul_unitigs",
                "representation": "unitigs",
                "phase_level": "UNPHASED_UNITIG_GRAPH",
                "ploidy": metadata["settings"]["ploidy"],
                "command": list(map(str, command)),
                "source_gfa": graphs[0].name,
                "unitig_count": count,
                "gfa_sha256": checksum(stage / "unitigs.gfa"),
                "fasta_sha256": checksum(stage / "unitigs.fa"),
                "unitig_policy": (
                    "Native p_utg IDs/order/sequences unchanged; no contig" " purge, polishing or scaffolding"
                ),
                "metadata": metadata,
            },
        )
        publish(stage, destination, temp)


def mapping_commands(fasta, reads, bam, threads):
    """返回HiFi比对与排序的管道命令列表；多个文库作为未配对长读段流，共享线程预算。"""
    mapper = max(1, (threads - 2) // 2)
    sorter = max(1, threads - mapper - 2)
    inputs = [Path(path).resolve() for path in reads]
    commands = [
        ["minimap2", "-t", mapper, "-I", "16G", "-ax", "map-hifi", fasta, *inputs],
        ["samtools", "sort", "-@", sorter, "-m", "512M", "-o", bam, "-"],
    ]
    if len(inputs) > 1:
        # Separate minimap2 positional files describe paired fragments. HiFi
        # readsets are independent molecules, including across SMRT cells.
        stream = (
            "import sys,gzip,shutil\n"
            "for p in sys.argv[1:]:\n"
            ' with open(p,"rb") as h: magic=h.read(2)\n'
            " opener=gzip.open if magic==bytes([31,139]) else open\n"
            ' with opener(p,"rb") as h: shutil.copyfileobj(h,sys.stdout.buffer,8388608)\n'
        )
        commands[0] = commands[0][: -len(inputs)] + ["-"]
        commands.insert(0, [sys.executable, "-c", stream, *inputs])
    return commands


def estimate_unitig_copies(args):
    """将HiFi回贴到原始unitig，按GPhase的popCNV路线估计拷贝数，写拷贝数与来源清单。"""
    dest = Path(args.out_dir)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fasta = Path(args.fasta).resolve()
    gphase = Path(sys.prefix) / "opt/GPhase/source"
    with assembly_scratch(dest) as temp:
        stage = Path(temp) / "coverage"
        stage.mkdir()
        bam = stage / "reads.bam"
        pipeline(mapping_commands(fasta, args.hifi, bam, args.threads), stage / "mapping.stdout")
        execute(["samtools", "quickcheck", "-v", bam])
        execute(["samtools", "index", "-@", args.threads, bam])
        (stage / "bam_files").mkdir()
        (stage / "bam_files/reads.bam").symlink_to("../reads.bam")
        (stage / "bam_files/reads.bam.bai").symlink_to("../reads.bam.bai")
        (stage / "group.list").write_text("reads\treads\n", encoding="utf-8")
        names = set()
        with (
            Path(str(fasta) + ".fai").open(encoding="utf-8") as src,
            (stage / "gene.list").open("w", encoding="utf-8") as out,
        ):
            for line in src:
                fields = line.split("\t")
                names.add(fields[0])
                out.write(f"{fields[0]}\t1\t{fields[1]}\t{fields[0]}\n")
        command = [
            sys.executable,
            gphase / "src/popCNV/popCNV.py",
            "-g",
            fasta,
            "-s",
            "1000",
            "-b",
            "bam_files",
            "-l",
            "gene.list",
            "-w",
            "popcnv",
            "--group",
            "group.list",
            "--wild",
            "0",
            "-t",
            args.threads,
        ]
        execute(command, cwd=stage)
        cn = stage / "popcnv/06.genes.round.cn"
        observed = set()
        with cn.open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                observed.add(fields[0])
                value = float(fields[1])
                if not (value >= 0 and value.is_integer()):
                    raise ValueError("Invalid or missing unitig copy-number estimate")
        if names != observed:
            raise ValueError("popCNV does not cover every input unitig")
        save(
            stage / "manifest.json",
            {
                "status": "PASS",
                "input_fasta_sha256": checksum(fasta),
                "copy_number_sha256": checksum(cn),
                "unitig_count": len(names),
                "command": list(map(str, command)),
                "limitations": [
                    "Read-depth copy numbers are estimates and do not validate chromosome phasing"
                ],
            },
        )
        publish(stage, dest, temp)


def relocate_internal_links(stage, destination):
    """将暂存目录内的绝对内部链接调整为发布后位置，外部输入链接不改动。"""
    stage = Path(stage)
    destination = Path(destination).resolve()
    for path in stage.rglob("*"):
        if path.is_symlink():
            target = path.readlink()
            if target.is_absolute() and target.is_relative_to(stage):
                new = destination / target.relative_to(stage)
                path.unlink()
                path.symlink_to(new)


def run_gphase(args):
    """核验原始图和拷贝数身份后运行GPhase原生分相与挂载，并发布原生结果及来源记录。

    倍性和线程预算来自样本参数，不复制二倍体结果模拟多倍体；完成后仍需独立相位验证。
    """
    meta = load(args.metadata)
    assembly = load(args.assembly)
    if checksum(args.fasta) != assembly["fasta_sha256"] or checksum(args.gfa) != assembly["gfa_sha256"]:
        raise ValueError("Primary unitig graph was modified after assembly")
    dest = Path(args.out_dir)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with assembly_scratch(dest) as temp:
        stage = Path(temp) / "gphase"
        stage.mkdir()
        command = [
            "gphase",
            "pipeline",
            "-f",
            Path(args.fasta).resolve(),
            "-g",
            Path(args.gfa).resolve(),
            "-c",
            Path(args.copies).resolve(),
            "-m",
            Path(args.bam).resolve(),
            "--n_chr",
            meta["settings"]["chromosomes_per_haplotype"],
            "--n_hap",
            meta["settings"]["ploidy"],
            "-p",
            "assembly",
            "--rescue",
            "--thread",
            max(1, args.threads - 2),
            "--processes",
            min(args.threads, meta["settings"]["scaffolding"]["processes"]),
            "-e",
            meta["settings"]["scaffolding"]["restriction_sites"],
        ]
        execute(command, cwd=stage)
        required = [
            "gphase_final.unitig.fasta",
            "gphase_final.unitig.scaffold.agp",
            "gphase_final.unitig.scaffold.fasta",
            "gphase_final.contig.fasta",
            "gphase_final.contig.scaffold.agp",
            "gphase_final.contig.scaffold.fasta",
            "gphase_final.ctg2utg.txt",
        ]
        for name in required:
            path = stage / "gphase_output/scaffold_hap" / name
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing GPhase output: {name}")
        if checksum(args.fasta) != assembly["fasta_sha256"] or checksum(args.gfa) != assembly["gfa_sha256"]:
            raise ValueError("GPhase modified its raw unitig inputs")
        save(
            stage / "manifest.json",
            {
                "status": "AUTOMATIC_POLYPLOID_CANDIDATE",
                "assembly_id": assembly["assembly_id"],
                "command": list(map(str, command)),
                "input_fasta_sha256": checksum(args.fasta),
                "input_gfa_sha256": checksum(args.gfa),
                "copy_number_sha256": checksum(args.copies),
                "human_review": "NOT_REVIEWED",
                "T2T_status": "NOT_ASSESSED",
            },
        )
        # GPhase makes absolute links between its stages. Preserve their meaning
        # when publishing the complete working directory from a temporary path.
        relocate_internal_links(stage, dest)
        publish(stage, dest, temp)
