"""Assembly-specific validation, stage definitions and workflow launcher."""

import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Keep direct invocation available outside the repository, including --help.
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows import __version__
from omics_workflows.common import workflow
from omics_workflows.common.resources import package_directory, source_boundary
from omics_workflows.common.workflow import (
    ConfigError,
    check_fields,
    check_value,
    digest,
    group_readsets,
    inspect_tools,
    load_configuration,
    load_defaults,
    merge,
    prepend_search_paths,
    require_available,
    select_samples,
    source_signature,
    table,
    validate_run_directory,
)

LAYOUT_VERSION = 3

STAGES = {
    "inputs": {
        "rule": "assembly/rules/01_inputs.smk",
        "entries": ["assembly/scripts/01_prepare_inputs.py"],
        "tools": {"seqkit", "samtools"},
    },
    "read_qc": {
        "rule": "assembly/rules/02_read_qc.smk",
        "entries": ["assembly/scripts/02_quality_reads.py"],
        "tools": {"hifiadapterfilt.sh", "NanoFilt", "NanoPlot", "fastp"},
    },
    "contigs": {
        "rule": "assembly/rules/03_assembly.smk",
        "entries": ["assembly/scripts/03_assemble_genome.py"],
        "tools": {"hifiasm", "samtools"},
    },
    "scaffolding": {
        "rule": "assembly/rules/03_assembly.smk",
        "entries": ["assembly/scripts/03_assemble_genome.py"],
        "tools": {
            "filter_bam",
            "bwa",
            "samblaster",
            "gphase",
            "samtools",
            "yahs",
            "minimap2",
            "mosdepth",
            "haphic",
        },
    },
    "finishing": {
        "rule": "assembly/rules/04_finishing.smk",
        "entries": ["assembly/scripts/04_finish_chromosomes.py"],
        "tools": {"quartet", "hifiasm", "samtools", "minimap2", "tidk"},
    },
    "polishing": {
        "rule": "assembly/rules/05_polishing.smk",
        "entries": ["assembly/scripts/05_polish_chromosomes.py"],
        "tools": {"yak", "samtools", "nextPolish2", "minimap2"},
    },
    "assessment": {
        "rule": "assembly/rules/06_assessment.smk",
        "entries": ["assembly/scripts/06_assess_assembly.py"],
        "tools": None,
    },
    "report": {
        "rule": "assembly/rules/07_report.smk",
        "entries": ["assembly/scripts/07_write_assembly_report.py"],
        "tools": set(),
    },
}


def validate_settings(settings: dict, path: str) -> None:
    """在Python中校验一组参数覆盖：允许的字段、实际类型、取值范围和枚举。

    适用于defaults、命名参数组合及单样本覆盖；未填写的值由模块默认参数补齐。
    数据是否存在、倍性与路线是否匹配等组合条件仍由resolve_config集中检查。
    """
    groups = {
        "assembly": {
            "backend",
            "strategy",
            "mode",
            "representation",
            "threads",
            "mem_mb",
            "bloom_filter_bits",
            "purge_level",
        },
        "read_qc": {
            "threads",
            "mem_mb",
            "ont_min_length",
            "ont_min_quality",
            "fastp_min_length",
            "fastp_quality",
            "hifi_adapter_length",
            "hifi_adapter_identity",
        },
        "scaffolding": {
            "enabled",
            "backend",
            "threads",
            "mem_mb",
            "processes",
            "restriction_sites",
            "enzyme_provenance",
            "correct_nrounds",
            "gap_size",
            "mapq",
            "nm_exclusive",
            "hic_reference",
            "min_retained_fraction",
            "plot_bin_kb",
            "review_agp",
            "review_manifest",
        },
        "finishing": {
            "enabled",
            "gap_closing",
            "telomere_completion",
            "source_delivery",
            "telomere_motif",
            "threads",
            "mem_mb",
            "min_anchor",
            "min_identity",
            "min_mapq",
            "min_support",
            "local_window",
            "max_span",
        },
        "polishing": {"enabled", "backend", "min_contig_length", "on_missing_wgs", "threads", "mem_mb"},
        "assessment": {
            "enabled",
            "busco",
            "kmer",
            "read_mapping",
            "busco_lineage",
            "threads",
            "mem_mb",
            "min_qv",
            "min_busco_complete",
            "terminal_window",
            "kmer_size",
            "min_mapq",
            "coverage_window",
            "quast_min_contig",
            "telomere_min_length",
            "keep_intermediates",
            "terminal_max_distance",
            "terminal_min_repeats",
        },
    }
    check_fields(settings, {*groups, "ploidy", "genome_size_bp", "chromosomes_per_haplotype"}, path=path)
    for group, fields in groups.items():
        if group in settings:
            check_fields(settings[group], fields, path=f"{path}.{group}")

    def check(name, expected_type, group=None, **limits):
        """仅检查本次覆盖中实际出现的字段；保留缺省值，错误信息附带完整参数路径。"""
        values = settings if group is None else settings.get(group, {})
        if name in values:
            location = f"{path}.{group}.{name}" if group else f"{path}.{name}"
            check_value(values[name], expected_type, location, **limits)

    for name in ("ploidy", "genome_size_bp", "chromosomes_per_haplotype"):
        check(name, (int, type(None)), minimum=1)

    # 各组的类型、范围直接列在检查函数中，查规则时无需跳转到额外文件。
    for group in groups:
        check("threads", int, group, minimum=1 if group in {"assembly", "read_qc"} else 4)
        check("mem_mb", int, group, minimum=512 if group in {"assembly", "read_qc"} else 1024)
    check("backend", str, "assembly", choices=("hifiasm",))
    check("strategy", (str, type(None)), "assembly", choices=(None, "single", "diploid", "polyploid"))
    check("mode", str, "assembly", choices=("hifi_only", "hifi_hic", "hifi_ul_hic"))
    check("representation", str, "assembly", choices=("primary_alternate", "phased", "unitigs"))
    check("bloom_filter_bits", int, "assembly", minimum=0, maximum=39)
    check("purge_level", (int, type(None)), "assembly", minimum=0, maximum=3)

    for name in ("ont_min_length", "fastp_min_length", "hifi_adapter_length"):
        check(name, int, "read_qc", minimum=1)
    check("ont_min_quality", (int, float), "read_qc", minimum=0)
    check("fastp_quality", int, "read_qc", minimum=0, maximum=40)
    check("hifi_adapter_identity", (int, float), "read_qc", minimum=0, maximum=100)

    for group in ("scaffolding", "finishing", "polishing", "assessment"):
        check("enabled", bool, group)
    check("backend", str, "scaffolding", choices=("haphic", "gphase"))
    for name in ("processes", "gap_size", "nm_exclusive", "plot_bin_kb"):
        check(name, int, "scaffolding", minimum=1)
    check("restriction_sites", str, "scaffolding", pattern=r"[ACGTN]+(,[ACGTN]+)*")
    check("enzyme_provenance", str, "scaffolding", min_length=1)
    check("correct_nrounds", int, "scaffolding", minimum=0, maximum=10)
    check("mapq", int, "scaffolding", minimum=0, maximum=60)
    check("hic_reference", str, "scaffolding", choices=("combined", "per_haplotype"))
    check("min_retained_fraction", (int, float), "scaffolding", minimum=0, maximum=1)
    for name in ("review_agp", "review_manifest"):
        check(name, (str, type(None)), "scaffolding")

    for name in ("gap_closing", "telomere_completion"):
        check(name, bool, "finishing")
    check("source_delivery", (str, type(None)), "finishing")
    check("telomere_motif", (str, type(None)), "finishing", pattern=r"[ACGT]{2,20}")
    check("min_anchor", int, "finishing", minimum=500)
    check("min_identity", (int, float), "finishing", minimum=0.95, maximum=1)
    check("min_mapq", int, "finishing", minimum=20, maximum=60)
    check("min_support", int, "finishing", minimum=2)
    for name in ("local_window", "max_span"):
        check(name, int, "finishing", minimum=1000)

    check("backend", str, "polishing", choices=("nextpolish2",))
    check("min_contig_length", int, "polishing", minimum=1000)
    check("on_missing_wgs", str, "polishing", choices=("skip", "error"))

    for name in ("busco", "kmer", "read_mapping", "keep_intermediates"):
        check(name, bool, "assessment")
    check("busco_lineage", (str, type(None)), "assessment")
    check("min_qv", (int, float), "assessment", minimum=0)
    check("min_busco_complete", (int, float), "assessment", minimum=0, maximum=100)
    check("terminal_window", int, "assessment", minimum=100)
    check("kmer_size", int, "assessment", minimum=15, maximum=31)
    check("min_mapq", int, "assessment", minimum=0, maximum=60)
    for name in ("coverage_window", "quast_min_contig", "telomere_min_length"):
        check(name, int, "assessment", minimum=1)
    check("terminal_max_distance", int, "assessment", minimum=0)
    check("terminal_min_repeats", int, "assessment", minimum=3)


def validate_config(config: dict) -> None:
    """检查项目配置结构，再分别校验默认参数、命名参数组合和单样本覆盖；不填充默认值。"""
    required = {"samples_table", "readsets_table"}
    check_fields(
        config,
        required | {"prepared_input_cache", "defaults", "analysis_profiles", "samples"},
        required=required,
    )
    for name in ("samples_table", "readsets_table"):
        check_value(config[name], str, name, min_length=1)
    if "prepared_input_cache" in config:
        check_value(config["prepared_input_cache"], (str, type(None)), "prepared_input_cache")
    if "defaults" in config:
        validate_settings(config["defaults"], "defaults")
    for scope in ("analysis_profiles", "samples"):
        if scope in config:
            check_fields(config[scope], path=scope)
            for name, settings in config[scope].items():
                validate_settings(settings, f"{scope}.{name}")


def copy_labels(strategy, ploidy):
    """根据组装策略与真实倍性返回reference或hap1至hapP标签；单套参考不冒充亲本单倍型。"""
    if strategy == "single":
        return ["reference"]
    count = 2 if strategy == "diploid" else ploidy
    if not isinstance(count, int) or count < 2:
        raise ValueError("Phased output requires an explicit ploidy >= 2")
    return [f"hap{index}" for index in range(1, count + 1)]


def normalize_legacy_settings(raw):
    """移除旧配置中明确关闭的genome_profile占位设置，返回配置副本。

    开启或格式异常的旧设置报错，防止用户请求的分析被静默忽略；不新增基因组特征评估步骤。
    """
    if not isinstance(raw, dict):
        return raw
    result = copy.deepcopy(raw)
    groups = [result.get("defaults")]
    for key in ("analysis_profiles", "samples"):
        if isinstance(result.get(key), dict):
            groups.extend(result[key].values())
    for settings in groups:
        if not isinstance(settings, dict) or "genome_profile" not in settings:
            continue
        value = settings["genome_profile"]
        if isinstance(value, dict) and (not value or set(value) == {"enabled"} and value["enabled"] is False):
            del settings["genome_profile"]
        else:
            raise ConfigError(
                "genome_profile is outside the assembly workflow; "
                "provide known sample ploidy and chromosome information instead. "
                "Only the disabled placeholder from old configurations is accepted."
            )
    return result


def resolve_config(config_path, repository, run_dir, selected, cores, mem_mb):
    """读取项目配置、样本表和测序数据表，返回完整组装运行计划。

    参数为配置路径、代码根目录、运行目录、所选样本及资源预算。按默认、参数组合、样本覆盖
    逐级合并设置，校验倍性、路线和数据要求；此时只查元数据，完整读段解析留给步骤01。
    """
    config_path, repository, run_dir = map(lambda p: Path(p).resolve(), (config_path, repository, run_dir))
    from omics_workflows.assembly.lib.reads import describe_prepared_input

    boundary = source_boundary(repository)
    validate_run_directory(
        boundary,
        run_dir,
        cores,
        mem_mb,
        LAYOUT_VERSION,
        historical_marker="layout_migration/status.json",
    )
    raw = load_configuration(config_path, validate_config, normalize_legacy_settings)
    defaults = load_defaults(repository / "assembly/defaults.yaml")
    sample_path = (config_path.parent / raw["samples_table"]).resolve()
    readset_path = (config_path.parent / raw["readsets_table"]).resolve()
    rows = table(sample_path, {"sample_id", "species", "profile"}, {"sample_id", "species", "profile"})
    samples, selected = select_samples(rows, selected, raw.get("samples", {}))
    reads = table(
        readset_path,
        {"sample_id", "readset_id", "library_id", "technology", "format", "file1"},
        {
            "sample_id",
            "readset_id",
            "library_id",
            "technology",
            "format",
            "file1",
            "file2",
            "md5_1",
            "md5_2",
            "prepared_manifest",
        },
    )
    by_sample = group_readsets(reads, samples, {"hifi", "ont", "hic", "wgs"})
    result = {
        "layout_version": LAYOUT_VERSION,
        "config_path": str(config_path),
        "repository": str(repository),
        "framework": {"name": "omics-workflows", "version": __version__, "workflow": "assembly"},
        "run_dir": str(run_dir),
        "budget": {"cores": cores, "mem_mb": mem_mb},
        "sample_order": selected,
        "samples": {},
    }
    used_paths = set()
    for sid in selected:
        row = samples[sid]
        profile = row["profile"]
        profiles = raw.get("analysis_profiles", {})
        if profile and profile not in profiles:
            raise ConfigError(f"{sid}: unknown analysis profile {profile}")
        # 参数优先级：模块默认值 < 项目defaults < 分析参数组合 < 单样本覆盖。
        settings = merge(defaults, raw.get("defaults", {}))
        settings = merge(settings, profiles.get(profile, {}))
        settings = merge(settings, raw.get("samples", {}).get(sid, {}))
        asm = settings["assembly"]
        joint = asm["mode"] in {"hifi_hic", "hifi_ul_hic"}
        explicit_strategy = asm["strategy"] is not None
        strategy = asm["strategy"] or (
            "polyploid"
            if settings["scaffolding"]["backend"] == "gphase"
            else "diploid" if joint and asm["representation"] == "phased" else "single"
        )
        asm["strategy"] = strategy
        polyploid = strategy == "polyploid"
        # 显式策略同时约束序列表示和挂载工具，避免产生相互矛盾的配置。
        if explicit_strategy:
            asm["representation"] = {
                "single": "primary_alternate",
                "diploid": "phased",
                "polyploid": "unitigs",
            }[strategy]
            settings["scaffolding"]["backend"] = "gphase" if polyploid else "haphic"
        if not joint and not (
            asm["mode"] == "hifi_only"
            and strategy == "single"
            and asm["representation"] == "primary_alternate"
        ):
            raise ConfigError(f"{sid}: supported data modes are hifi_only, hifi_hic and hifi_ul_hic")
        if joint and (not settings["ploidy"] or not settings["chromosomes_per_haplotype"]):
            raise ConfigError(
                f"{sid}: chromosome assembly requires explicit ploidy=2 for diploid, "
                "known ploidy for other strategies, and chromosomes_per_haplotype; no automatic fallback"
            )
        if strategy == "diploid" and settings["ploidy"] != 2:
            raise ConfigError(f"{sid}: diploid strategy requires ploidy=2")
        if polyploid and (
            not isinstance(settings["ploidy"], int)
            or settings["ploidy"] < 3
            or asm["mode"] != "hifi_ul_hic"
            or asm["representation"] != "unitigs"
        ):
            raise ConfigError(
                f"{sid}: polyploid assembly strategy: GPhase requires ploidy>=3, hifi_ul_hic and unitigs"
            )
        if settings["ploidy"] and settings["ploidy"] > 2 and not polyploid:
            raise ConfigError(f"{sid}: polyploid assembly strategy must be selected explicitly for ploidy>2")
        settings["output_copies"] = len(copy_labels(strategy, settings["ploidy"]))
        scaffold = settings["scaffolding"]
        if joint != scaffold["enabled"]:
            raise ConfigError(f"{sid}: chromosome data modes require scaffolding.enabled=true")
        if scaffold["hic_reference"] == "per_haplotype" and strategy != "diploid":
            raise ConfigError(f"{sid}: per_haplotype Hi-C reference requires diploid phased HapHiC assembly")
        if polyploid and asm["purge_level"] not in (None, 0):
            raise ConfigError(f"{sid}: raw unitig route requires purge_level=0 or null")
        if polyploid and (scaffold["gap_size"] != 100 or scaffold["correct_nrounds"]):
            raise ConfigError(f"{sid}: GPhase requires gap_size=100 and correct_nrounds=0")
        if polyploid and (scaffold["review_agp"] or scaffold["review_manifest"]):
            raise ConfigError(f"{sid}: reviewed polyploid AGP import is not implemented")
        if joint and (cores < 4 or mem_mb < 32000):
            raise ConfigError(f"{sid}: chromosome assembly requires >=4 cores and >=32000 MB")
        if any(c.isspace() for c in str(run_dir)):
            raise ConfigError("Native assembly tools require a whitespace-free run directory")
        if bool(scaffold["review_agp"]) != bool(scaffold["review_manifest"]):
            raise ConfigError(f"{sid}: review_agp and review_manifest must be supplied together")
        finish, polish, assess = (settings[k] for k in ("finishing", "polishing", "assessment"))
        if not joint and any(assess[key] for key in ("busco", "kmer", "read_mapping")):
            raise ConfigError(
                f"{sid}: advanced assessment is unavailable for the basic contig statistics only endpoint"
            )
        if assess["busco"] and not assess["busco_lineage"]:
            raise ConfigError(f"{sid}: BUSCO needs an existing offline lineage directory")
        if not joint and (finish["enabled"] or polish["enabled"]):
            raise ConfigError(
                f"{sid}: finishing and polishing are not implemented for the contig-only endpoint"
            )
        if not finish["enabled"] and (
            finish["gap_closing"] or finish["telomere_completion"] or finish["source_delivery"]
        ):
            raise ConfigError(f"{sid}: enable finishing before requesting finishing operations")
        if finish["enabled"] and not finish["telomere_motif"]:
            raise ConfigError(f"{sid}: finishing requires an explicit telomere_motif")
        if finish["local_window"] < finish["min_anchor"]:
            raise ConfigError(f"{sid}: local_window must be >= min_anchor")
        if assess["terminal_max_distance"] >= assess["terminal_window"]:
            raise ConfigError(f"{sid}: terminal_max_distance must be smaller than terminal_window")
        active_groups = {"assembly"}
        if joint:
            active_groups.update(["scaffolding", "assessment"])
        if finish["enabled"]:
            active_groups.add("finishing")
        if polish["enabled"]:
            active_groups.add("polishing")
        for group_name in ("assembly", "scaffolding", "read_qc", "finishing", "polishing", "assessment"):
            group = settings[group_name]
            if group_name in active_groups and group["mem_mb"] > mem_mb:
                raise ConfigError(f"{sid}: module memory exceeds task budget")
            if group_name not in active_groups:
                group["mem_mb"] = min(group["mem_mb"], mem_mb)
            group["threads"] = min(group["threads"], cores)
        scaffold["processes"] = min(scaffold["processes"], scaffold["threads"])
        for group, field in [
            (scaffold, "review_agp"),
            (scaffold, "review_manifest"),
            (finish, "source_delivery"),
            (assess, "busco_lineage"),
        ]:
            if group[field]:
                resolved = (config_path.parent / group[field]).resolve()
                if not resolved.exists() or resolved.is_relative_to(run_dir):
                    raise ConfigError(f"{sid}: {field} must be an existing external path")
                group[field] = str(resolved)
        if finish["source_delivery"] and not Path(finish["source_delivery"]).is_file():
            raise ConfigError(f"{sid}: source_delivery must point to manifest.json")
        if assess["busco_lineage"] and not (Path(assess["busco_lineage"]) / "dataset.cfg").is_file():
            raise ConfigError(f"{sid}: BUSCO lineage requires dataset.cfg")
        if asm["mem_mb"] > mem_mb:
            raise ConfigError(f'{sid}: assembly requires {asm["mem_mb"]} MB; ' f"task budget is {mem_mb} MB")
        # Bloom storage alone is 2**f bits; retain at least 512 MB for other state.
        bloom_mb = (2 ** asm["bloom_filter_bits"]) / 8 / 1_000_000 if asm["bloom_filter_bits"] else 0
        if bloom_mb + 512 > asm["mem_mb"]:
            raise ConfigError(f"{sid}: bloom_filter_bits exceeds the assembly memory budget")
        asm["threads"] = min(asm["threads"], cores)
        active = {}
        technologies = (
            {"hifi", "hic"} | ({"ont"} if asm["mode"] == "hifi_ul_hic" else set()) if joint else {"hifi"}
        )
        required_technologies = set(technologies)
        if joint:
            technologies.add("wgs")
        for read in by_sample.get(sid, []):
            if read["technology"] not in technologies:
                continue
            rid = read["readset_id"]
            if read["technology"] in {"hic", "wgs"}:
                if read["format"] != "fastq" or not read.get("file2"):
                    raise ConfigError(f"{sid}/{rid}: Hi-C/WGS requires paired FASTQ file1 and file2")
            elif (
                read["format"] not in {"fastq", "bam"}
                or read.get("file2")
                or (read["technology"] == "ont" and read["format"] != "fastq")
            ):
                raise ConfigError(f"{sid}/{rid}: HiFi input must be single-read FASTQ or unpaired BAM")
            path = (readset_path.parent / read["file1"]).resolve()
            if not path.is_file() or path.stat().st_size == 0:
                raise ConfigError(f"{sid}/{rid}: missing or empty input {path}")
            if path.is_relative_to(run_dir):
                raise ConfigError(f"{sid}/{rid}: source input must be outside the run directory")
            if path in used_paths:
                raise ConfigError(f"Duplicate active source path (would duplicate reads): {path}")
            used_paths.add(path)
            md5 = read.get("md5_1", "").lower()
            if md5 and not re.fullmatch("[0-9a-f]{32}", md5):
                raise ConfigError(f"{sid}/{rid}: invalid md5_1")
            active[rid] = {
                **read,
                "file1": str(path),
                "expected_md5": md5 or None,
                "source_signature": source_signature(path),
            }
            if read.get("prepared_manifest"):
                if not raw.get("prepared_input_cache"):
                    raise ConfigError(
                        f"{sid}/{rid}: prepared_manifest requires an external prepared_input_cache"
                    )
                try:
                    manifest = (readset_path.parent / read["prepared_manifest"]).resolve()
                    cache = (config_path.parent / raw["prepared_input_cache"]).resolve()
                    active[rid]["prepared_manifest"] = str(manifest)
                    active[rid]["prepared_reuse"] = describe_prepared_input(
                        active[rid], manifest, cache, run_dir, boundary
                    )
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    raise ConfigError(f"{sid}/{rid}: {exc}") from exc
            if read["technology"] in {"hic", "wgs"}:
                mate = (readset_path.parent / read["file2"]).resolve()
                md5_2 = read.get("md5_2", "").lower()
                if not mate.is_file() or not mate.stat().st_size or mate.is_relative_to(run_dir):
                    raise ConfigError(f"{sid}/{rid}: missing or empty external mate {mate}")
                if mate in used_paths:
                    raise ConfigError(f"Duplicate active source path (would duplicate reads): {mate}")
                if md5_2 and not re.fullmatch("[0-9a-f]{32}", md5_2):
                    raise ConfigError(f"{sid}/{rid}: invalid md5_2")
                used_paths.add(mate)
                active[rid].update(
                    file2=str(mate), expected_md5_2=md5_2 or None, source_signature_2=source_signature(mate)
                )
        if not active:
            raise ConfigError(f"{sid}: hifi_only requires at least one HiFi readset")
        missing = required_technologies - {r["technology"] for r in active.values()}
        if missing:
            raise ConfigError(f"{sid}: required read technologies missing: {sorted(missing)}")
        if (
            polish["enabled"]
            and polish["on_missing_wgs"] == "error"
            and not any(r["technology"] == "wgs" for r in active.values())
        ):
            raise ConfigError(f"{sid}: NextPolish2 requires ordinary WGS short reads, not Hi-C")
        physical_inputs = {}
        for rid, spec in active.items():
            paths = {spec["file1"]}
            if spec.get("file2"):
                paths.add(spec["file2"])
            if spec.get("prepared_reuse"):
                paths.add(spec["prepared_reuse"]["identity"]["prepared_fastq"])
            for path in paths:
                if path in physical_inputs and physical_inputs[path] != rid:
                    raise ConfigError(f"{sid}: duplicate raw/prepared source path across readsets")
                physical_inputs[path] = rid
        hic_libraries = [r["library_id"] for r in active.values() if r["technology"] == "hic"]
        if len(hic_libraries) != len(set(hic_libraries)):
            raise ConfigError(
                f"{sid}: combine lanes of each Hi-C library before registering; "
                "duplicate marking must cover the complete physical library"
            )
        result["samples"][sid] = {
            "sample_id": sid,
            "species": row["species"],
            "profile": profile,
            "settings": settings,
            "readsets": active,
            "unused_readsets": [
                r["readset_id"] for r in by_sample.get(sid, []) if r["technology"] not in technologies
            ],
        }
    strategies = {sample["settings"]["assembly"]["strategy"] for sample in result["samples"].values()}
    result["workflow_route"] = next(iter(strategies)) if len(strategies) == 1 else "mixed"
    return result


def required_tools(plan, target):
    """按所选样本路线、数据技术和分析终点返回所需工具名称集合，不给未启用的分析强加依赖。"""
    tools = {"samtools", "hifiasm", "snakemake"}
    full = target in ("report", "t2t")
    for sample in plan["samples"].values():
        settings = sample["settings"]
        technologies = {read["technology"] for read in sample["readsets"].values()}
        if "hifi" in technologies:
            tools.update(["NanoPlot", "hifiadapterfilt.sh", "blastn", "bamtools", "pigz"])
        if "ont" in technologies:
            tools.update(["NanoPlot", "NanoFilt", "pigz", "gzip"])
        if technologies & {"hic", "wgs"}:
            tools.add("fastp")
        if any(read.get("prepared_reuse") for read in sample["readsets"].values()):
            tools.add("seqkit")
        if settings["scaffolding"]["enabled"] and target != "contigs":
            tools.update(["bwa", "samblaster", "haphic", "filter_bam"])
            if settings["assembly"]["strategy"] == "polyploid":
                tools.update(["gphase", "mosdepth", "minimap2", "yahs"])
        if not full or not settings["scaffolding"]["enabled"]:
            continue
        tools.update(
            [
                "tidk",
                "quast.py",
                "busco",
                "minimap2",
                "mosdepth",
                "gt",
                "LTR_FINDER_parallel",
                "LTR_retriever",
            ]
        )
        if "wgs" in technologies:
            tools.update(["meryl", "merqury.sh"])
        if settings["finishing"]["enabled"] and settings["finishing"]["gap_closing"]:
            tools.add("quartet")
        if settings["polishing"]["enabled"] and "wgs" in technologies:
            tools.update(["yak", "nextPolish2"])
    return tools


def source_locks(prefix):
    """根据环境前缀返回GPhase和T2T扩展的源码锁规范，供通用校验函数使用。"""
    return [
        {
            "root": prefix / "opt/GPhase",
            "files_key": "source_files",
            "files_base": "source",
            "label": "GPhase",
        },
        {"root": prefix / "opt/T2T", "files_key": "files", "label": "T2T"},
    ]


def verify_source_locks(prefix, tools):
    """使用组装模块的扩展规范校验源码锁，并更新本次工具记录。"""
    workflow.verify_source_locks(tools, source_locks(prefix))


def environment_tools(plan, target, strict):
    """检查当前环境中的组装依赖和扩展锁，返回工具记录字典。

    根据分析终点检查数据库与端粒设置；strict决定缺失依赖是报错还是仅供预演报告。
    """
    prefix = Path(sys.prefix).resolve()
    prepend_search_paths([prefix / "opt/GPhase/bin", prefix / "opt/T2T/bin"], [prefix / "opt/GPhase/python"])

    def fallback(name):
        """定位HapHiC源码目录中的filter_bam；其他工具或文件不存在时返回None。"""
        haphic = shutil.which("haphic") if name == "filter_bam" else None
        candidate = Path(haphic).resolve().parent / "utils/filter_bam" if haphic else None
        return str(candidate) if candidate and candidate.is_file() else None

    tools, missing = inspect_tools(
        required_tools(plan, target),
        prefix,
        versions={name: ["--version"] for name in ("hifiasm", "samtools", "tidk", "fastp", "mosdepth")},
        fallback=fallback,
    )
    haphic = tools.get("haphic", {})
    if haphic.get("path"):
        lock = Path(haphic["path"]).resolve().parent / "SOURCE_LOCK.json"
        if lock.is_file():
            haphic["source_lock"] = json.loads(lock.read_text(encoding="utf-8"))
    verify_source_locks(prefix, tools)
    if target in ("report", "t2t"):
        for sample in plan["samples"].values():
            settings = sample["settings"]
            if settings["scaffolding"]["enabled"]:
                if not settings["assessment"]["enabled"]:
                    missing.append(
                        f'{sample["sample_id"]}: final seven-category assessment cannot be disabled'
                    )
                if not settings["assessment"]["busco_lineage"]:
                    missing.append(f'{sample["sample_id"]}: assessment.busco_lineage')
                if not settings["finishing"]["telomere_motif"]:
                    missing.append(f'{sample["sample_id"]}: finishing.telomere_motif')
    require_available(missing, strict)
    return tools


def verify_cached_outputs(plan):
    """为通用缓存检查器提供组装各步骤的清单路径和专用核验函数，不凭目录存在认定完成。"""
    from omics_workflows.assembly.lib.postassembly import assembly_bundle, read_bundle

    def chromosomes(path, _sample):
        """校验染色体交付清单、序列集合及校验值，异常直接阻止复用该步骤。"""
        assembly_bundle(path)

    rules = {name: ("manifest.json", None) for name in ("03_assembly/contigs", "07_report")}
    rules.update({name: ("manifest.json", read_bundle) for name in ("01_inputs", "02_read_qc")})
    rules.update(
        {
            name: ("manifest.json", chromosomes)
            for name in ("03_assembly/delivery", "04_finishing", "05_polishing")
        }
    )
    rules["06_assessment"] = ("assessment.json", None)
    workflow.verify_cached_outputs(plan, rules)


def sample_metadata(sample, stage, tools, code):
    """为一个样本和一个步骤构建参数、工具与代码来源字典，只纳入该步骤实际相关的设置。"""
    settings = sample["settings"]
    identity = {key: settings[key] for key in ("ploidy", "output_copies", "chromosomes_per_haplotype")}
    selected = {"assembly": settings["assembly"], "scaffolding": settings["scaffolding"]}
    reads = {}
    if stage == "inputs":
        selected = {"assembly": {"strategy": settings["assembly"]["strategy"]}}
        reads = {
            rid: {**record, "samtools": tools["samtools"], "seqkit": tools.get("seqkit")}
            for rid, record in sample["readsets"].items()
        }
    elif stage == "read_qc":
        selected = {"read_qc": settings["read_qc"]}
    elif stage == "contigs":
        selected["scaffolding"] = {"enabled": settings["scaffolding"]["enabled"]}
    elif stage == "scaffolding":
        selected["assembly"] = {"strategy": settings["assembly"]["strategy"]}
    elif stage == "finishing":
        selected = {
            "finishing": settings["finishing"],
            "assessment": {
                key: settings["assessment"][key]
                for key in ("terminal_window", "terminal_max_distance", "terminal_min_repeats")
            },
        }
    elif stage == "polishing":
        selected = {"polishing": settings["polishing"]}
    elif stage == "assessment":
        selected = {
            "assessment": settings["assessment"],
            "scaffolding": {
                key: value
                for key, value in settings["scaffolding"].items()
                if key not in ("review_agp", "review_manifest")
            },
            "finishing": {"telomere_motif": settings["finishing"]["telomere_motif"]},
        }
    return {
        "sample_id": sample["sample_id"],
        "settings": {**identity, **selected},
        "readsets": reads,
        "tools": tools,
        "code": code,
    }


def write_plan_metadata(plan, repository, run_dir):
    """把组装步骤定义和样本记录构造函数交给公共模块，写入各步骤参数JSON。"""
    workflow.write_plan_metadata(plan, repository, run_dir, STAGES, sample_metadata)


def main(argv=None, prog=None):
    """解析组装命令，按需切换分析环境，在包资源有效期内执行预检或分析；返回退出码。"""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argument_parser(
        description="基因组组装：检查输入与环境，按配置执行组装、挂载、完善、纠错、评估和报告。",
        prog=prog,
        epilog=(
            "分析终点（--target）：\n"
            "  contigs      完成数据准备、质控和初步组装，输出contig或原始unitig图。\n"
            "  chromosomes  完成步骤03染色体挂载交付，不执行步骤04至07。\n"
            "  report       默认值；执行所选分析路线并生成样本及批次报告。\n"
            "               HiFi-only路线仅提供contig统计；染色体路线执行七类最终评估。\n"
            "  t2t          染色体路线完整后处理及报告的兼容目标，不表示T2T认证。\n\n"
            "路径与分析参数：\n"
            "  --config、--run-dir、--env-prefix的相对路径按当前工作目录解析。\n"
            "  YAML内的表格路径相对YAML；读段路径相对readsets.tsv。\n"
            "  组装策略、已知倍性、数据库和各工具参数通过YAML配置。\n\n"
            "项目目录内的预检示例：\n"
            "  omics-workflows assembly --config config/assembly.yaml \\\n"
            "    --run-dir runs/assembly001 --cores 8 --mem-mb 64000 --check"
        ),
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="类型：YAML文件路径；必填。分析配置；其中样本表、测序数据表的相对路径以此文件所在目录为基准。",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="类型：目录路径；必填。分析结果与运行记录目录；须位于代码目录之外且路径不含空白字符。恢复运行时使用原目录。",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        help="类型：样本ID字符串列表，以空格分隔；默认选择samples.tsv中的全部样本。只分析指定样本。",
    )
    parser.add_argument(
        "--cores",
        required=True,
        type=int,
        help="类型：整数，至少1；必填。流程可用CPU总数，染色体组装至少4。应与实际申请资源一致；本参数不提交集群作业。",
    )
    parser.add_argument(
        "--mem-mb",
        required=True,
        type=int,
        help="类型：整数，单位MB，至少1024；必填。流程可用总内存，染色体组装至少32000；本参数不向集群申请内存。",
    )
    parser.add_argument(
        "--target",
        choices=["contigs", "chromosomes", "t2t", "report"],
        default="report",
        help="类型：字符串枚举；默认report。选择分析终点，各取值含义见下方说明。",
    )
    parser.add_argument(
        "--env-prefix",
        type=Path,
        help="类型：环境目录路径；默认使用当前Python环境。切换到已有分析环境，不安装软件；框架仍使用本次调用的安装版本。",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。生成运行计划并预演步骤依赖，不执行分析；缺失依赖会列出。",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="类型：布尔开关；默认关闭，无需赋值。仅检查配置、文件存在性和软件依赖，不创建运行目录；不做全量读段解析或校验。与--dry-run同时使用时本选项优先。",
    )
    args = parser.parse_args(argv)
    if args.env_prefix:
        prefix = args.env_prefix.resolve()
        executable = prefix / "bin/python"
        if not executable.is_file():
            parser.error(f"Environment Python not found: {executable}")
        env = dict(
            os.environ,
            PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""),
            CONDA_PREFIX=str(prefix),
            PYTHONNOUSERSITE="1",
            PYTHONDONTWRITEBYTECODE="1",
            PYTHONPATH=os.pathsep.join(
                filter(None, [str(Path(__file__).resolve().parents[2]), os.environ.get("PYTHONPATH")])
            ),
        )
        if Path(sys.prefix).resolve() != prefix or os.environ.get("CONDA_PREFIX") != str(prefix):
            os.execve(str(executable), [str(executable), "-m", "omics_workflows.cli", "assembly", *argv], env)
        os.environ.update(env)

    with package_directory() as repository:
        return execute(args, parser, repository)


def execute(args, parser, repository):
    """以包资源根目录执行配置解析、工具检查和Snakemake启动，返回退出码。

    记录框架版本及代码指纹；--check仅预检，--dry-run构建分析依赖，实际执行前核验已有结果。
    """

    try:
        plan = resolve_config(args.config, repository, args.run_dir, args.samples, args.cores, args.mem_mb)

        if args.target in ("chromosomes", "t2t") and not all(
            sample["settings"]["scaffolding"]["enabled"] for sample in plan["samples"].values()
        ):
            raise ConfigError(f"--target {args.target} requires a chromosome route for every sample")
        tools = environment_tools(plan, args.target, strict=not args.dry_run or args.check)
        snakemake = Path(sys.prefix) / "bin/snakemake"
        if not snakemake.is_file():
            raise ConfigError("snakemake missing from the unified environment")
    except (ConfigError, OSError, ImportError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"PRECHECK FAILED: {exc}\n")
    code_paths = [
        repository / "__init__.py",
        repository / "cli.py",
        *sorted((repository / "common").rglob("*.py")),
        *sorted((repository / "assembly").rglob("*.py")),
        *sorted((repository / "assembly").rglob("*.smk")),
        repository / "assembly/Snakefile",
        repository / "assembly/defaults.yaml",
    ]
    plan["code"] = {
        str(p.relative_to(repository)): hashlib.sha256(p.read_bytes()).hexdigest() for p in code_paths
    }

    plan["tools"] = tools
    plan["target"] = args.target
    plan["plan_id"] = digest(plan)
    print(
        f'PRECHECK PASS: samples={",".join(plan["sample_order"])} '
        f"cores={args.cores} mem_mb={args.mem_mb} "
        f'modes={",".join(plan["samples"][s]["settings"]["assembly"]["mode"] for s in plan["sample_order"])} '
        "executor=local",
        flush=True,
    )
    if args.check:
        return 0

    def prepare():
        """取得运行锁后核验历史结果并写入步骤参数，避免多个启动器同时改写同一运行目录。"""
        verify_cached_outputs(plan)
        write_plan_metadata(plan, repository, Path(plan["run_dir"]))

    try:
        return workflow.execute_workflow(
            plan,
            repository / "assembly/Snakefile",
            "assembly_" + args.target,
            prepare,
            args.dry_run,
        )
    except (ConfigError, OSError, ValueError) as error:
        parser.exit(2, f"RUN PREPARATION FAILED: {error}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
