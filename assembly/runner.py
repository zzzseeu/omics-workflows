"""Assembly-specific validation, stage definitions and workflow launcher."""

import argparse
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
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# pylint: disable=wrong-import-position
from common import workflow
from common.workflow import (
    ConfigError,
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


def copy_labels(strategy, ploidy):
    """Keep an unphased reference distinct from biological haplotypes."""
    if strategy == "single":
        return ["reference"]
    count = 2 if strategy == "diploid" else ploidy
    if not isinstance(count, int) or count < 2:
        raise ValueError("Phased output requires an explicit ploidy >= 2")
    return [f"hap{index}" for index in range(1, count + 1)]


def normalize_legacy_settings(raw):
    """Accept an old disabled profiling placeholder without retaining a stage.

    Earlier templates carried this unimplemented switch. Removing a disabled
    placeholder preserves those configurations; enabled or malformed settings
    must fail instead of silently dropping a requested analysis.
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
    """Metadata-only preflight; full checksums/FASTQ parsing belong in compute rules."""
    config_path, repository, run_dir = map(lambda p: Path(p).resolve(), (config_path, repository, run_dir))
    from assembly.lib.reads import describe_prepared_input

    validate_run_directory(
        repository,
        run_dir,
        cores,
        mem_mb,
        LAYOUT_VERSION,
        historical_marker="layout_migration/status.json",
        migration_tool="docs/legacy/assembly_migrate_layout.py",
    )
    raw = load_configuration(config_path, repository / "assembly/schema.yaml", normalize_legacy_settings)
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
        "schema_version": 1,
        "layout_version": LAYOUT_VERSION,
        "config_path": str(config_path),
        "repository": str(repository),
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
                        active[rid], manifest, cache, run_dir, repository
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
    workflow.verify_source_locks(tools, source_locks(prefix))


def environment_tools(plan, target, strict):
    prefix = Path(sys.prefix).resolve()
    prepend_search_paths([prefix / "opt/GPhase/bin", prefix / "opt/T2T/bin"], [prefix / "opt/GPhase/python"])

    def fallback(name):
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
    from assembly.lib.postassembly import assembly_bundle, read_bundle

    def chromosomes(path, _sample):
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
    workflow.write_plan_metadata(plan, repository, run_dir, STAGES, sample_metadata)


def main(argv=None, prog=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__, prog=prog)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--samples", nargs="+", help="Default: all samples in samples.tsv")
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--mem-mb", required=True, type=int)
    parser.add_argument("--target", choices=["contigs", "chromosomes", "t2t", "report"], default="report")
    parser.add_argument("--env-prefix", type=Path, help="Optional unified environment prefix")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--check", action="store_true", help="Metadata/dependency preflight only")
    args = parser.parse_args(argv)
    repository = Path(__file__).resolve().parents[1]
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
        )
        if Path(sys.prefix).resolve() != prefix or os.environ.get("CONDA_PREFIX") != str(prefix):
            os.execve(str(executable), [str(executable), str(repository / "run"), "assembly", *argv], env)
        os.environ.update(env)
    sys.path.insert(0, str(repository))

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
        repository / "run",
        *sorted((repository / "common").rglob("*.py")),
        *sorted((repository / "assembly").rglob("*.py")),
        *sorted((repository / "assembly").rglob("*.smk")),
        repository / "assembly/Snakefile",
        repository / "assembly/schema.yaml",
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
        verify_cached_outputs(plan)
        write_plan_metadata(plan, repository, Path(plan["run_dir"]))

    try:
        return workflow.execute_workflow(
            plan, repository / "assembly/Snakefile", "assembly_" + args.target, prepare, args.dry_run
        )
    except (ConfigError, OSError, ValueError) as error:
        parser.exit(2, f"RUN PREPARATION FAILED: {error}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
