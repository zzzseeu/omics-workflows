#!/usr/bin/env python3
"""Migrate a completed layout-1 run without reassembling or fabricating engine metadata.

Default: print a move plan only. --apply performs a journaled, resumable migration.
The resulting completed run is historical/read-only to the current analysis launcher;
new analyses use a new run directory. Original commands, logs and code hashes remain evidence.
"""
import argparse
import fcntl
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def stat_record(path):
    stat = path.stat()
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "inode": stat.st_ino, "device": stat.st_dev}


def mappings(run):
    pairs = [("reports", "10_reports"), ("logs/batch_report.log", "logs/10_reports/batch_report.log")]
    for sample in load(run / "config/plan.json")["sample_order"]:
        base = f"samples/{sample}"
        pairs += [
            (f"{base}/{a}", f"{base}/{b}")
            for a, b in [
                ("inputs", "01_inputs"),
                ("qc/reads", "02_read_qc"),
                ("assembly", "04_assembly"),
                ("qc/assembly", "09_assessment"),
                ("report", "10_report"),
                ("logs/prepare", "logs/01_inputs"),
                ("logs/read_qc", "logs/02_read_qc"),
                ("logs/assembly", "logs/04_assembly"),
                ("logs/index", "logs/04_assembly/index"),
                ("logs/assessment", "logs/09_assessment"),
                ("logs/report", "logs/10_report"),
                ("benchmarks/prepare", "benchmarks/01_inputs"),
                ("benchmarks/read_qc", "benchmarks/02_read_qc"),
                ("benchmarks/assembly", "benchmarks/04_assembly"),
                ("benchmarks/assessment", "benchmarks/09_assessment"),
            ]
        ]
    return pairs


def translated(text, run, pairs):
    for old, new in sorted(pairs, key=lambda pair: len(pair[0]), reverse=True):
        for prefix in [str(run) + "/", ""]:
            source, target = prefix + old, prefix + new
            if text == source or text.startswith(source + "/"):
                return target + text[len(source) :]
    return text


def rewrite(value, run, pairs):
    if isinstance(value, str):
        return translated(value, run, pairs)
    if isinstance(value, list):
        return [rewrite(x, run, pairs) for x in value]
    if isinstance(value, dict):
        # These are historical facts, not current artifact pointers.
        keep = {
            "command",
            "code",
            "assembly_id",
            "plan_id",
            "source_signature",
            "source_md5",
            "source_sha256",
            "started_utc",
            "finished_utc",
        }
        return {k: v if k in keep else rewrite(v, run, pairs) for k, v in value.items()}
    return value


def check_artifact_links(run):
    for sid in load(run / "config/plan.json")["sample_order"]:
        report = load(run / f"samples/{sid}/10_report/manifest.json")
        assert Path(report["assembly_manifest"]).is_file(), report["assembly_manifest"]
        for item in report["sets"]:
            for key in ["fasta", "fai"]:
                path = Path(item[key])
                assert path.is_file(), str(path)
                expected = item.get("sha256" if key == "fasta" else "fai_sha256")
                if expected:
                    assert sha256(path) == expected, f"Content differs from analysis manifest: {path}"
        assert (run / f"samples/{sid}/10_report/report.html").is_file()
    batch = run / "10_reports/batch.json"
    if batch.is_file():
        for sample in load(batch)["samples"]:
            assert Path(sample["manifest"]).is_file(), sample["manifest"]


def migrate(run, apply=False):
    run = Path(run).resolve()
    state = load(run / "run_status.json")
    if state.get("status") != "COMPLETED" or state.get("exit_code") != 0:
        raise ValueError("Only successfully completed runs may be migrated")
    status_path = run / "layout_migration/status.json"
    if not status_path.exists():
        if load(run / "config/plan.json").get("layout_version", 1) >= 2:
            raise ValueError("Run already uses numbered output paths")
        pairs = [(a, b) for a, b in mappings(run) if (run / a).exists()]
        if not pairs:
            raise ValueError("No layout-1 outputs found")
        for old, new in pairs:
            if (run / new).exists():
                raise ValueError(f"Migration target already exists: {new}")
        if not apply:
            return {"status": "PLAN_ONLY", "moves": pairs}
    if not apply:
        return load(status_path)
    with (run / ".launcher.lock").open("a", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if status_path.exists():
            journal = load(status_path)
            pairs = journal["moves"]
            if journal["status"] == "MIGRATED":
                check_artifact_links(run)
                return journal
        else:
            files, metadata = {}, {}
            for path in run.rglob("*"):
                if not path.is_file() or path.is_symlink() or ".snakemake" in path.parts:
                    continue
                rel = str(path.relative_to(run))
                if rel.startswith("layout_migration/") or rel == ".launcher.lock":
                    continue
                files[rel] = stat_record(path)
                if path.suffix in {".fa", ".fai", ".gfa", ".bin"}:
                    # Hash published FASTA/indices/GFAs; binary caches are protected by
                    # inode/size/mtime and are never opened for writing during migration.
                    if path.suffix != ".bin":
                        files[rel]["sha256"] = sha256(path)
                if path.suffix in {".json", ".tsv", ".md", ".html"}:
                    metadata[rel] = sha256(path)
            journal = {
                "status": "MIGRATING",
                "from_layout": 1,
                "to_layout": 2,
                "run_dir": str(run),
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "moves": pairs,
                "files_before": files,
                "metadata_before": metadata,
                "analysis_status": state,
                "preserve_analysis_identity": True,
            }
            save(status_path, journal)
        # Every move is atomic within this run; journal lets an interrupted migration resume.
        for old, new in pairs:
            source, target = run / old, run / new
            if source.exists() and target.exists():
                raise ValueError(f"Both source and target exist: {old}, {new}")
            if source.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                source.replace(target)
            elif not target.exists():
                raise ValueError(f"Neither source nor target exists: {old}, {new}")
        modified = []
        for rel in journal["metadata_before"]:
            path = run / translated(rel, run, pairs)
            archive = run / "layout_migration/original_metadata" / rel
            if not archive.exists():
                archive.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, archive)
            if path.suffix == ".json":
                old_value = load(archive)
                new_value = rewrite(old_value, run, pairs)
                if new_value != old_value:
                    save(path, new_value)
                    modified.append(rel)
            elif path.suffix == ".tsv":
                text = archive.read_text(encoding="utf-8")
                lines = [
                    "\t".join(translated(cell, run, pairs) for cell in line.split("\t"))
                    for line in text.splitlines()
                ]
                replacement = "\n".join(lines) + ("\n" if text.endswith("\n") else "")
                if replacement.splitlines() != text.splitlines():
                    path.write_text(replacement, encoding="utf-8")
                    modified.append(rel)
        # Remove only empty containers left by the split of qc/reads and qc/assembly.
        for sid in load(run / "config/plan.json")["sample_order"]:
            leftover = run / f"samples/{sid}/qc"
            if leftover.is_dir() and not any(leftover.iterdir()):
                leftover.rmdir()
        preserved = 0
        for rel, before in journal["files_before"].items():
            path = run / translated(rel, run, pairs)
            if rel not in modified:
                assert stat_record(path) == {k: v for k, v in before.items() if k != "sha256"}, rel
            if "sha256" in before:
                assert sha256(path) == before["sha256"], rel
                preserved += 1
        check_artifact_links(run)
        # Original run times/commands and .snakemake records remain unchanged; the new
        # launcher refuses this historical run rather than claiming the new code ran it.
        assert load(run / "run_status.json") == journal["analysis_status"]
        journal.update(
            status="MIGRATED",
            finished_utc=datetime.now(timezone.utc).isoformat(),
            metadata_updated=modified,
            protected_hashes_verified=preserved,
            original_snakemake_metadata="preserved unchanged; historical layout-1 evidence",
            reuse_policy="historical completed run; new computation requires a new run directory",
        )
        save(status_path, journal)
        return journal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    result = migrate(args.run_dir, args.apply)
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in {"files_before", "metadata_before"}},
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
