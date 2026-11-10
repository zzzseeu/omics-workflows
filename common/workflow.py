"""Workflow-independent configuration, provenance and local execution utilities.

Callers supply schemas, stage paths, tools and validation functions. This module
contains no assay-specific configuration keys or analysis routes.
"""

import ast
import copy
import csv
import fcntl
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


class ConfigError(ValueError):
    pass


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n"


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def write_stable(path, value):
    path = Path(path)
    content = canonical(value)
    if path.exists() and path.read_text(encoding="utf-8") == content:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(content, encoding="utf-8")
    tmp.replace(path)


def merge(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def identifier(value, label):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value or ""):
        raise ConfigError(f"{label}: use a nonempty ASCII identifier (letters, digits, _, ., -)")
    return value


def table(path, required, allowed):
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)) or not required <= set(fields):
                raise ConfigError(f"{path}: required unique columns: {sorted(required)}")
            if set(fields) - allowed:
                raise ConfigError(f"{path}: unknown columns: {sorted(set(fields) - allowed)}")
            result = []
            for line, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise ConfigError(f"{path}:{line}: malformed TSV row")
                if not any(v.strip() for v in row.values()):
                    continue
                result.append({k: v.strip() for k, v in row.items()})
            return result
    except OSError as exc:
        raise ConfigError(str(exc)) from exc


def source_signature(path):
    stat = Path(path).stat()
    # st_dev identifies a local mount and changes between HPC clients of the
    # same shared filesystem. It must not invalidate a cross-node resume.
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "inode": stat.st_ino,
    }


def load_configuration(path, schema_path, normalize=None):
    """Load and validate configuration using a caller-provided JSON schema."""
    import jsonschema
    import yaml

    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if normalize:
            data = normalize(data)
        schema = yaml.safe_load(Path(schema_path).read_text(encoding="utf-8"))
        jsonschema.Draft202012Validator(schema).validate(data)
    except (OSError, yaml.YAMLError, jsonschema.ValidationError) as error:
        raise ConfigError(f"Configuration validation failed: {error}") from error
    return data


def load_defaults(path):
    import yaml

    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def validate_run_directory(
    repository, run_dir, cores, mem_mb, layout_version, historical_marker=None, migration_tool=None
):
    repository, run_dir = Path(repository).resolve(), Path(run_dir).resolve()
    if cores < 1 or mem_mb < 1024:
        raise ConfigError("Execution budget requires cores >= 1 and mem_mb >= 1024")
    if run_dir == repository or run_dir.is_relative_to(repository):
        raise ConfigError("Run directory must be outside the workflow code repository")
    if historical_marker and (run_dir / historical_marker).exists():
        raise ConfigError(
            "This completed run was migrated and is preserved as a historical result. "
            "Read its numbered outputs; use a new --run-dir for new computation."
        )
    previous = run_dir / "config/plan.json"
    if previous.is_file():
        old = json.loads(previous.read_text(encoding="utf-8"))
        if old.get("layout_version", 1) != layout_version and (run_dir / "samples").is_dir():
            advice = (
                f"Use {migration_tool} for historical runs, or choose a new --run-dir. "
                if migration_tool
                else "Choose a new --run-dir. "
            )
            raise ConfigError(
                "Legacy output layout detected. "
                + advice
                + "No legacy outputs or Snakemake metadata were overwritten."
            )


def select_samples(rows, selected=None, overrides=()):
    samples = {}
    for row in rows:
        sid = identifier(row["sample_id"], "sample_id")
        if sid in samples:
            raise ConfigError(f"Duplicate sample_id: {sid}")
        samples[sid] = row
    if not samples:
        raise ConfigError("samples.tsv is empty")
    unknown = set(overrides) - set(samples)
    if unknown:
        raise ConfigError(f"Sample overrides refer to unknown samples: {sorted(unknown)}")
    selected = selected or list(samples)
    if len(selected) != len(set(selected)) or set(selected) - set(samples):
        raise ConfigError(f"Duplicate or unknown selected samples: {selected}")
    return samples, selected


def group_readsets(rows, samples, technologies):
    groups, keys = {}, set()
    for row in rows:
        sid = row["sample_id"]
        if sid not in samples:
            raise ConfigError(f"Readset refers to unknown sample: {sid}")
        rid = identifier(row["readset_id"], "readset_id")
        identifier(row["library_id"], "library_id")
        if (sid, rid) in keys:
            raise ConfigError(f"Duplicate readset_id for {sid}: {rid}")
        keys.add((sid, rid))
        if row["technology"] not in technologies:
            raise ConfigError(f"{sid}/{rid}: unknown technology")
        groups.setdefault(sid, []).append(row)
    return groups


def script_sources(repository, entry):
    """Follow local Python imports so changes invalidate only dependent steps."""
    repository = Path(repository).resolve()
    pending = [repository / entry]
    found = set()
    while pending:
        path = pending.pop().resolve()
        if path in found:
            continue
        found.add(path)
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom):
                if node.level:
                    base = path.parents[node.level - 1]
                else:
                    base = repository
                module = base.joinpath(*(node.module or "").split("."))
                candidates = [module.with_suffix(".py"), module / "__init__.py"]
                candidates.extend(module / (alias.name + ".py") for alias in node.names)
                candidates.extend(
                    parent / "__init__.py" for parent in module.parents if parent.is_relative_to(repository)
                )
            elif isinstance(node, ast.Import):
                candidates = [
                    repository.joinpath(*alias.name.split(".")).with_suffix(".py") for alias in node.names
                ]
            else:
                continue
            for candidate in candidates:
                if candidate.is_file() and candidate.resolve().is_relative_to(repository.resolve()):
                    pending.append(candidate)
    return sorted(found)


def stage_sources(repository, rule, entries):
    repository = Path(repository).resolve()
    paths = {repository / rule}
    for entry in entries:
        paths.update(script_sources(repository, entry))
    return {
        str(path.relative_to(repository)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(paths)
    }


def prepend_search_paths(executables=(), modules=()):
    for key, paths in (("PATH", executables), ("PYTHONPATH", modules)):
        additions = [str(path) for path in paths if Path(path).is_dir()]
        existing = os.environ.get(key, "").split(os.pathsep)
        os.environ[key] = os.pathsep.join(dict.fromkeys([*additions, *existing]))


def inspect_tools(required, prefix, versions=None, fallback=None):
    """Locate tools inside the supplied environment and record their identities."""
    prefix = Path(prefix).resolve()
    tools, missing = {}, []
    for name in sorted(required):
        executable = shutil.which(name)
        if not executable and fallback:
            executable = fallback(name)
        if (
            not executable
            or not os.access(executable, os.X_OK)
            or not Path(executable).resolve().is_relative_to(prefix)
        ):
            missing.append(name)
            tools[name] = {"path": executable, "version": "MISSING"}
            continue
        path = Path(executable).resolve()
        tools[name] = {"path": str(path), "sha256": file_checksum(path)}
        if name in (versions or {}):
            version = subprocess.run(
                [str(path), *versions[name]], capture_output=True, timeout=30, check=False
            )
            lines = (version.stdout + version.stderr).decode(errors="replace").splitlines()
            tools[name]["version"] = lines[0] if lines else "UNREPORTED"
    return tools, missing


def require_available(missing, strict=True):
    if missing:
        message = "Required dependencies/inputs missing: " + ", ".join(missing)
        if strict:
            raise ConfigError(message)
        print("DRY-RUN ONLY: " + message, flush=True)


def file_checksum(path):
    digest_value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest_value.update(chunk)
    return digest_value.hexdigest()


def verify_source_locks(tools, specifications):
    """Validate caller-specified source locks only for extensions in use."""
    for spec in specifications:
        root = Path(spec["root"]).resolve()
        users = [
            name
            for name, record in tools.items()
            if record.get("path") and Path(record["path"]).resolve().is_relative_to(root)
        ]
        if not users:
            continue
        lock = root / "SOURCE_LOCK.json"
        if not lock.is_file():
            raise ConfigError(f"Missing installed source lock: {lock}")
        data = json.loads(lock.read_text(encoding="utf-8"))
        for relative, expected in data[spec["files_key"]].items():
            path = root / spec.get("files_base", "") / relative
            if not path.is_file() or file_checksum(path) != expected:
                raise ConfigError(f'Locked {spec["label"]} source changed: {relative}')
        value = file_checksum(lock)
        for name in users:
            tools[name]["source_lock_sha256"] = value


def verify_cached_outputs(plan, rules):
    """Validate timestamped results with caller-supplied manifest validators."""
    for sid in plan["sample_order"]:
        base = Path(plan["run_dir"]) / "samples" / sid
        for stage, (filename, validator) in rules.items():
            folder = base / stage
            if not (folder / ".snakemake_timestamp").is_file():
                continue
            manifest = folder / filename
            try:
                data = json.loads(manifest.read_text(encoding="utf-8"))
                if data["sample_id"] != sid:
                    raise ValueError("Cached result belongs to another sample")
                for record in data.get("files", []):
                    path = folder / record["file"]
                    if not path.is_file():
                        raise ValueError(f"Missing declared artifact: {path}")
                    if (
                        record.get("signature") != source_signature(path)
                        and file_checksum(path) != record["sha256"]
                    ):
                        raise ValueError(f"Artifact checksum mismatch: {path}")
                if validator:
                    validator(manifest, sid)
            except (OSError, ValueError, KeyError) as error:
                raise ConfigError(f"Invalid cached output in {folder}: {error}; use a new run-dir") from error


def write_plan_metadata(plan, repository, run_dir, stages, build_metadata):
    for stage, definition in stages.items():
        code = stage_sources(repository, definition["rule"], definition["entries"])
        names = definition["tools"]
        tools = {key: value for key, value in plan["tools"].items() if names is None or key in names}
        for sid, sample in plan["samples"].items():
            write_stable(
                Path(run_dir) / f"config/samples/{sid}.{stage}.json",
                build_metadata(sample, stage, tools, code),
            )


def execute_workflow(plan, snakefile, target, prepare=None, dry_run=False):
    """Launch one local Snakemake workflow with locking and signal propagation."""
    snakemake = Path(sys.prefix) / "bin/snakemake"
    run_dir = Path(plan["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / ".launcher.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ConfigError(f"Run directory already has an active launcher: {run_dir}") from None
        if prepare:
            prepare()
        write_stable(run_dir / "config/plan.json", plan)
        command = [
            str(snakemake),
            "--snakefile",
            str(snakefile),
            "--directory",
            str(run_dir),
            "--configfile",
            str(run_dir / "config/plan.json"),
            "--executor",
            "local",
            "--cores",
            str(plan["budget"]["cores"]),
            "--resources",
            f'mem_mb={plan["budget"]["mem_mb"]}',
            "--rerun-incomplete",
            "--printshellcmds",
            target,
        ]
        if dry_run:
            command.append("--dry-run")
        status = {
            "plan_id": plan["plan_id"],
            "samples": plan["sample_order"],
            "target": plan["target"],
            "executor": "local",
            "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(),
        }
        status_path = run_dir / ("dryrun_status.json" if dry_run else "run_status.json")
        write_stable(status_path, {**status, "status": "PLANNING" if dry_run else "RUNNING"})
        env = dict(
            os.environ,
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            NUMEXPR_NUM_THREADS="1",
            PYTHONNOUSERSITE="1",
            PYTHONDONTWRITEBYTECODE="1",
        )
        handlers, interrupted = {}, []
        try:
            # Explicit process cleanup below owns the pipe/signal lifetime.
            # pylint: disable-next=consider-using-with
            process = subprocess.Popen(command, env=env, start_new_session=True)

            def forward_signal(signum, _frame):
                interrupted.append(signum)
                try:
                    os.killpg(process.pid, signum)
                except ProcessLookupError:
                    pass

            for signum in (signal.SIGINT, signal.SIGTERM):
                handlers[signum] = signal.signal(signum, forward_signal)
            interruption_started = None
            while True:
                try:
                    rc = process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    if interrupted:
                        if interruption_started is None:
                            interruption_started = time.monotonic()
                        elif time.monotonic() - interruption_started > 15:
                            # Some Snakemake worker shutdown races otherwise leave
                            # the outer allocation and launcher lock alive forever.
                            import psutil

                            try:
                                descendants = psutil.Process(process.pid).children(recursive=True)
                            except psutil.NoSuchProcess:
                                descendants = []
                            for child in reversed(descendants):
                                try:
                                    child.kill()
                                except psutil.NoSuchProcess:
                                    pass
                            try:
                                os.killpg(process.pid, signal.SIGKILL)
                            except ProcessLookupError:
                                pass
                            rc = process.wait()
                            break
            if interrupted:
                rc = 128 + interrupted[-1]
        except OSError:
            rc = 1
            raise
        finally:
            for signum, handler in handlers.items():
                signal.signal(signum, handler)
            write_stable(
                status_path,
                {
                    **status,
                    "exit_code": locals().get("rc", 1),
                    "status": (
                        ("PLANNED" if dry_run else "COMPLETED") if locals().get("rc", 1) == 0 else "FAILED"
                    ),
                    "finished_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
        return rc
