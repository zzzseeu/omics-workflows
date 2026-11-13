"""Workflow-independent configuration, provenance and local execution utilities.

Callers supply stage paths, tools and Python validation functions. This module
contains no assay-specific configuration keys or analysis routes.
"""

import ast
import copy
import csv
import fcntl
import hashlib
import json
import math
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
    """将可序列化对象转换为排序、缩进固定且以换行结束的JSON字符串，供记录和指纹计算使用。"""
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n"


def digest(value):
    """对对象的规范JSON表示计算SHA-256，返回十六进制字符串，作为计划或参数的稳定标识。"""
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def write_stable(path, value):
    """将对象写入JSON路径；内容相同时保留原文件时间，否则先写临时文件再替换，避免无关重算。"""
    path = Path(path)
    content = canonical(value)
    if path.exists() and path.read_text(encoding="utf-8") == content:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(content, encoding="utf-8")
    tmp.replace(path)


def merge(base, override):
    """递归合并两个参数字典，返回独立副本；override覆盖base，嵌套映射逐项合并而不修改输入。"""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def identifier(value, label):
    """校验非空ASCII标识符，仅允许字母、数字、下划线、点和连字符；返回原值，非法时抛ConfigError。"""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value or ""):
        raise ConfigError(f"{label}: use a nonempty ASCII identifier (letters, digits, _, ., -)")
    return value


def table(path, required, allowed):
    """读取TSV为字典列表，并检查必需列、允许列、重复列和行完整性；异常转换为ConfigError。"""
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
    """返回文件大小、纳秒时间和inode组成的轻量签名。

    不读取文件内容；省略跨计算节点可能变化的设备编号，签名本身不替代内容校验值。
    """
    stat = Path(path).stat()
    # st_dev identifies a local mount and changes between HPC clients of the
    # same shared filesystem. It must not invalidate a cross-node resume.
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "inode": stat.st_ino,
    }


def check_fields(value: dict, allowed=None, required=(), path="config") -> None:
    """检查配置映射、字符串键、必填项和未知字段；allowed为None时允许模块自定义名称。"""
    if not isinstance(value, dict):
        raise ConfigError(f"{path} must be a mapping")
    if any(not isinstance(key, str) for key in value):
        raise ConfigError(f"{path} keys must be strings")
    missing = set(required) - value.keys()
    if missing:
        raise ConfigError(f"{path} missing required fields: {', '.join(sorted(missing))}")
    unknown = value.keys() - set(allowed) if allowed is not None else set()
    if unknown:
        raise ConfigError(f"{path} unknown fields: {', '.join(sorted(unknown))}")


def check_value(
    value, expected_type, path, *, minimum=None, maximum=None, choices=None, pattern=None, min_length=None
) -> None:
    """在Python中检查一个值的实际类型及约束，失败时用完整字段路径报错。

    expected_type为类型或类型元组；不会把字符串、小数或布尔值隐式转换成整数。
    数值须有限，范围包含边界；字符串可检查长度和完整正则匹配。
    """
    types = expected_type if isinstance(expected_type, tuple) else (expected_type,)
    valid_type = isinstance(value, types) and (not isinstance(value, bool) or bool in types)
    if not valid_type:
        expected = " or ".join(kind.__name__ for kind in types)
        raise ConfigError(f"{path} must be {expected}, got {type(value).__name__}")
    if value is not None:
        if isinstance(value, float) and not math.isfinite(value):
            raise ConfigError(f"{path} must be finite")
        if minimum is not None and value < minimum:
            raise ConfigError(f"{path} must be >= {minimum}")
        if maximum is not None and value > maximum:
            raise ConfigError(f"{path} must be <= {maximum}")
        if min_length is not None and len(value) < min_length:
            raise ConfigError(f"{path} must contain at least {min_length} characters")
        if pattern is not None and re.fullmatch(pattern, value) is None:
            raise ConfigError(f"{path} has an invalid format")
    if choices is not None and value not in choices:
        raise ConfigError(f"{path} must be one of {choices!r}")


def load_configuration(path, validate, normalize=None):
    """读取YAML并调用模块的Python校验函数，返回原配置字典，不读取独立校验文件。

    validate(config)检查成功时返回None，非法配置抛ConfigError；normalize用于兼容旧字段。
    默认参数合并与跨字段分析条件由模块负责，此函数不修改参数类型或默认值。
    """
    import yaml

    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if normalize:
            data = normalize(data)
        check_fields(data)
        validate(data)
    except (OSError, UnicodeError, yaml.YAMLError, ConfigError) as error:
        raise ConfigError(f"Configuration validation failed: {path}: {error}") from error
    return data


def load_defaults(path):
    """读取默认参数YAML路径，返回解析后的对象；不在这里混入样本覆盖设置。"""
    import yaml

    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def validate_run_directory(
    repository, run_dir, cores, mem_mb, layout_version, historical_marker=None
):
    """校验运行目录、资源预算与输出布局版本，失败时抛ConfigError。

    禁止把结果写入代码目录；历史完成结果或不兼容布局要求使用新运行目录，不移动或覆盖已有文件。
    """
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
            raise ConfigError(
                "Legacy output layout detected. Choose a new --run-dir. "
                "No legacy outputs or Snakemake metadata were overwritten."
            )


def select_samples(rows, selected=None, overrides=()):
    """校验样本ID、重复项和覆盖设置，返回样本索引及所选样本顺序；未指定时选择全部样本。"""
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
    """按样本归类测序数据行，检查样本引用、数据组ID和调用者允许的技术类型；返回分组字典。"""
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
    """追踪入口脚本及其包内Python依赖，返回排序后的路径列表。

    通过AST读取导入关系，不执行被分析代码；限制在给定代码根目录内，不扫描第三方库。
    """
    repository = Path(repository).resolve()
    pending = [repository / entry]
    found = set()

    def module_path(name):
        """将绝对导入名映射到代码根目录内的路径，兼容普通目录与统一包命名空间。"""
        parts = name.split(".") if name else []
        if parts[:1] == [repository.name] and (repository / "__init__.py").is_file():
            parts = parts[1:]
        return repository.joinpath(*parts)

    while pending:
        path = pending.pop().resolve()
        if path in found:
            continue
        found.add(path)
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom):
                if node.level:
                    base = path.parents[node.level - 1]
                    module = base.joinpath(*(node.module or "").split("."))
                else:
                    module = module_path(node.module)
                candidates = [module.with_suffix(".py"), module / "__init__.py"]
                candidates.extend(module / (alias.name + ".py") for alias in node.names)
                candidates.extend(
                    parent / "__init__.py" for parent in module.parents if parent.is_relative_to(repository)
                )
            elif isinstance(node, ast.Import):
                candidates = []
                for alias in node.names:
                    module = module_path(alias.name)
                    candidates.extend([module.with_suffix(".py"), module / "__init__.py"])
            else:
                continue
            for candidate in candidates:
                if candidate.is_file() and candidate.resolve().is_relative_to(repository.resolve()):
                    pending.append(candidate)
    return sorted(found)


def stage_sources(repository, rule, entries):
    """汇总步骤规则及其入口依赖的SHA-256，返回以包内相对路径为键的来源字典。"""
    repository = Path(repository).resolve()
    paths = {repository / rule}
    for entry in entries:
        paths.update(script_sources(repository, entry))
    return {
        str(path.relative_to(repository)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(paths)
    }


def prepend_search_paths(executables=(), modules=()):
    """将存在的可执行文件或Python模块目录加入进程搜索路径，并按首次出现顺序去重。"""
    for key, paths in (("PATH", executables), ("PYTHONPATH", modules)):
        additions = [str(path) for path in paths if Path(path).is_dir()]
        existing = os.environ.get(key, "").split(os.pathsep)
        os.environ[key] = os.pathsep.join(dict.fromkeys([*additions, *existing]))


def inspect_tools(required, prefix, versions=None, fallback=None):
    """检查指定工具是否位于目标环境中，返回工具记录字典和缺失名称列表。

    记录真实路径和文件SHA-256；仅对versions中指定的工具调用版本命令，fallback可补充定位方式。
    """
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
    """处理缺失依赖列表：严格模式抛ConfigError，预演模式打印缺项；列表为空时不做操作。"""
    if missing:
        message = "Required dependencies/inputs missing: " + ", ".join(missing)
        if strict:
            raise ConfigError(message)
        print("DRY-RUN ONLY: " + message, flush=True)


def file_checksum(path):
    """分块读取文件并返回SHA-256十六进制值，避免将大文件一次载入内存。"""
    digest_value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest_value.update(chunk)
    return digest_value.hexdigest()


def verify_source_locks(tools, specifications):
    """仅校验本次使用到的源码扩展锁，并把锁文件指纹补入工具记录。

    锁目录、文件清单键名和相对根目录均由模块提供；内容不一致时拒绝运行。
    """
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
    """检查已完成步骤的结果清单和文件，并调用模块提供的专用校验器。

    轻量签名变化时再核对内容指纹；发现错样本、文件丢失或内容变化时拒绝复用缓存。
    """
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
                    # 签名变化后再核对完整内容，兼顾大文件检查成本与来源可靠性。
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
    """按步骤和样本写入参数JSON；输入stages定义规则、入口和工具范围。

    build_metadata负责生成模块专属内容；仅写入与该步骤有关的软件和代码，减少无关重算。
    """
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
    """根据运行计划启动本地Snakemake，并返回退出码。

    plan至少包含run_dir、sample_order、plan_id、target和CPU/内存budget；
    prepare在取得运行锁后执行。记录计划与状态，将中断传给整个子进程组，失败不得记为完成。
    """
    snakemake = Path(sys.prefix) / "bin/snakemake"
    run_dir = Path(plan["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / ".launcher.lock").open("a", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ConfigError(f"Run directory already has an active launcher: {run_dir}") from None
        # 先取得运行锁，再生成参数与检查缓存，避免并发启动器相互覆盖。
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
                """记录收到的中断信号并转发到分析进程组；子进程已退出时忽略查找失败。"""
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
