# 开发与维护

框架代码、规则和模板随Python包发布；样本数据、项目配置、CSUB脚本、日志和结果由各项目保存。新增流程使用现有Snakemake local executor，集群资源由项目在外部申请。

## 1. 接入前先确定分析接口

先写清分析目标和交付，再编写规则。每个模块至少确定以下内容：

| 内容 | 需要明确的问题 |
|---|---|
| 输入 | 接受什么数据和格式；需要哪些样本信息、参考序列或数据库；相对路径以哪里为基准 |
| 方法 | 每步的生物学目标、工具、参数和适用条件；哪些步骤可选 |
| 输出 | 核心文件、格式、指标定义和来源清单；每步目录编号与规则、脚本一致 |
| 校验 | 错样本、文件损坏、缺数据、工具失败如何处理；何时允许复用结果 |
| 资源 | 各步骤线程和内存需求；是否存在并发子进程 |
| 验收 | 软件执行完成、质量检查完成和生物学结论成立分别需要什么证据 |

模块不必使用样本表或七个步骤。组装、转录组、变异检测各自定义数据要求；公共模块只处理运行共性。

## 2. 准备目录和公共接口

在 `src/omics_workflows/workflows/<name>/` 新建模块。以下 `seqstats` 是接入教程示例，当前正式注册的流程仍为assembly。它使用seqkit统计一个序列文件，只有一个步骤。

```text
src/omics_workflows/workflows/seqstats/
├── __init__.py
├── runner.py
├── Snakefile
├── defaults.yaml
├── rules/01_stats.smk
├── scripts/
│   ├── __init__.py
│   └── 01_stats.py
└── examples/config.yaml
```

Python导入目录添加 `__init__.py`；规则和示例数据子目录作为包资源。文件定位相对安装包，不能写死服务器、项目目录或环境路径。实际数据路径在预检时解析为绝对路径，写入来源记录。

公共能力位于 [common/workflow.py](../src/omics_workflows/common/workflow.py)、[resources.py](../src/omics_workflows/common/resources.py)和[arguments.py](../src/omics_workflows/common/arguments.py)：

| 接口 | 输入与返回/作用 |
|---|---|
| `argument_parser(description, prog, epilog)` | 建立统一帮助格式；调用者逐参数填写类型、含义、单位及默认行为 |
| `package_directory()` | 在上下文中提供安装包根目录；生命周期须覆盖Snakemake执行 |
| `load_configuration(path, validate, normalize=None)` | 读取YAML，先按需兼容旧字段，再调用模块的Python校验函数；返回配置字典 |
| `merge(base, override)` | 递归合并副本，override优先，不修改原字典 |
| `validate_run_directory(...)` | 检查代码边界、资源和layout_version，不自动迁移历史结果 |
| `inspect_tools(required, prefix, versions)` | 返回工具记录和缺失列表；只接受指定环境内工具 |
| `stage_sources(root, rule, entries)` | 记录规则、入口及其包内依赖的SHA-256 |
| `write_plan_metadata(plan, root, run_dir, stages, builder)` | 按步骤写JSON；builder接收sample、stage、tools、code并返回模块参数 |
| `verify_cached_outputs(plan, rules)` | rules将编号目录映射到清单名及可选校验器；校验器接收清单路径和样本ID |
| `execute_workflow(plan, snakefile, target, prepare, dry_run)` | 加锁、执行prepare、写计划、启动Snakemake、记录状态并返回退出码 |

参数检查直接写在模块的 `validate_config(config: dict) -> None` 函数中：检查通过返回None，非法配置抛ConfigError。公共函数 `check_fields` 可检查必填项和未知键，`check_value` 可检查实际类型、范围及可选值；也可以使用普通Python条件语句编写专用检查。类型注解本身不执行校验。默认值继续保存在defaults.yaml，模块在读取并检查用户配置后合并默认值。

`execute_workflow` 必需的plan字段：`run_dir`、`sample_order`、`plan_id`、`target`、`budget.cores`、`budget.mem_mb`。规则通常还使用 `samples`、`tools`；建议同时记录 `layout_version`、配置来源、框架版本和代码指纹。`plan_id` 在计划填充完成后计算。

STAGES的键是参数文件中的步骤名，如stats；`rule`和`entries`是相对包根目录的路径；`tools`是该步骤使用的软件集合，空集合表示无外部工具，None表示记录全部工具。它与输出目录编号是不同概念，例如stats对应 `01_stats`。

## 3. 实现一个可运行的最小模块

下面各代码块对应完整文件。先创建两个空 `__init__.py`。示例在runner.py中定义参数检查，不需要额外的校验描述文件。本例使用当前激活的分析环境，提供与现有流程一致的配置、资源、预检和预演参数。

### `runner.py`

runner负责检查输入、构建计划和启动分析，具体统计交给步骤脚本。新增模块的入口必须是 `main(argv=None, prog=None)` 并返回退出码。

```python
"""序列统计模块：配置预检、来源记录及步骤启动。"""

import sys
from pathlib import Path

from omics_workflows import __version__
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.resources import package_directory, source_boundary
from omics_workflows.common.workflow import (
    ConfigError,
    check_fields,
    check_value,
    digest,
    execute_workflow,
    identifier,
    inspect_tools,
    load_configuration,
    load_defaults,
    merge,
    require_available,
    source_signature,
    stage_sources,
    validate_run_directory,
    verify_cached_outputs,
    write_plan_metadata,
)

BASE = "workflows/seqstats"
ENTRY = f"{BASE}/scripts/01_stats.py"
STAGES = {"stats": {"rule": f"{BASE}/rules/01_stats.smk", "entries": [ENTRY], "tools": {"seqkit"}}}


def validate_config(config: dict) -> None:
    """直接检查字段、类型和取值；成功返回None，错误抛ConfigError，不修改配置。"""
    required = {"sample_id", "sequence_file"}
    check_fields(config, required | {"threads"}, required=required)
    for name in ("sample_id", "sequence_file"):
        check_value(config[name], str, name, min_length=1)
    if "threads" in config:
        check_value(config["threads"], int, "threads", minimum=1)


def build_metadata(sample, _stage, tools, code):
    """生成一个样本的统计参数及来源记录，供步骤脚本读取。"""
    return {**sample, "tools": tools, "code": code}


def run_analysis(args):
    """校验配置并生成计划，返回预检或Snakemake执行的退出码。"""
    with package_directory() as package:
        boundary = source_boundary(package)
        run_dir = args.run_dir.resolve()
        validate_run_directory(boundary, run_dir, args.cores, args.mem_mb, 1)
        config_path = args.config.resolve()
        raw = load_configuration(config_path, validate_config)
        settings = merge(load_defaults(package / BASE / "defaults.yaml"), raw)
        sid = identifier(settings["sample_id"], "sample_id")
        source = (config_path.parent / settings["sequence_file"]).resolve()
        if not source.is_file() or source.is_relative_to(run_dir) or source.is_relative_to(boundary):
            raise ConfigError("sequence_file must be an existing file outside code and run directories")
        tools, missing = inspect_tools(
            {"snakemake", "seqkit"}, Path(sys.prefix), versions={"seqkit": ["version"]}
        )
        require_available(missing, strict=not args.dry_run or args.check)
        sample = {
            "sample_id": sid,
            "sequence_file": str(source),
            "input_signature": source_signature(source),
            "threads": min(settings["threads"], args.cores),
        }
        plan = {
            "layout_version": 1,
            "run_dir": str(run_dir),
            "config_path": str(config_path),
            "sample_order": [sid],
            "samples": {sid: sample},
            "tools": tools,
            "framework": {"version": __version__, "workflow": "seqstats"},
            "budget": {"cores": args.cores, "mem_mb": args.mem_mb},
            "target": "stats",
            "code": stage_sources(package, f"{BASE}/rules/01_stats.smk", [ENTRY, f"{BASE}/runner.py"]),
        }
        plan["plan_id"] = digest(plan)
        print(f"PRECHECK PASS: sample={sid}", flush=True)
        if args.check:
            return 0

        def prepare():
            """取得运行锁后核验缓存，并写入各样本的步骤参数。"""
            verify_cached_outputs(plan, {"01_stats": ("manifest.json", None)})
            write_plan_metadata(plan, package, run_dir, STAGES, build_metadata)

        return execute_workflow(
            plan, package / BASE / "Snakefile", "seqstats_result", prepare=prepare, dry_run=args.dry_run
        )


def main(argv=None, prog=None):
    """提供统一模块入口；分析依赖使用当前激活的环境。"""
    parser = argument_parser("统计单个FASTA/FASTQ文件的序列数量和长度。", prog=prog)
    parser.add_argument("--config", type=Path, required=True, help="类型：YAML路径；必填。序列统计配置文件。")
    parser.add_argument(
        "--run-dir", type=Path, required=True, help="类型：目录路径；必填。代码目录之外的结果目录。"
    )
    parser.add_argument("--cores", type=int, required=True, help="类型：正整数；必填。流程可用CPU总数。")
    parser.add_argument(
        "--mem-mb", type=int, required=True, help="类型：整数，单位MB，至少1024；必填。流程内存预算。"
    )
    parser.add_argument(
        "--check", action="store_true", help="类型：布尔开关；默认关闭。仅检查配置和软件，不创建结果目录。"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="类型：布尔开关；默认关闭。写运行计划并预演，不执行统计。"
    )
    args = parser.parse_args(argv)
    try:
        return run_analysis(args)
    except (ConfigError, OSError, ValueError) as error:
        parser.exit(2, f"PRECHECK FAILED: {error}\n")
        return 2
```

### `defaults.yaml`

只保存默认值；字段类型和取值规则在上面的validate_config函数中维护。

```yaml
threads: 1  # 正整数；默认统计线程数，运行时不超过--cores。
```

### `Snakefile`

Snakefile定义分析终点并包含编号规则。脚本和依赖源码均作为规则输入，代码变化才能触发相应步骤更新。

```python
"""序列统计示例的分析终点和规则引用。"""
from pathlib import Path
import sys
ROOT = Path(workflow.basedir).parents[1]
sys.path.insert(0, str(ROOT.parent))
from omics_workflows.common.workflow import script_sources
PYTHON = sys.executable
ENTRY = ROOT / "workflows/seqstats/scripts/01_stats.py"
HELPERS = [str(path) for path in script_sources(ROOT, ENTRY) if path != ENTRY]

rule seqstats_result:
    input: expand("samples/{sample}/01_stats", sample=config["sample_order"])

include: "rules/01_stats.smk"
```

### `rules/01_stats.smk`

规则声明输入、输出目录、线程、内存和日志。外部命令参数使用Snakemake的`:q`处理路径；分析脚本内部用参数列表执行工具。

```python
rule sequence_statistics:
    input:
        source=lambda w: config["samples"][w.sample]["sequence_file"],
        metadata="config/samples/{sample}.stats.json",
        code=str(ENTRY),
        helpers=HELPERS
    output: directory("samples/{sample}/01_stats")
    threads: lambda w: config["samples"][w.sample]["threads"]
    resources: mem_mb=1024
    log: "samples/{sample}/logs/01_stats.log"
    shell:
        "{PYTHON:q} {input.code:q} --metadata {input.metadata:q} "
        "--out-dir {output:q} > {log:q} 2>&1"
```

### `scripts/01_stats.py`

步骤脚本只处理一个明确目标。先核验输入身份，在临时文件中运行工具，确认成功后发布结果；异常保持非零退出。

```python
"""调用seqkit统计单个序列文件并交付统计表和来源清单。"""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

# 从脚本自身的包内位置定位导入目录，安装到任意位置均可运行。
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
# pylint: disable=wrong-import-position
from omics_workflows.common.arguments import argument_parser
from omics_workflows.common.workflow import file_checksum, source_signature, write_stable


def main():
    """读取步骤参数，核验输入身份，成功后发布TSV及manifest.json。"""
    parser = argument_parser("统计FASTA/FASTQ序列数量和长度。")
    parser.add_argument(
        "--metadata", type=Path, required=True, help="类型：JSON路径；必填。样本、工具和输入来源记录。"
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True, help="类型：目录路径；必填。本步骤的输出目录。"
    )
    args = parser.parse_args()
    meta = json.loads(args.metadata.read_text(encoding="utf-8"))
    source = Path(meta["sequence_file"])
    before = source_signature(source)
    if before != meta["input_signature"]:
        raise ValueError("Input changed after precheck")
    input_sha256 = file_checksum(source)
    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()) and not (output / "manifest.json").is_file():
        raise ValueError("Refusing to overwrite an unmanaged nonempty directory")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=output, prefix=".stats-", delete=False
        ) as handle:
            temporary = Path(handle.name)
            subprocess.run(
                [meta["tools"]["seqkit"]["path"], "stats", "-T", "-j", str(meta["threads"]), str(source)],
                stdout=handle,
                check=True,
            )
        with temporary.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) != 1 or int(rows[0]["num_seqs"].replace(",", "")) < 1:
            raise ValueError("Expected statistics for one nonempty sequence file")
        if source_signature(source) != before:
            raise ValueError("Input changed during statistics")
        stats = output / "statistics.tsv"
        temporary.replace(stats)
        write_stable(
            output / "manifest.json",
            {
                "sample_id": meta["sample_id"],
                "status": "COMPLETED",
                "input_sha256": input_sha256,
                "tools": meta["tools"],
                "code": meta["code"],
                "files": [
                    {
                        "file": "statistics.tsv",
                        "sha256": file_checksum(stats),
                        "signature": source_signature(stats),
                    }
                ],
            },
        )
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
```

### `examples/config.yaml`

真实项目可按该配置填写输入，示例数据保存在项目目录的data下。每个key后写清类型及含义。

```yaml
sample_id: rice01  # 字符串；样本ID。
sequence_file: ../data/rice01.fa  # 字符串路径；相对本YAML所在目录，支持FASTA/FASTQ及其压缩文件。
threads: 2  # 正整数；按实际资源填写。
```

## 4. 注册命令、资源和依赖

在 [cli.py](../src/omics_workflows/cli.py) 的WORKFLOWS中增加一项，保留已有assembly注册：

```python
"seqstats": ("omics_workflows.workflows.seqstats.runner", "序列数量与长度统计"),
```

在 [pyproject.toml](../pyproject.toml) 增加模块资源，避免源码能运行而wheel缺文件：

```toml
"omics_workflows.workflows.seqstats" = [
    "Snakefile", "*.yaml", "rules/*.smk", "examples/*.yaml",
]
```

该条目应加入现有 `[tool.setuptools.package-data]` 表。Python包由src自动发现；运行资源须显式声明。模块新增的直接Python依赖写入 `project.optional-dependencies.<name>`，通用依赖写入project.dependencies；seqkit等生信软件列入 [environment.yaml](../env/environment.yaml) 并保存实际部署锁，不在分析中自动安装软件或数据库。本例所需seqkit已在统一环境需求中。

为每个模块添加一份 `docs/workflows/<name>.md`，说明输入、方法、输出格式、质量检查和运行示例，并更新[文档导航](README.md)。[Assembly手册](workflows/assembly.md)仅描述组装方法，不作为公共接口的定义。

## 5. 安装和验证新模块

在开发环境中安装框架；模块使用的生信工具须已部署到当前环境：

```bash
python -m pip install -e '.[assembly,dev]'
omics-workflows list
omics-workflows seqstats --help
```

在仓库外准备project/config/seqstats.yaml和project/data/rice01.fa。进入project目录后运行：

```bash
omics-workflows seqstats --config config/seqstats.yaml \
  --run-dir runs/seqstats001 --cores 2 --mem-mb 4096 --check
omics-workflows seqstats --config config/seqstats.yaml \
  --run-dir runs/seqstats001 --cores 2 --mem-mb 4096 --dry-run
omics-workflows seqstats --config config/seqstats.yaml \
  --run-dir runs/seqstats001 --cores 2 --mem-mb 4096
```

本站实质计算和合成数据生成使用CSUB；上述实际运行命令放入项目作业，申请资源与命令预算一致。本例输出：

```text
runs/seqstats001/
├── config/plan.json
├── config/samples/rice01.stats.json
├── samples/rice01/01_stats/
│   ├── statistics.tsv
│   └── manifest.json
├── samples/rice01/logs/01_stats.log
└── run_status.json
```

还须执行以下验证，不能只确认--help可启动：

| 检查 | 应验证的行为 |
|---|---|
| 配置与输入 | 拒绝未知键、非法类型、重复样本、错误格式和错误数据引用；相对路径可正确解析 |
| 方法与交付 | 用小型已知输入核对实际结果、计数、文件格式和样本身份 |
| 失败传播 | 工具非零退出时整个步骤失败，不生成新的COMPLETED清单 |
| 来源及恢复 | 记录输入指纹、软件/数据库版本和参数；相同命令重跑不重复计算；损坏结果不能直接复用 |
| 修改范围 | 只改某一步参数时，未受影响的前序步骤保持有效 |
| 安装 | 构建wheel，在仓库外安装并运行，确认规则、配置和模板齐全 |

正式模块按自身数据增加专用结果校验函数，并传给verify_cached_outputs；通用校验器只检查清单中声明的文件，不能理解某个指标的生物学含义。示例的一行序列统计不能替代测序数据质控或完整生物学分析。

## 6. 代码说明和发布检查

每个生产函数写docstring：小函数说明用途和返回值；复杂函数补充输入、输出、副作用、失败条件和关键约束。注释解释分析选择、坐标口径、来源核验和资源分配，不逐行复述语句。参数含义写入--help；YAML示例逐key注明类型、单位、可选值和覆盖关系。

```bash
python -m pylint --rcfile=.pylintrc --enable=missing-function-docstring src env run.py
python -m pylint --rcfile=.pylintrc tests/test_packaging.py
python -m unittest discover -s tests -t . -p 'test_*.py'
python -m build
python -m pip install --force-reinstall --no-deps dist/omics_workflows-0.1.0-py3-none-any.whl
```

上述wheel安装命令适用于依赖已准备好的分析环境；`--force-reinstall`用于更新已有同版本框架，`--no-deps`保留已部署依赖。新环境先按环境定义准备所需软件。

版本由包内 `__version__` 定义；更新版本时同时检查安装示例。wheel包含代码和运行资源，源码包另含环境定义、文档和测试。发布前从仓库外验证安装版本的命令、资源、预检、预演及小型运行，不依赖开发安装掩盖漏打包文件。

## 7. 维护现有结果

组装继续使用layout_version=3及01至07编号目录；新增模块独立维护自己的布局版本。更改结果格式时明确兼容规则，不自动移动或覆盖旧运行；必要时使用新run-dir。[历史组装迁移工具](legacy/assembly_migrate_layout.py)仅处理旧01至10布局。

每步来源记录只纳入相关工具、代码和参数。完善/纠错前的坐标不能直接用于修改后的序列；程序完成、质量评估完成与生物学验收分别记录。真实样本、任务日志和验证证据保存在项目目录。
