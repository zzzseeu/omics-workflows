# 第一次运行 Assembly：单样本 HiFi 组装

本例从一份 PacBio HiFi FASTQ 或 FASTQ.gz 开始，完成一个样本的 contig 组装和基础统计。需要流程代码、可用的 `omics_workflows` 环境，以及自己的 HiFi 数据。染色体级组装的数据选择见[主手册](assembly.md#1-选择组装路线)。

## 1. 确认软件与数据位置

在 Bash 中填写四个绝对路径。`OMICS_ENV_PREFIX` 是环境目录，其中应有 `bin/python`；`PROJECT_ROOT` 使用仓库外的新项目目录，读段保留在原位置。

```bash
WORKFLOW_ROOT=/absolute/path/to/omics-workflows
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/new_assembly_project
HIFI_FASTQ=/absolute/path/to/sample.hifi.fastq.gz

test -x "$WORKFLOW_ROOT/run"
test -x "$OMICS_ENV_PREFIX/bin/python"
test -s "$HIFI_FASTQ"
```

三项检查均应返回成功。已有 micromamba 时可用 `micromamba env list` 查找环境路径；环境不存在时，先按[环境创建与验证](../environment.md#创建与验证)完成安装和自检。`--env-prefix` 只使用已有环境，不安装软件。

## 2. 创建配置和两张表

本例只使用一个样本 `sample01`。以下命令会创建配置文件，仅在新的项目目录执行一次；恢复组装时不要重复执行。

```bash
mkdir -p "$PROJECT_ROOT/config" "$PROJECT_ROOT/jobs" "$PROJECT_ROOT/logs"

cat > "$PROJECT_ROOT/config/config.yaml" <<'YAML'
schema_version: 1
samples_table: samples.tsv
readsets_table: readsets.tsv
defaults:
  assembly:
    mode: hifi_only
    representation: primary_alternate
    threads: 8
    mem_mb: 64000
  scaffolding:
    enabled: false
samples:
  sample01:
    ploidy: 2
YAML

printf 'sample_id\tspecies\tprofile\nsample01\tExample_species\t\n' \
  > "$PROJECT_ROOT/config/samples.tsv"
printf 'sample_id\treadset_id\tlibrary_id\ttechnology\tformat\tfile1\nsample01\thifi01\tlib01\thifi\tfastq\t%s\n' \
  "$HIFI_FASTQ" > "$PROJECT_ROOT/config/readsets.tsv"
```

运行前检查三处：`samples.tsv` 的物种名改为实际物种；YAML 的 `ploidy` 按已知材料填写 1 或 2，未知填 `null`；`readsets.tsv` 的末列须为自己的 HiFi 文件。四倍体不要直接修改本例的倍性，应使用[四倍体组装配置](assembly-advanced.md#1-四倍体组装)。

本例使用 8 核、64,000 MB，只作调用示范。真实基因组按大小、杂合度和测序量调整线程、内存，以及下面实际运行命令和外层作业申请的资源。

## 3. 检查配置和待运行步骤

```bash
"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/config.yaml" \
  --run-dir "$PROJECT_ROOT/runs/run001" \
  --cores 8 --mem-mb 64000 --target report --check
```

成功时终端输出 `PRECHECK PASS`，并列出 `samples=sample01`、`modes=hifi_only`。若报错，先修正路径、文件、配置或环境，不继续执行。将末尾 `--check` 改为 `--dry-run`，可检查将运行的分析步骤；这两种检查都没有实际组装读段。

## 4. 执行组装

普通服务器去掉上面命令末尾的 `--check` 即可运行。本站 HPC 须在项目 `jobs/assembly.csub.sh` 保存以下脚本；修改三行路径变量和两个日志路径，保持它们与前面的项目一致。

```bash
#!/usr/bin/env bash
#CSUB -J assembly_sample01
#CSUB -n 8
#CSUB -R span[hosts=1]
#CSUB -M 64G
#CSUB -o /absolute/path/to/new_assembly_project/logs/assembly.%J.out
#CSUB -e /absolute/path/to/new_assembly_project/logs/assembly.%J.err
set -Eeuo pipefail

WORKFLOW_ROOT=/absolute/path/to/omics-workflows
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/new_assembly_project

exec "$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/config.yaml" \
  --run-dir "$PROJECT_ROOT/runs/run001" \
  --cores 8 --mem-mb 64000 --target report
```

在终端提交，保存返回的作业号：

```bash
csub < "$PROJECT_ROOT/jobs/assembly.csub.sh"
```

整个组装只提交一次，不再单独提交各步骤。队列有额外要求时，在作业脚本中补充相应设置。不要同时在普通终端和计算任务中对同一运行目录启动组装。

## 5. 查看进度和结果

HPC 的终端输出保存在 `logs/assembly.<作业号>.out`、`.err`；具体步骤日志在 `runs/run001/samples/sample01/logs/`。作业开始执行后，可查看：

```bash
cat "$PROJECT_ROOT/runs/run001/run_status.json"
```

`COMPLETED` 且退出码为 0 表示所选分析完成。排队期间运行目录可能还不存在；失败时先查看作业日志和相应步骤日志。完成后应查看以下文件：

| 文件（相对于 `runs/run001/`） | 用途 |
|---|---|
| `samples/sample01/04_assembly/hifiasm/primary.fa` | 主组装 contig 序列 |
| `samples/sample01/04_assembly/hifiasm/alternate.fa` | 替代 contig，可为空 |
| `samples/sample01/09_assessment/basic.tsv` | 序列条数、总长度、N50、最长 contig、N 碱基数等 |
| `samples/sample01/10_report/report.html` | 可离线打开的样本报告 |
| `10_reports/batch.html` | 本次分析的样本汇总 |

本例没有使用 Hi-C，结果不是两套完整分相组装，也不会生成染色体挂载和接触图。完成计算后仍须结合预期基因组大小等判断组装质量；基础统计不等于 BUSCO、QV 或 T2T 验收。

## 6. 中断后恢复

先确认原作业及子进程全部退出，再重新执行原组装命令，或重新提交同一份作业脚本。保留原配置、数据、运行目录、`.snakemake` 和日志，不重新创建配置或覆盖样本表。原方案继续使用 `run001`；更改组装方案时另建运行目录。
