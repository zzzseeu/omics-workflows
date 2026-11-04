# 快速开始

本页说明所有流程共用的调用方式。具体输入、参数和产物见 [流程手册导航](README.md)。当前可用模块是 `assembly`。

## 1. 划分仓库与项目

omics-workflows 提供统一入口、通用规则和分析实现。项目保存真实样本、数据清单、调用参数、执行脚本和结果：

```text
project/
├── config/           # 实际配置和样本表
├── scripts/          # 调用统一入口的项目启动脚本
├── jobs/             # 申请资源的项目投递脚本
├── logs/
└── runs/             # 分析结果、来源记录和 .snakemake 状态
```

原始数据可以引用其他路径，无需搬进项目。运行目录必须在代码仓库外，原始输入也必须在运行目录外。

## 2. 查看入口

```bash
/path/to/omics-workflows/run --help
/path/to/omics-workflows/run assembly --help
```

`--help` 不需要数据。实际运行使用统一的 `omics_workflows` 环境：激活后直接调用，或使用 `--env-prefix` 指定前缀。环境准备见 [环境管理](environment.md)。

## 3. 准备项目配置

根据所选流程手册，将对应模板复制到项目并填写真实输入。Assembly 的完整模板位于 [模块示例目录](../assembly/examples/README.md)；参数定义见 [Assembly 手册](workflows/assembly.md)。

先做只读元数据预检，再查看 DAG：

```bash
WORKFLOW_ROOT=/absolute/path/to/omics-workflows
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/project

"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/assembly.yaml" \
  --run-dir "$PROJECT_ROOT/runs/run001" \
  --cores 8 --mem-mb 64000 --check

"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/assembly.yaml" \
  --run-dir "$PROJECT_ROOT/runs/run001" \
  --cores 8 --mem-mb 64000 --dry-run
```

`--check` 不创建运行目录；`--dry-run` 写入解析计划但不执行分析。二者都不能替代实际运行时的 MD5 和数据全量解析。

## 4. 执行和恢复

普通服务器去掉检查开关运行同一命令。HPC 在项目 `jobs/` 中申请单节点资源，外层一次提交完整任务，再调用同一接口。规则内部使用 Snakemake local executor。本站的实质计算和合成数据生成在计算节点执行。

`--cores`、`--mem-mb` 应与实际可用资源一致。它们约束内部调度预算，不能代替外层资源申请；每步估算也不能保证实际峰值不超限。

恢复时沿用原配置和 run-dir；比较不同方案时使用新的 run-dir。保留 `.snakemake` 和来源记录，同一运行目录只允许一个启动器。具体流程的恢复边界见对应手册。

## 5. 用小型案例熟悉接口

以下命令使用合成数据，仅验证调用方法。在计算任务中首次生成一个仓库外的新项目：

```bash
DEMO_PROJECT=/absolute/path/to/new_demo_project
"$OMICS_ENV_PREFIX/bin/python" "$WORKFLOW_ROOT/assembly/examples/prepare_demo.py" \
  --project-dir "$DEMO_PROJECT"

"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$DEMO_PROJECT/config/config.yaml" \
  --run-dir "$DEMO_PROJECT/runs/demo001" \
  --cores 4 --mem-mb 4096 --target report
```

生成器要求项目目录不存在，不覆盖已有数据。恢复案例只运行第二条命令。可用 `--samples sample_A` 选择样本；省略时运行样本表全部样本。

生成器写入的小型测试参数不能用于真实基因组。真实项目须按流程手册重新准备数据、倍性和资源预算。
