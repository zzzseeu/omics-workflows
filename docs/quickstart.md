# 快速开始

所有流程使用同一调用接口，当前可用模块是 `assembly`。具体输入、参数和产物见 [流程手册导航](README.md)。

## 1. 划分仓库与项目

omics-workflows 保存统一入口、通用规则和分析代码；项目保存真实样本、数据清单、调用参数、执行脚本和结果：

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

按所选流程的手册，将模板复制到项目并填写真实输入。Assembly 的完整模板见 [模块示例目录](../assembly/examples/)，参数定义见 [Assembly 手册](workflows/assembly.md)。

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

`--check` 不创建运行目录。`--dry-run` 会写入解析后的计划，但不执行分析。实际运行时仍需校验 MD5 并全量解析数据。

## 4. 执行和恢复

在普通服务器上，去掉检查开关即可运行同一命令。在 HPC 上，由项目 `jobs/` 中的脚本申请单节点资源，一次提交完整任务，再调用同一接口。各规则由 Snakemake local executor 执行。本站的实质计算和合成数据生成均在计算节点执行。

`--cores`、`--mem-mb` 应与实际可用资源一致，内部调度按这两个预算分配资源。在 HPC 上仍需由外层任务申请资源；各步骤声明的用量是估算值，实际峰值可能超限。

恢复时沿用原配置和 run-dir，比较不同方案时使用新的 run-dir。保留 `.snakemake` 和来源记录，同一运行目录只允许一个启动器。各流程的恢复条件见对应手册。

## 5. 用小型案例熟悉接口

以下命令用合成数据验证调用方法。首次使用时，在计算任务中生成一个仓库外的新项目：

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

生成器要求项目目录尚不存在，以免覆盖已有数据。恢复分析时，只运行第二条命令。可用 `--samples sample_A` 选择样本，省略时运行样本表中的全部样本。

真实基因组分析须按流程手册准备数据、填写倍性并估算资源，不能沿用生成器写入的小型测试参数。
