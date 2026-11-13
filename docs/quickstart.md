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

## 2. 安装与查看入口

各流程按自身依赖选择环境，可以独立准备，也可复用兼容的已有环境。先按流程手册和[环境管理](environment.md)准备软件，再使用所选环境的 Python 安装框架；以下以当前组装模块为例：

```bash
python -m pip install '/path/to/omics-workflows[assembly]'
omics-workflows --version
omics-workflows list
omics-workflows assembly --help
```

环境已具备全部依赖时，可加 `--no-deps` 安装框架，避免重新解析依赖。离线部署使用预先构建的 wheel。`--help` 不要求生信软件已部署；实际运行时检查所选流程的依赖。

安装后命令可从任意目录执行；`python -m omics_workflows.cli` 等效。在仓库根目录也可直接运行该模块；开发时可用 `python -m pip install -e '.[assembly]'` 安装。按环境名称调用可用 `conda run -n omics_workflows omics-workflows assembly ...`，micromamba 同样支持 `run -n`；`--env-prefix` 按路径选择已有环境，不负责安装或更新环境。

## 3. 准备项目配置

按所选流程的手册，将模板复制到项目并填写真实输入。Assembly 的完整模板见 [模块示例目录](../omics_workflows/assembly/examples/)，参数定义见 [Assembly 手册](workflows/assembly.md)。

先做只读元数据预检，再查看 DAG：

```bash
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/project

omics-workflows assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/assembly.yaml" \
  --run-dir "$PROJECT_ROOT/runs/run001" \
  --cores 8 --mem-mb 64000 --check

omics-workflows assembly \
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
python -m omics_workflows.assembly.examples.prepare_demo \
  --project-dir "$DEMO_PROJECT"

omics-workflows assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$DEMO_PROJECT/config/config.yaml" \
  --run-dir "$DEMO_PROJECT/runs/demo001" \
  --cores 4 --mem-mb 4096 --target report
```

生成器要求项目目录尚不存在，以免覆盖已有数据。恢复分析时，只运行第二条命令。可用 `--samples sample_A` 选择样本，省略时运行样本表中的全部样本。

真实基因组分析须按流程手册准备数据、填写倍性并估算资源，不能沿用生成器写入的小型测试参数。
