# 开发与维护

## 1. 代码与项目的边界

仓库提供统一入口、流程规则、通用分析代码、参数约束和测试。项目管理真实样本配置、启动/CSUB 脚本、数据、日志和结果。所有流程使用统一环境，内部用 Snakemake local executor 执行，项目决定是否在外层提交作业。

当前源码布局见 [仓库 README](../README.md)。公共入口为 [run](../run)，运行管理位于 [common/runner.py](../common/runner.py)，Assembly 的 Snakefile、配置解析、默认值、schema、规则和分析脚本集中在 [assembly/](../assembly/)。

## 2. 新增一个流程

| 位置 | 需要提供或修改的内容 |
|---|---|
| `<流程名>/` | Snakefile、configuration.py、defaults.yaml、schema.yaml、rules/、scripts/；配置模板统一放 examples/ |
| `run` | 登记模块及其入口 |
| `common/runner.py` | 首次接入第二个模块时参数化；把组装专用配置、工具、目标、计划写出和代码指纹逻辑交给模块 |
| `env/` | 按需增加依赖并更新锁文件；保持统一环境 |
| `tests/<流程名>/`、接口测试 | 配置约束、最小真实工具闭环、失败/恢复及既有流程回归 |
| `docs/workflows/<流程名>.md` | 按下一节模板写流程手册，并在导航中登记 |

当前入口仅注册 assembly，公共运行器仍含 assembly 专用逻辑，尚未实现自动模块发现。新增模块前，应按上表将接口参数化。源代码指纹应按模块和步骤限定失效范围，并通过测试确认新增模块不会触发既有分析重算。

模块内部遵循以下约定：

- `examples/<案例>/` 保存一套完整匹配的配置和表格；模块顶层不再重复放置 `*.example.*`。`defaults.yaml` 是运行时默认值，`schema.yaml` 是参数约束，两者不是示例。
- `rules/NN_步骤.smk` 与输出阶段编号对应，一个业务阶段一个文件，可包含多个规则。`Snakefile` 只负责目标、依赖辅助函数和阶段接入；未实现阶段不创建空规则文件。
- 简单工具调用写在规则的 `shell` 中；需要用 Python 或其他语言解析数据、转换格式、画图、校验多个产物或处理事务时，将代码放在 `scripts/`。脚本按“动作 + 对象”命名，例如 `validate_hic_pairs.py`、`plot_hic_contacts.py`，每个脚本有一个明确入口，避免将多项职责藏在综合任务分发器中。
- 脚本可复用 `*_utils.py` 中的辅助函数。规则须将入口脚本及其本地辅助依赖声明为输入；每阶段参数文件只携带相关代码和工具指纹，完整来源记录仍保存在运行计划中。
- 工具适配器若还负责多进程管道、失败现场保留、校验或完整目录发布，可以保留在脚本中，但须写明这些职责。不要仅为调用一条工具命令增加 Python 包装。

Assembly 各阶段的分析目的和产物见 [Assembly 手册](workflows/assembly.md#3-步骤与输出)，具体实现见 [规则](../assembly/rules/)和[分析脚本](../assembly/scripts/)。真实项目的启动脚本和投递脚本仍放在项目中。

各流程自行定义分析单位、样本选择方式、输入模型和目标。矩阵或群体分析不必沿用组装的 readsets 表。环境切换、资源预算、运行锁、状态、日志和来源记录等稳定功能由各流程共用。

## 3. 文档框架与新增手册模板

公共文档固定为导航、快速开始、环境和开发维护，只在公共约定变化时修改。每个已接入流程在 workflows/ 下维护使用手册，可按以下结构组织；复杂方法和进阶操作另页说明，由主手册链接：

```text
# 流程名与分析目的
## 1. 分析路线与适用范围
## 2. 输入准备
## 3. 分析步骤、目的与输出
## 4. 配置与运行示例
## 5. 结果解读与质量评估
## 6. 恢复与排错
```

手册说明输入要求、已实现范围、关键选择及结果解释。完整默认值和 schema 链接到源码，完整可运行案例链接到模块 examples/。公共调用语法、环境创建和 HPC 规则只链接公共文档，不在每份手册重复。

能力表只登记已验证的路径。未来设计、未实现模块、合成数据验证和真实生物学验收需分别标明，真实项目的状态放在项目文档中。历史设计及部署记录移入维护档案，当前导航只保留一个档案入口。新手册通过评审和最小闭环验证后，登记到 [文档导航](README.md)。

## 4. 验证与发布

新增或调整代码后，按影响范围依次验证：配置和接口测试 → 仓库外最小分析 → 失败或中断后的恢复 → 受影响既有流程回归。本站实质测试由维护项目的计算任务执行，不应为验证重构而重跑真实组装。

检查源码和文档路径、Markdown 链接、示例配置与 CLI 参数是否一致。保留验证代码、环境锁、来源指纹和外部验证记录。分析计算完成、质量评估完成和人工审阅完成应有独立状态。

四倍体扩展的历史基线通过 50 项单元＋4 项接口检查（77961575，DONE）。四倍体原始图/CNV/Hi-C/恢复、GPhase 作者完整案例、实际输出坐标及四份导出/重比对/接触图/报告分别通过 77961563、77961516、77961567（均 DONE）。[完整证据](../../omics-resources/deployments/polyploid_extension_20261103/validation.json) 分别记录真实三倍体作者案例与合成四拷贝验证；这些结果不代表水稻材料的生物学相位或 T2T 认证。阶段化源码的二倍体、失败恢复、审阅及依赖隔离历史见 [维护记录](../../project/omics_workflows_maintenance/docs/assembly_stage_refactor_20261103.md)，案例生成验收为 77958207（DONE）。

07 至 09 初版扩展通过 62 项单元＋4 项接口检查（77967316，DONE）。冻结源码后的验证包括：二倍体缺 WGS 闭环（77967317，DONE）、四拷贝带 WGS 闭环与已知 SNP 纠错（77967318，DONE）、离线 BUSCO 及质量失败传播（77967315，DONE）、既有 HiFi-only 多样本/校验失败/中断恢复回归（77967343，DONE）。两个完善闭环均检查无重算恢复和发布后的证据链接。[验证记录](../../omics-resources/deployments/t2t_modules_20261104/validation.json) 保存固定源码、工具锁、失败迭代与生物学边界。这些验证使用合成数据，只验证功能，不代表水稻 T2T 或独立相位认证。测试代码为 `tests/assembly/t2t_integration.py`、`t2t_busco_integration.py` 和 `test_postassembly.py`。

输入准备复用升级通过 73 项单元＋4 项接口检查（77967638，DONE）；原生 BAM/FASTQ、三种校验模式、并发单次校验、源文件时间戳保护、QC 复用、恢复及变更阻断通过 77967654（DONE）。原流程回归、启用复用的二倍体及四拷贝 T2T 功能闭环分别通过 77967655、77967656、77967657（均 DONE）。[复用验证记录](../../omics-resources/deployments/input_reuse_20261104/validation.json) 保存冻结源码和证据。真实项目约 789 GiB 的旧准备文件目前只核对了元数据，首次使用时才由计算任务流式校验；准备状态不表示真实数据校验已完成。

## 5. Assembly 旧结果布局迁移

Assembly 的配置 schema_version 为 1，输出 layout_version 为 2。调整源码目录不等于更改输出布局版本。新输出保持 01 至 10 编号，其他流程可自行定义稳定的阶段编号，无需凑满十步。

已完成旧结果使用 [migrate_layout.py](../assembly/migrate_layout.py)，默认仅打印计划：

```bash
"$OMICS_ENV_PREFIX/bin/python" "$WORKFLOW_ROOT/assembly/migrate_layout.py" \
  --run-dir /path/to/completed/run
```

确认后在项目计算任务中添加 --apply，执行目录移动及完整性核验。工具要求旧运行 COMPLETED/exit_code=0，使用同一运行锁，拒绝处理运行中任务和目标冲突。

迁移将 inputs、qc/reads、assembly、qc/assembly、report 分别映射至 01_inputs、02_read_qc、04_assembly、09_assessment、10_report，批次 reports 改为 10_reports；FASTA/FAI/GFA 前后校验。原命令、时间、代码指纹和 `.snakemake` 证据保留，journal 位于 layout_migration/。

迁移后的运行作为历史结果保存，当前入口拒绝向其写入新的分析计划。查看结果直接读取更新后的 manifest；新分析使用新 run-dir。运行中的任务继续使用其原有冻结源码和环境。
