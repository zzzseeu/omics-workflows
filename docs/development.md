# 开发与维护

仓库保存通用代码、规则、参数约束和测试；项目保存真实样本、原始数据、启动/CSUB 脚本、日志和结果。内部使用 Snakemake local executor，项目决定是否外层提交任务。

## 代码结构

- [run](../run) 选择流程；[common/workflow.py](../common/workflow.py) 提供通用配置检查、软件检查、来源记录、运行锁、信号处理和 Snakemake 启动。
- [assembly/runner.py](../assembly/runner.py) 提供组装专属的数据要求、倍性与策略约束、软件清单和步骤定义，调用公共模块执行检查和启动。公共模块不导入 assembly。
- [规则](../assembly/rules/)与输出编号一致。[scripts/](../assembly/scripts/) 只保留 `01_` 至 `07_` 七个步骤入口；[lib/](../assembly/lib/) 保存合并后的共享实现，不再提供逐工具命令行。
- 每步参数清单只记录相关工具和源码。改变纠错或评估参数不应重新执行未受影响的初始组装。
- 输入完整性、程序完成和生物学验收分开；命令失败不得转换为跳过或成功。
- 向外部工具传参使用列表，避免 shell 插值。新产物先在独立目录生成、核验，再发布；不替换不受流程管理的用户目录。

三种组装策略及各步固定输入输出只在 [Assembly 手册](workflows/assembly.md)说明，不另维护设计文档或进阶手册。

`common/`、`assembly/`、`assembly/scripts/`、`assembly/lib/` 和测试目录以 `__init__.py` 明确声明为 Python 包。规则、示例和文档目录不作为包导入。

## 格式检查与测试

开发依赖见 [development.txt](../env/development.txt)，pylint 配置见 [.pylintrc](../.pylintrc)。复杂度和公共 API 文档约定单独审阅，导入、命名、编码、格式和错误检查保持开启。

```bash
python -m pip install -r env/development.txt
python -m pylint --rcfile=.pylintrc assembly common env run
python -m unittest discover -s tests -t . -p 'test_*.py'
```

统一测试命令从仓库根运行，`-t .` 指定项目根目录，避免把 `tests/assembly/` 识别为顶层 assembly 包。

按变更范围检查输入破损、样本错配、双端顺序、单套与分相计数、相位标签、关闭/缺数据状态、原生结果解析、输出清理与来源校验。小型工具检查与真实组装结论分别记录；本站实际计算和合成组装在 CSUB 计算任务中进行。

发布前运行 pylint、相关测试、三种策略的 DAG 预演及 Markdown 链接检查。缺依赖的 dry-run 仅证明步骤依赖可构建，不表示实际工具已可运行。

## 布局与历史结果

新运行 layout_version=3，输出固定为 01_inputs、02_read_qc、03_assembly、04_finishing、05_polishing、06_assessment、07_report。旧运行和旧迁移记录保持只读，不自动移动或覆盖。继续计算使用新的 run-dir；已有交付可按完整清单核验后作为输入。

[历史目录迁移工具](legacy/assembly_migrate_layout.py) 仅用于历史 01…10 布局的旧版迁移，不会将旧结果转换为当前七步布局，也不会把历史验证当作新代码的验证。

## 新增流程

提供流程 runner、Snakefile、defaults.yaml、schema.yaml、编号规则、分析脚本、examples 和测试。runner 向公共模块提供字段约束、软件清单、步骤路径和结果核验函数；公共模块接收这些定义，不写入物种、倍性或具体工具选择。在统一入口注册流程，新增一份手册并更新导航。
