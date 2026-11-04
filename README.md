# omics-workflows

提供可复用的组学分析接口与流程实现。实际样本配置、项目启动脚本、CSUB 投递脚本、数据、日志和结果归各项目管理。

## 目录

```text
omics-workflows/
├── run          # 统一入口
├── assembly/    # 组装模块：规则、分析脚本、配置与示例
├── common/      # 运行管理
├── env/         # 统一环境定义、包锁和源码锁
├── tests/
└── docs/        # 公共说明与各流程手册
```

所有流程共用一个 `omics_workflows` micromamba 环境。普通服务器直接调用接口；HPC 在项目中为完整任务外层提交一次作业，内部使用 Snakemake local executor。

## 使用入口

```bash
./run --help
./run assembly --help
```

当前已接入 assembly。完整配置示例位于 `assembly/examples/`，默认值及字段约束位于模块的 defaults.yaml 和 schema.yaml。能力与交付边界见流程手册。

## 文档入口

从 [文档导航](docs/README.md) 开始，按需要阅读：

- [快速开始](docs/quickstart.md)：项目如何准备配置、调用和恢复。
- [Assembly 手册](docs/workflows/assembly.md)：模式、输入、关键参数、输出和审阅。
- [环境管理](docs/environment.md)：统一环境的使用、安装和维护。
- [开发与维护](docs/development.md)：新增流程、文档模板、验证和旧结果迁移。

公共说明集中维护；新增流程在 `docs/workflows/` 下增加手册并登记导航。项目实时运行状态和历史设计不混入当前使用说明。
