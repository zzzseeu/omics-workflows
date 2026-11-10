# omics-workflows

omics-workflows 提供可复用的组学分析流程和统一调用接口。实际样本配置、启动脚本、CSUB 投递脚本、数据、日志和结果由各项目保存。

## 目录

```text
omics-workflows/
├── run          # 统一入口
├── assembly/    # 组装模块：启动器、规则、分析脚本、配置与示例
├── common/      # 通用配置检查、来源记录与任务运行
├── env/         # 统一环境定义、包锁和源码锁
├── tests/
└── docs/        # 公共说明与各流程手册
```

所有流程共用一个 `omics_workflows` micromamba 环境。普通服务器可直接调用接口；在 HPC 上，由项目提交一个完整分析作业，内部使用 Snakemake local executor 执行各步骤。

## 使用入口

```bash
./run --help
./run assembly --help
```

当前支持 assembly。完整配置示例位于 `assembly/examples/`；模块的 defaults.yaml 定义默认值，schema.yaml 定义字段约束。分析模式和交付范围见流程手册。

## 文档入口

从 [文档导航](docs/README.md) 开始，按需要阅读：

- [快速开始](docs/quickstart.md)：项目如何准备配置、调用和恢复。
- [Assembly 手册](docs/workflows/assembly.md)：模式、输入、关键参数、输出和审阅。
- [环境管理](docs/environment.md)：统一环境的使用、安装和维护。
- [开发与维护](docs/development.md)：新增流程、文档模板、验证和旧结果迁移。

新增流程时，在 `docs/workflows/` 下添加手册并更新导航。公共文档集中说明使用约定，项目记录实时运行状态，历史设计另行归档。
