# omics-workflows

omics-workflows 是可安装的生信分析流程框架，通过同一命令调用各分析模块。当前提供基因组组装模块；后续流程在包内注册，复用配置检查、来源记录和运行管理。样本数据、项目配置、CSUB 投递脚本、日志和结果保存在各项目中。

## 安装与调用

支持 Linux，要求 Python 3.12 或更高版本；当前验收使用 Python 3.12。先按[环境管理](docs/environment.md)准备分析软件，再把框架安装到该环境：

```bash
# 在目标环境中，从源码安装框架及组装模块的 Python 依赖
python -m pip install '/path/to/omics-workflows[assembly]'

omics-workflows --version
omics-workflows list
omics-workflows assembly --help
```

安装框架不会部署 hifiasm、HapHiC、GPhase 或数据库。若环境已具备全部 Python 依赖，可用 `python -m pip install --no-deps /path/to/omics-workflows` 安装框架。普通服务器直接运行；HPC 上由项目提交一个完整作业，内部使用 Snakemake local executor。

```bash
omics-workflows assembly \
  --config /path/to/project/config/assembly.yaml \
  --run-dir /path/to/project/runs/run001 \
  --cores 8 --mem-mb 64000 --check
```

安装后可在任意目录调用，也可使用 `python -m omics_workflows`。保留源码调用方式 `/path/to/omics-workflows/run.py assembly ...`，参数相同。`--env-prefix` 切换到已有分析环境，不创建或更新环境；框架自带脚本和规则仍来自本次调用的安装版本。

代码、规则和模板按包内相对位置定位，不写死项目或安装路径。配置中的表格路径相对 YAML，读段路径相对 readsets.tsv；运行时解析绝对路径用于执行与来源记录。

## 目录

```text
omics-workflows/
├── pyproject.toml
├── src/omics_workflows/
│   ├── cli.py                # 命令入口与流程注册
│   ├── common/               # 通用运行能力和资源读取
│   └── workflows/assembly/   # 组装规则、脚本、参数和示例
├── run.py                       # 源码调用兼容入口
├── env/                      # 生信环境定义、包锁和源码锁
├── tests/
└── docs/
```

## 文档

- [快速开始](docs/quickstart.md)：安装后如何准备项目、调用和恢复。
- [Assembly 手册](docs/workflows/assembly.md)：输入输出、步骤、工具和文件格式。
- [环境管理](docs/environment.md)：分析软件环境的部署与验证。
- [开发与维护](docs/development.md)：模块接入、构建、pylint 和验证。
