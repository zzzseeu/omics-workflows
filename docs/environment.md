# 环境管理

框架可调用不同分析流程，每个流程按自身的 Python、生信工具和数据库需求选择环境。不同流程可以分别使用独立环境；依赖兼容时也可以复用同一环境。环境名称和安装位置由使用项目确定。

框架的通用 Python 依赖定义在 [pyproject.toml](../pyproject.toml) 中，模块的专用依赖和数据库要求见各自手册。`env/` 只保留 [environment.yaml](../env/environment.yaml)，作为可选依赖参考；使用前核对其内容是否匹配所选流程。实际分析不读取该文件，也不创建、安装或更新环境。

## 安装到所选环境

先准备满足目标流程要求的环境，再安装框架。以下 `analysis_env` 是示例环境名，可替换为项目自己的名称：

```bash
conda run -n analysis_env python -m pip install /path/to/omics-workflows
conda run -n analysis_env omics-workflows list
```

模块还有专用 Python 依赖时，按模块手册安装对应的 optional-dependencies。全部 Python 依赖已具备时，可加 `--no-deps`；更新已有同版本框架时可加 `--force-reinstall`。安装框架不部署外部生信工具或数据库。

开发环境可在仓库根目录执行 `python -m pip install -e '.[dev]'`，开发工具版本由 pyproject.toml 的 dev 选项管理。

## 选择运行环境

名称由 conda/micromamba 解析，框架使用所选环境的 Python 和分析工具。以下 `workflow_name` 表示已注册的流程名称，替换后查看其参数：

```bash
conda run -n analysis_env omics-workflows workflow_name --help
micromamba run -n analysis_env omics-workflows workflow_name --help
```

也可先激活环境，再运行 `omics-workflows workflow_name ...`。各流程使用不同环境时，分别在其环境中安装和调用框架即可。

流程提供 `--env-prefix` 时，可通过环境路径选择已有环境；具体参数见该流程的 --help。当前实现会切换分析环境的 Python 和工具，框架脚本与规则仍来自本次调用的安装版本。它不会创建或更新目标环境。

Snakemake 本身支持按规则声明 Conda 环境或容器；当前框架的运行器尚未启用这些部署方式。现有运行方式是选择已准备好的环境后执行。Snakemake 容器支持见[官方说明](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)。

## 准备与记录环境

需要新建环境时，按目标流程的手册选择软件、版本及数据库。若仓库中的参考定义适合所选流程，可在仓库根目录使用：

```bash
micromamba create -n analysis_env -f env/environment.yaml
micromamba run -n analysis_env python -m pip install .
```

使用项目保存实际解析的软件包清单、源码工具来源、数据库版本和验证记录，数据与部署产物均位于仓库外。软件环境建好后，按模块提供的预检和小型案例确认其可运行；环境检查通过与实际分析结果的质量验收分别记录。

HPC 上由项目自己的 CSUB 作业申请资源，内部使用 Snakemake local executor。环境部署和实质验证放在计算任务中。升级环境前确认当前任务不再依赖该版本；已有任务续跑宜沿用原环境及框架版本。

## 新流程开发

通用 Python 依赖写入 pyproject.toml 的 dependencies，模块专用 Python 依赖写入对应 optional-dependencies。流程专用生信工具、版本、数据库和特殊部署要求写入该流程手册；有依赖冲突时单独准备环境，无需安装其他流程的工具。

各模块在运行前检查自身依赖，并记录实际使用的软件来源。框架公共模块负责可复用的环境检查与运行能力，新流程按[开发与维护](development.md)接入和验证。
