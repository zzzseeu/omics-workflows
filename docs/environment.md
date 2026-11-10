# 环境管理

所有流程共用一个 `omics_workflows` micromamba 环境。项目指定环境前缀；流程规则不安装软件。环境、安装日志和验证产物均位于仓库外。

## 定义与锁文件

| 文件 | 用途 |
|---|---|
| [environment.yaml](../env/environment.yaml) | 统一软件需求 |
| [linux-64.lock](../env/linux-64.lock) | 历史基础环境包锁；不包含七步流程新增依赖 |
| [versions.tsv](../env/versions.tsv) | 关键软件版本速查 |
| [haphic.lock.json](../env/haphic.lock.json) | HapHiC 官方源码固定版本 |

environment.yaml 包含当前软件需求；linux-64.lock 和 versions.tsv 保存历史部署基线，不能当作本次新增工具已安装或已验证的证据。创建新环境后保存实际解析包锁；源码工具还需要对应源码锁。

tidk 0.3.1 属于最终端粒检查的基础依赖。T2T 完善的可选依赖另见 [t2t.lock.json](../env/t2t.lock.json)，不属于基础 conda 包锁。

## 使用已有环境

激活环境后调用统一入口，或传入 `--env-prefix /path/to/envs/omics_workflows`。该选项用于切换到已有环境，不创建或更新环境。调用案例见 [快速开始](quickstart.md)。

## 创建与验证

由项目自己的启动脚本或计算任务调用环境管理工具。[create.sh](../env/create.sh) 拒绝覆盖现有前缀，默认按完整 environment.yaml 创建环境，再执行 [HapHiC 源码安装](../env/install_haphic.py) 和 [离线自检](../env/validate.py)：

```bash
export MICROMAMBA_EXE=/path/to/micromamba
export OMICS_ENV_PREFIX=/path/to/envs/omics_workflows
export ASSEMBLY_DEPLOYMENT_DIR=/path/to/project/results/environment_install
bash /path/to/omics-workflows/env/create.sh
```

MICROMAMBA_EXE 未设置时从 PATH 查找。保留 ASSEMBLY_ENV_PREFIX 别名；未指定部署目录时必须提供仓库外的 OMICS_RESOURCE_ROOT。可用 ASSEMBLY_ENV_SPEC 显式覆盖安装规范，ASSEMBLY_VALIDATE=0 跳过自检；跳过后不能报告验证通过。

已有环境的单独自检：

```bash
/path/to/envs/omics_workflows/bin/python /path/to/omics-workflows/env/validate.py \
  --output /path/to/project/results/environment_validation
```

本站的软件安装、自检和实质分析均放在计算任务中。当前 HapHiC 未部署可选 MKL，使用上游支持的 dense-matrix 路径；高度碎片化或大型组装应依据实际资源记录评估预算。

## 新流程增加依赖

先检查已有工具，只补充缺少的依赖。更新主定义、精确包锁、版本清单及所需源码锁，并对既有流程做回归验证，避免为新增模块直接降级已验证的基础依赖。所有流程继续共用一个环境，等当前分析不再使用该环境时再升级。

环境安装成功、功能验收通过、真实样本质量合格是不同结论。历史安装和验收来源见 [维护档案](../../project/omics_workflows_maintenance/docs/README.md)。

## GPhase 多倍体扩展

GPhase 与已运行的二倍体流程共用 `omics_workflows` 的同一个 Python，不另建环境。为保留运行中任务的现有依赖，新增 Python 包放在该前缀的 `opt/GPhase/python/`，原始源码放在 `opt/GPhase/source/`，新增原生工具放在 `opt/GPhase/bin/`。多倍体入口使用这些附加搜索路径，既有 site-packages、conda history 和 linux-64 基础锁保持不变。

- [gphase.lock.json](../env/gphase.lock.json)：上游 commit、源码逐文件指纹、附加包版本、mosdepth/YaHS 来源与二进制指纹。
- [gphase.requirements.txt](../env/gphase.requirements.txt)：附加 Python 包的固定版本。
- [gphase_budget.patch](../env/gphase_budget.patch)：将外层线程预算传给 GPhase 的 haplotype 聚类，避免上游默认 8×8 并行。

当前锁定 GPhase `a9221035274e7397fc05dc5be693270c2758528c`、YaHS 1.2.2（`9b2cc15740f03e84dffe248fab04445c994f17ab`）、mosdepth 0.3.11。statsmodels 使用 0.15.0 以兼容现有 pandas；不降低公用 pandas/NumPy 版本。popCNV GUI 和 Pore-C/PPL 不用于本路线。

新部署须在计算任务中完成以下操作：按锁文件获取 GPhase 源码并应用预算补丁，安装附加包，按固定 commit 编译 YaHS 并校验 mosdepth 二进制，建立指向该前缀内部的 gphase/mosdepth/yahs 入口，生成与安装位置对应的 `opt/GPhase/SOURCE_LOCK.json`。入口会检查源码指纹。

安装后还需执行完整作者案例及本仓库四拷贝、拷贝数、比对、交付和绘图测试。能显示 help 只能确认入口可启动。

screed 1.1.3 的旧构建元数据需要兼容的 setuptools。可在同一前缀的临时构建目录放置 setuptools 80.9.0、setuptools-scm 8.3.1、wheel 0.45.1，再以 `--no-build-isolation --no-deps --target .../opt/GPhase/python` 安装固定依赖；构建目录不加入分析运行时搜索路径。不要整体升级正在使用的基础环境。部署来源以现有锁文件和项目中的实际验证记录为准。


## T2T 完善工具

04 至 06 使用的工具共用同一个 `omics_workflows` 前缀，附加目录为 `opt/T2T/`；基础 conda 包和现有 conda history 保持不变。固定来源和逐文件指纹见 [t2t.lock.json](../env/t2t.lock.json)，安装器为 [install_t2t.py](../env/install_t2t.py)。

| 工具 | 固定来源 | 用途 |
|---|---|---|
| quarTeT 1.3.1 | a44d8c87e3f4bcefaf43dc4134e88b6a21890265 | 经证据筛查的局部 GapFiller fill；不使用参考基因组强制排序 |
| tidk 0.3.1 | Bioconda hef77e2b_0，包 SHA256 固定 | 指定 motif 的端粒分布 |
| NextPolish2 0.2.2 | Bioconda h74ec884_0，包 SHA256 固定 | HiFi＋普通 WGS k-mer 纠错，显式 reference phase 模式 |
| Inspector tag v1.3 | 2d813d96eb38adec9b8fb25388a3f4f1a2c63573 | 最终序列的读段错误定位；上游 help 版本字段较旧，以 commit 为准 |
| hifiasm/minimap2/samtools/yak/fastp | 既有统一环境 | 局部组装、比对、WGS QC 和 k-mer |
| Merqury/meryl、BUSCO | 既有统一环境 | QV、k-mer 完整性及逐 copy 的保守基因完整性 |

Inspector 的附加 statsmodels 0.15.0/patsy 1.0.2 安装在 `opt/T2T/python/`，仅 Inspector 包装器加入自己的模块搜索路径。原生二进制包装器只为该进程补充同前缀的共享库搜索路径。官方 tidk 0.3.2 Linux 发布包要求较新的 glibc，因此本站选择经过计算节点验证的 0.3.1 Bioconda 构建。

首次部署先在可联网主机缓存，再在计算任务中安装：

```bash
"$OMICS_ENV_PREFIX/bin/python" env/install_t2t.py --cache /external/t2t_sources --download-only
# 以下安装命令放在项目的计算任务中；不调用 conda solve。
"$OMICS_ENV_PREFIX/bin/python" env/install_t2t.py --cache /external/t2t_sources
```

安装器拒绝覆盖已锁定的扩展，保留基础环境中版本一致的 tidk。新位置需重新安装，不能直接复制包含绝对包装器路径的前缀。安装不等于功能验收；新布局的原生组件检查入口为 `tests/assembly/native_seven_stage_smoke.py`，输出须位于仓库外。已有真实作业继续用冻结源码，不自动启用新模块。

BUSCO 需要项目明确提供适合物种的离线 lineage 目录，预检要求其中包含 dataset.cfg。规则使用 `--offline`，不会隐式下载数据库。NextPolish2 的官方许可限学术和非商业使用。未提供普通 WGS 时，流程会明确记录纠错和 Merqury 的缺失状态，不能用 Hi-C 代替普通 WGS。
