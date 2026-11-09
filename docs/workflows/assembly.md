# Assembly：基因组组装

使用 PacBio HiFi 组装 contig，结合 ONT 超长读段（UL）和 Hi-C 获得分相的染色体级组装，并可继续补洞、端粒延伸、碱基纠错和质量评估。多样本分别组装。

[输入准备](#2-准备输入) · [分析步骤](#3-步骤与输出) · [运行示例](#4-运行示例) · [结果解读](#5-查看结果与评估质量) · [进阶操作](assembly-advanced.md)

第一次使用可先按[单样本 HiFi 组装教程](assembly-first-run.md)完成配置、检查、提交和结果查看，再选择其他组装路线。

## 1. 选择组装路线

| 数据与材料 | `assembly.mode` | `assembly.representation` | 组装方法与主要结果 |
|---|---|---|---|
| HiFi；单倍体、二倍体或未知倍性 | `hifi_only` | `primary_alternate` | hifiasm；primary/alternate contig |
| HiFi + Hi-C；二倍体 | `hifi_hic` | `phased` | hifiasm 分相、HapHiC 挂载；两套染色体候选序列 |
| HiFi + ONT UL + Hi-C；二倍体 | `hifi_ul_hic` | `phased` | 增加 UL 辅助组装；两套染色体候选序列 |
| HiFi + ONT UL + Hi-C；四倍体 | `hifi_ul_hic` | `unitigs` | hifiasm 原始 unitig 图、GPhase 分相挂载；四份染色体候选序列 |

染色体级组装须填写倍性和每套染色体数。本项目的四倍体路线须指定 `scaffolding.backend: gphase`，详见[四倍体组装](assembly-advanced.md#1-四倍体组装)。本流程尚未接入其他已知多倍体路线，二倍体与四倍体应分别运行；这里说明的是本项目可用配置，不代表 hifiasm 软件全部能力。

基因组特征估计（03）尚未接入。HiFi-only 只提供 contig 及基础统计，补洞、纠错和扩展评估仅用于染色体级组装。

## 2. 准备输入

准备 `omics_workflows` 环境，见[环境管理](../environment.md)，再复制[配置模板](../../assembly/examples/README.md)，填写以下三份文件。

| 文件 | 内容 |
|---|---|
| `config.yaml` | 样本表和读段表的位置、组装路线、倍性、染色体数、线程与内存 |
| `samples.tsv` | `sample_id`、`species`、`profile`；每行一个样本 |
| `readsets.tsv` | 每行一组测序数据，记录所属样本、物理文库、测序技术和文件位置 |

### 测序数据要求

| `technology` | 格式 | 用途 |
|---|---|---|
| `hifi` | 未配对 BAM、FASTQ 或 FASTQ.gz | 主组装；补洞和后续纠错时提供长读段证据 |
| `ont` | 已完成 basecalling 的 UL FASTQ 或 FASTQ.gz | 跨越长重复，辅助组装、补洞和端粒延伸 |
| `hic` | 成对 FASTQ 或 FASTQ.gz | 辅助分相、染色体挂载和接触图检查 |
| `wgs` | 普通基因组双端短读段 FASTQ 或 FASTQ.gz，可选 | NextPolish2 碱基纠错、Merqury QV 和 k-mer 完整性评估 |

Hi-C 不能代替普通 WGS。启用 BUSCO 时，另需适合物种的离线 lineage 目录，其中包含 `dataset.cfg`。

### 样本与读段表

HiFi-only 示例的 `samples.tsv` 如下。`profile` 对应 YAML 的 `analysis_profiles`，也可留空。

```tsv
sample_id	species	profile
sample_A	Example species A	hifi_contigs
sample_B	Example species B	hifi_contigs
```

`readsets.tsv` 至少包含以下六列。同一样本可以有多组 HiFi/ONT 数据，`readset_id` 在该样本内须唯一。

```tsv
sample_id	readset_id	library_id	technology	format	file1
sample_A	cell01	library_A01	hifi	bam	../data/sample_A.cell01.bam
sample_A	cell02	library_A02	hifi	fastq	../data/sample_A.cell02.fastq.gz
sample_B	cell01	library_B01	hifi	fastq	../data/sample_B.cell01.fastq.gz
```

Hi-C/WGS 增加 `file2` 列填写 R2，HiFi/ONT 留空。压缩 FASTQ 的 `format` 仍填 `fastq`。可选列 `md5_1`、`md5_2` 填交付 MD5，`prepared_manifest` 用于[复用已准备读段](assembly-advanced.md#2-复用已准备的-hifiont-数据)。

同一物理 Hi-C 文库的多 lane 须先合并，再登记一行，保证对完整文库去重；不同文库分行，路径不能重复。ID 使用英文字母、数字、下划线、点或短横线，首字符为字母或数字。

TSV 路径相对于 YAML，读段路径相对于 `readsets.tsv`，也可用绝对路径。原始数据放在运行目录外，运行目录放在代码仓库外。

## 3. 步骤与输出

下列路径均相对于运行目录的 `samples/<sample_id>/`。各路线只执行所需步骤，03 的编号保留。

| 步骤 | 目的与处理方法 | 主要输出 |
|---|---|---|
| 01 输入检查与准备 | 核对文件校验值、Hi-C 配对关系；将 HiFi/ONT 整理为 FASTQ，或复用已验证的准备结果 | `01_inputs/`：长读段 `*.fastq.gz`、输入校验与来源记录 |
| 02 长读段质控 | 用 seqkit 统计读段数、碱基数和读长，检查 FASTQ；复用数据沿用已通过的统计 | `02_read_qc/`：各组读段的 QC 结果与汇总 |
| 04 contig/unitig 组装 | hifiasm 组装；二倍体联合路线用 Hi-C 分相，UL 路线加入超长读段；四倍体保留原始 unitig 图 | `04_assembly/hifiasm/`：FASTA、GFA、FASTA 索引 |
| 05 Hi-C 比对与连接证据 | BWA 比对，按物理文库去重和配对过滤；四倍体另用 HiFi 比对和 popCNV 估计 unitig 拷贝数 | `05_linkage/`：Hi-C BAM、比对统计；四倍体另有拷贝数结果 |
| 06 分相与染色体挂载 | HapHiC（二倍体）或 GPhase（四倍体）聚类、排序和定向，核对 AGP 与 FASTA 并绘制 Hi-C 接触图 | `06_scaffolding/delivery/`：染色体候选、未挂载序列和 AGP；`contact_maps/`：接触图 |
| 07 补洞与端粒延伸，可选 | 用局部长读段组装和 quarTeT 补洞；仅接受有唯一锚定和跨越读段支持的修改 | `07_finishing/`：`genome.fa`、`edits.tsv` 和支持证据 |
| 08 碱基纠错，可选 | 用 HiFi 和普通 WGS 运行 NextPolish2，检查纠错前后各拷贝的一致性；关闭或缺 WGS 时记录跳过原因 | `08_polishing/`：`genome.fa`、各份 `hap*.fa` 和纠错状态 |
| 09 组装质量评估 | 基础统计包括长度、N50 和缺口等；后处理可选 tidk 端粒检查、Merqury、Inspector 和逐拷贝 BUSCO | `09_assessment/basic.tsv`；扩展评估为 `t2t/assessment.json`、`t2t/chromosomes.tsv` |
| 10 报告 | 汇总序列、统计、来源和待检查问题 | `10_report/report.html`、`report.md`、`manifest.json`；后处理报告位于 `10_report/t2t/` |

07 至 09 的后处理须显式启用，见[补洞、纠错与扩展评估](assembly-advanced.md#3-补洞纠错与扩展评估)。全部参数见[默认值](../../assembly/defaults.yaml)和[参数定义](../../assembly/schema.yaml)。

## 4. 运行示例

### HiFi-only：从项目配置开始

将以下路径改为实际位置，复制模板后编辑 YAML 和两张 TSV。模板含两个示例样本，样本表、读段表和 YAML 的 `samples` 设置须一起调整。

```bash
WORKFLOW_ROOT=/absolute/path/to/omics-workflows
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/project

mkdir -p "$PROJECT_ROOT/config/hifi_only"
cp "$WORKFLOW_ROOT/assembly/examples/hifi_only/"*.{yaml,tsv} \
  "$PROJECT_ROOT/config/hifi_only/"

"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/hifi_only/config.yaml" \
  --run-dir "$PROJECT_ROOT/runs/hifi_only_001" \
  --cores 8 --mem-mb 64000 --target report --check
```

`--check` 核对样本、路径、参数和软件，不创建运行目录。改为 `--dry-run` 可查看分析步骤，去掉检查选项即执行分析；完整文件校验和 FASTQ 解析在实际分析时完成。

### HiFi + UL + Hi-C：二倍体染色体级组装

复制另一套模板，填写 `ploidy: 2`、`chromosomes_per_haplotype` 和 Hi-C 酶切位点，按上一节方法检查后执行实际组装命令：

```bash
mkdir -p "$PROJECT_ROOT/config/hifi_ul_hic"
cp "$WORKFLOW_ROOT/assembly/examples/hifi_ul_hic/"*.{yaml,tsv} \
  "$PROJECT_ROOT/config/hifi_ul_hic/"

# 复制后先编辑 YAML 和两张 TSV，完成检查后再执行以下分析。
"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/hifi_ul_hic/config.yaml" \
  --run-dir "$PROJECT_ROOT/runs/chromosomes_001" \
  --cores 32 --mem-mb 240000 --target chromosomes
```

没有 UL 时改用 `assembly.mode: hifi_hic`，保留 HiFi 和 Hi-C。酶切位点默认 `GATC`，须按建库协议填写，并在 `scaffolding.enzyme_provenance` 记录依据。

| `--target` | 得到什么 |
|---|---|
| `contigs` | contig 或四倍体 unitig 图及索引，止于组装 |
| `chromosomes` | 挂载结果、接触图及样本报告；所有选中样本须启用挂载，不继续补洞纠错 |
| `report`，默认 | 所选路线的结果和样本报告，以及运行目录下的 `10_reports/` 批次汇总；包含已启用的后处理 |
| `t2t` | 补洞、纠错和扩展评估及独立报告；所有选中样本须启用 `finishing.enabled`，可接续已完成的染色体组装 |

省略 `--samples` 时处理样本表全部样本，用 `--samples sample_A` 可只运行一个样本。配置按“模块默认值 → YAML `defaults` → 样本 `profile` → YAML `samples`”依次覆盖。

资源按基因组大小、杂合度和测序量调整，每步内存不能超过申请总量。染色体路线至少需 4 核、32,000 MB，实际可能更多。hifiasm 默认 `bloom_filter_bits: 37`，不能沿用小型合成数据的 `0`。

普通服务器直接运行；本站以 `csub < 项目投递脚本.sh` 申请单节点资源，一次提交完整分析，由 Snakemake 按线程和内存安排各步骤。实质计算和合成数据生成均在计算节点执行。试运行见[小型合成案例](../quickstart.md#5-用小型案例熟悉接口)。

## 5. 查看结果与评估质量

先看 `samples/<sample_id>/10_report/report.html` 和 `09_assessment/basic.tsv`，再按上表找序列。常用结果包括：

- HiFi-only：`04_assembly/hifiasm/primary.fa`、`alternate.fa`；alternate 可以为空，两者不等于两套完整单倍型。
- 染色体级：`06_scaffolding/delivery/chromosomes.fa`、分开的 `hap*.chromosomes.fa`、`unplaced.fa`、AGP 和 `SHA256SUMS`。查看 `06_scaffolding/contact_maps/` 下的 `contact_map.pdf`、`separate_plots.pdf`，检查连接、方向和分相混合。
- 补洞纠错后：使用 `08_polishing/genome.fa` 与同版本 `hap*.fa`；逐染色体评估见 `09_assessment/t2t/chromosomes.tsv`，报告见 `10_report/t2t/report.html`。07 修改序列后，旧 AGP 不再表示新坐标。

`COMPLETED` 只表示选定分析完成。未审阅的染色体结果和 `CANDIDATE_NEEDS_REVIEW` 均需继续验收；缺失的质量评估不能视为通过。高 N50、染色体数吻合或无 N，均不能单独证明组装正确或达到 T2T，还需结合 Hi-C、独立相位、污染和复杂重复区域检查，见[人工审阅](assembly-advanced.md#4-hi-c-挂载与人工审阅)。

hap/copy 标签不代表父母本或亚基因组，染色体编号按组内长度排序；标准染色体命名须另做共线性分析。

## 6. 恢复与排错

恢复时仅重新执行原 `run assembly` 命令，保留配置、样本选择、运行目录、`.snakemake` 和日志，无需再复制模板。同一目录只运行一个分析。新方案另建目录；更换数据、参数或软件版本可能重跑相关步骤，进行中的组装继续使用原代码和环境。

| 问题 | 先检查什么 |
|---|---|
| 缺输入、未知样本或参数 | TSV 列名、样本 ID、相对路径和 YAML；缺少必需测序数据会报错 |
| MD5 不符或 FASTQ 无法解析 | 交付校验值、文件完整性、测序文件格式 |
| 内存不足或 Bloom 设置不合理 | 外层作业申请、每步内存、基因组规模及已有资源统计 |
| Hi-C 读对保留很少或挂载数量不符 | 比对和去重统计、接触图、原始组装图及未挂载序列；相似单倍型可参考[分别比对](assembly-advanced.md#4-hi-c-挂载与人工审阅) |
| 缺普通 WGS | 纠错按 `polishing.on_missing_wgs` 跳过或报错，Merqury 记录缺失；不可改用 Hi-C |
| 中断后提示运行目录被占用 | 确认旧任务及子进程全部退出，再按实际情况解除 Snakemake 锁并恢复 |

日志和资源统计分别在样本目录的 `logs/`、`benchmarks/`。核对 `run_status.json` 与本次样本报告，避免误读旧结果。旧输出布局迁移见[开发与维护](../development.md#5-assembly-旧结果布局迁移)。
