# Assembly：基因组组装

从已知材料背景的 HiFi、可选 ONT 超长读段和 Hi-C 出发，完成单套参考或分相染色体组装，按需完善和纠错，再运行最终质量检查。样本真实倍性、每份染色体数、物种端粒重复序列和 Hi-C 建库信息由项目提供，不在组装流程中重新推断。

## 1. 分析路线与输入

| `assembly.strategy` | 组装与挂载方式 | 输出目标 |
|---|---|---|
| `single` | hifiasm primary → HapHiC；Hi-C 用于挂载，不执行 Hi-C 分相 | 单套参考 `reference` |
| `diploid` | hifiasm Hi-C 分相 → HapHiC | `hap1`、`hap2` |
| `polyploid` | hifiasm 原始 unitig 图 → HiFi/minimap2/popCNV 拷贝数估计 → GPhase 原生分相与挂载 | `hap1` 至 `hapP`，P 为实际倍性且至少为 3 |

`single` 表示输出一套参考，不改变材料倍性；二倍体水稻仍填 `ploidy: 2`。primary 可能包含不同同源染色体的片段，不宣称连续亲本相位。alternate 仅作为替代等位序列保留，不当作第二份完整单倍型。GPhase 内部使用 YaHS/HapHiC 模块，不在其最终结果后重复运行完整 HapHiC。

| 数据 | 文件格式 | 要求及用途 |
|---|---|---|
| HiFi | 未配对 BAM、FASTQ、FASTQ.gz | 所有路线必需，主组装及长读段支持 |
| ONT UL | 已完成 basecalling 的 FASTQ/FASTQ.gz | `hifi_ul_hic` 必需；多倍体采用此组合 |
| Hi-C | R1/R2 FASTQ 或 FASTQ.gz | 染色体组装必需；登记完整物理文库和酶切信息 |
| 普通 WGS | R1/R2 FASTQ 或 FASTQ.gz | 可选；NextPolish2、QV 和 k-mer 完整率。Hi-C 不能替代 |
| BUSCO lineage | 离线数据库目录，含 `dataset.cfg` | 完整最终评估必需，选择适合物种的数据库 |

单套及二倍体可选 `assembly.mode: hifi_hic` 或 `hifi_ul_hic`；多倍体用 `hifi_ul_hic`。仅有 HiFi 时保留 `hifi_only` contig 终点，不将其报告为染色体组装或完整七类最终评估。

统一环境的使用和安装见[环境管理](../environment.md)。新增依赖已列入 [environment.yaml](../../env/environment.yaml)；当前旧环境可能缺少 NanoPlot、NanoFilt、HiFiAdapterFilt 或 LTR 工具，先运行 `--check`。`--dry-run` 允许展示缺依赖的 DAG，但会明确列出缺项，不能据此认为具备实际运行条件。

## 2. 配置文件

准备以下三份文件。原始文件和运行目录均放在代码仓库外；原始文件也不能放在运行目录内。上游原生脚本不保证处理空白字符路径，运行目录使用无空白字符的路径。

| 文件 | 内容与格式 |
|---|---|
| `config.yaml` | YAML：表格位置、策略、已知倍性、染色体数、数据库、工具参数和资源 |
| `samples.tsv` | 制表符分隔：`sample_id`、`species`、`profile`；一行一个样本，profile 可空 |
| `readsets.tsv` | 制表符分隔：`sample_id`、`readset_id`、`library_id`、`technology`、`format`、`file1`、`file2` |

例如 `samples.tsv`：

```tsv
sample_id	species	profile
rice01	Oryza sativa	""
```

`readsets.tsv` 中 HiFi/ONT 的 `file2` 留空。压缩 FASTQ 的 format 仍填 `fastq`：

```tsv
sample_id	readset_id	library_id	technology	format	file1	file2
rice01	hifi01	H1	hifi	bam	../data/hifi.bam	""
rice01	hic01	C1	hic	fastq	../data/hic.R1.fastq.gz	../data/hic.R2.fastq.gz
rice01	wgs01	W1	wgs	fastq	../data/wgs.R1.fastq.gz	../data/wgs.R2.fastq.gz
```

同一样本可登记多个 HiFi/ONT 数据组。readset_id 在样本内唯一；同一物理 Hi-C 文库的多 lane 须先合并，避免分开标记重复。不同数据组不能重复引用同一原始文件。可选 `md5_1`、`md5_2` 为交付校验值。表格路径相对 YAML，读段路径相对 readsets.tsv，也可用绝对路径。

单套水稻配置示例：

```yaml
schema_version: 1
samples_table: samples.tsv
readsets_table: readsets.tsv
defaults:
  ploidy: 2
  chromosomes_per_haplotype: 12
  assembly:
    strategy: single
    mode: hifi_hic
    threads: 8
    mem_mb: 64000
    bloom_filter_bits: 0
  scaffolding:
    enabled: true
    threads: 8
    mem_mb: 64000
    processes: 4
    restriction_sites: GATC
    enzyme_provenance: replace_with_actual_library_protocol
  finishing:
    enabled: false
    gap_closing: false
    telomere_completion: false
    telomere_motif: TTTAGGG
  polishing:
    enabled: false
    on_missing_wgs: skip
  assessment:
    enabled: true
    busco_lineage: /path/to/appropriate_lineage
    threads: 8
    mem_mb: 64000
```

二倍体分相改为 `strategy: diploid`，倍性保持 2。多倍体改为 `strategy: polyploid`、`mode: hifi_ul_hic` 并登记 ONT UL，ploidy 填实际 P；流程按 P 决定输出份数，不能通过复制二倍体结果得到多倍体。

仅在需要修改组装时开启 finishing，再选择 gap_closing、telomere_completion；两项独立。碱基纠错由 polishing.enabled 独立控制，不要求先开启补洞。最终七类评估独立于这两个开关，完整染色体报告不能关闭它。

完整参数及类型见 [defaults.yaml](../../assembly/defaults.yaml) 和 [schema.yaml](../../assembly/schema.yaml)。配置模板见 [单套参考](../../assembly/examples/single_hic/config.yaml)、[二倍体分相](../../assembly/examples/hifi_ul_hic/config.yaml)和[多倍体分相](../../assembly/examples/polyploid/config.yaml)；替换样本数据、建库信息和数据库路径后再检查。

## 3. 七个步骤的实现与输入输出

`assembly/scripts/` 只保留七个带编号的步骤入口，共享实现位于 `assembly/lib/`。下列路径相对 `samples/<sample_id>/`，入口、目录编号与步骤一一对应。JSON 清单记录文件校验值、来源和参数，供后续步骤核对。

| 步骤与实现入口 | 输入 | 实现方式与工具 | 核心输出及格式 |
|---|---|---|---|
| **1 数据检查与准备** `01_prepare_inputs.py` | 样本、文库信息及原始 BAM/FASTQ | 校验文件和 MD5/SHA-256；SAMtools 转 HiFi BAM；pysam 全量解析及双端 ID/顺序检查；多倍体另查同类型长读段 ID 唯一性 | `01_inputs/`：准备后 FASTQ.gz、`manifest.json`、`reads.tsv` |
| **2 测序数据质控** `02_quality_reads.py` | 步骤 1 的读段与清单 | WGS/Hi-C 用 fastp；HiFi 用 NanoPlot 报告＋HiFiAdapterFilt；ONT 用 NanoPlot 报告＋NanoFilt；保留双端配对 | `02_read_qc/`：质控后 FASTQ.gz、HTML 质量报告及必要报告资源、`reads.tsv`、`manifest.json` |
| **3 基因组组装与染色体挂载** `03_assemble_genome.py` | 质控后读段、策略、真实倍性、每份染色体数及 Hi-C 建库信息 | 按第 1 节三选一。单套和二倍体用 HapHiC，多倍体用 GPhase；Hi-C 按物理文库标记重复；验证 FASTA/AGP 与输出份数 | `03_assembly/contigs/`：核心 FASTA/GFA；`delivery/`：染色体、各序列集合、未挂载序列、组成序列、AGP、命名表、清单及接触图 |
| **4 组装完善，可选** `04_finish_chromosomes.py` | 染色体交付、清洁 HiFi/ONT 和端粒 motif | minimap2/SAMtools 提取局部读段，hifiasm 局部组装；候选验证后用 quarTeT GapFiller 补洞，或接受有支持的末端延伸，tidk 检查 | `04_finishing/`：`genome.fa`、各集合 FASTA、`manifest.json`、`edits.tsv`；有尝试时保留统计及接受位点的必要证据 |
| **5 碱基纠错，可选** `05_polish_chromosomes.py` | 上一步组装、清洁 HiFi 与普通 WGS | yak 构建 21/31-mer 计数库，minimap2/SAMtools 回贴 HiFi，NextPolish2 纠错，再核对序列名称、方向、拷贝及前后对应 | `05_polishing/`：`genome.fa`、各集合 FASTA、`manifest.json`；纠错时保留前后 PAF |
| **6 最终质量评估** `06_assess_assembly.py` | 最终组装交付、质控后读段清单、BUSCO lineage、端粒 motif | 固定执行下表七类检查；针对最终版本重做回贴和 Hi-C 图，不修改序列 | `06_assessment/`：`assessment.json`、`metrics.tsv`、`chromosomes.tsv`、`telomeres.tsv`、必要原生统计及接触图 |
| **7 报告与交付** `07_write_assembly_report.py` | 最终组装及同版本评估结果 | 核对校验值，汇总指标、缺失检查和待审阅事项，不重复分析 | `07_report/`：样本 `report.html`、交付 `manifest.json`；运行目录 `07_report/summary.tsv`、`summary.html` 为批次汇总 |

关闭完善或纠错时，在对应编号目录链接原序列，记录 DISABLED，不再物理复制同一份基因组。缺少 WGS 时，纠错按 on_missing_wgs 跳过或报错，Merqury 标记缺数据；不会使用 Hi-C 代替。

### 最终七类质量检查

| 必查内容 | 固定工具 | 主要指标与交付 |
|---|---|---|
| 连续性与缺口 | QUAST | 总长度、数量、contig/scaffold N50、L50、GC、N；注明统计长度和分割口径 |
| 基因完整性 | BUSCO genome 模式 | 各集合 C/S/D/F/M、数据库版本、原生 summary 和 full_table |
| 碱基准确度与序列完整性 | meryl＋Merqury | QV、k-mer 完整率、拷贝数频谱；分别保留染色体、全部交付及各单倍型的统计口径 |
| 读段支持 | minimap2＋SAMtools＋mosdepth | HiFi 比对率、split-read 比例、覆盖广度、深度和窗口统计 |
| 植物重复区质量 | LTRharvest＋LTR_FINDER_parallel＋[LTR_retriever/LAI](https://github.com/oushujun/LTR_retriever) | 全局和窗口 LAI；完整及总 LTR 比例低于原论文适用范围时标记不适用 |
| 染色体挂载与结构 | 交付统计＋BWA/samblaster/SAMtools/HapHiC | 染色体数和长度、非 N 碱基挂载比例、未挂载长度、最终版本全局及逐染色体 Hi-C 图 |
| 端粒 | tidk ends | 每条染色体两端状态、阵列坐标/长度、距端点距离；按集合汇总检出末端数和双端检出数 |

单套路线对 reference 评估，分相路线逐份评估并保留合并结果。单份单倍型相对全样本 k-mer 的完整率不能当作整个基因组完整率。BUSCO 重复不直接等于冗余。HiFi 对全部染色体及未挂载竞争序列比对，覆盖统计明确 MAPQ 与序列范围。

tidk 的 present、absent、wrong_strand、not_terminal 状态直接保留；未检出不等于生物学缺失。无 N、两端端粒、较高 QV 或 BUSCO 均不能单独证明 T2T。流程不以自定 QV/BUSCO 阈值生成统一合格结论。

Inspector 错误定位、独立相位验证和污染筛查属于另行开展的专项，不自动加入七项默认路线。分相结果报告独立相位尚未验证，单套参考不适用相位误差。Hi-C 图保留人工审阅状态；T2T 目标还需复核着丝粒、rDNA 和所有接受的完善位点。

### 文件格式与精简原则

| 文件 | 格式与含义 |
|---|---|
| FASTQ.gz | 四行读段格式，经 gzip 压缩；双端保持 R1/R2 对应关系 |
| FASTA / `.fa.fai` | 染色体或 contig 序列及 SAMtools 索引 |
| GFA | hifiasm 组装图；多倍体 primary unitig 的 ID、顺序和序列保持原样 |
| AGP | 组成序列的顺序、方向和连接关系；与 components.fa 配套；补洞/纠错后旧 AGP 不直接代表新坐标 |
| `chromosome_names.tsv` | `sequence_id`、`source_object`、`haplotype`、`status`、`length`，单套集合登记为 reference |
| `edits.tsv` | 完善前后区段坐标和支持读段数；坐标从 0 开始、右端不包含 |
| PAF | 纠错前后序列比对，供拷贝和方向一致性核对 |
| `metrics.tsv` | `category`、`sequence_set`、`metric`、`value`、`unit`、`status` |
| `chromosomes.tsv` | 序列 ID、所属集合、长度、N/缺口、模糊碱基、左右端粒状态 |
| `telomeres.tsv` | tidk 原生逐末端 TSV，坐标从 0 开始、右端不包含；保留 strand 和 window_limited |
| JSON 清单 | 样本、策略、真实倍性、输出份数、文件校验、参数、工具来源及状态；不靠文件名猜测单倍型 |

步骤成功后默认清理组装缓存、临时 FASTA/FASTQ、回贴和 Hi-C BAM、BWA/GT 索引、yak/meryl 数据库及完整工具工作目录。保留用于解释指标的原生统计、必要图和接受修改的证据，不同时生成 Markdown/HTML/JSON 三份同内容报告。失败现场暂存于对应步骤父目录的 failed_attempts，排错后可清理。最终评估可用 assessment.keep_intermediates 或脚本 --keep-intermediates 保留排错现场。

步骤 3 的 chromosome FASTA、各集合 FASTA、unplaced.fa、components.fa、AGP 和核心 GFA 属于必要交付；单套输出 reference，分相输出 hap1…hapP。阶段日志在 `samples/<sample_id>/logs/01_inputs.log` 至 `07_report.log`，03 的组装和挂载分别记录。

## 4. 运行与恢复

在统一环境中运行，或用 --env-prefix 指定环境。HPC 上由项目提交一个 CSUB 作业，内部以 Snakemake local executor 执行；--cores、--mem-mb 与申请资源一致。

```bash
WORKFLOW_ROOT=/absolute/path/to/omics-workflows
OMICS_ENV_PREFIX=/absolute/path/to/envs/omics_workflows
PROJECT_ROOT=/absolute/path/to/project

"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/config.yaml" \
  --run-dir "$PROJECT_ROOT/runs/assembly001" \
  --cores 8 --mem-mb 64000 --target report --check
```

将 --check 改为 --dry-run 查看步骤依赖；去掉检查选项实际运行。--check 不创建运行目录，完整文件校验与 FASTQ 解析在实际步骤执行。--samples rice01 可选择样本。

| target | 终点 |
|---|---|
| `contigs` | 完成输入、质控和 contig/原始 unitig 图 |
| `chromosomes` | 完成染色体挂载交付 |
| `report` | 完整执行所选步骤、最终评估和批次报告 |
| `t2t` | 染色体路线完整后处理及报告的兼容目标；不表示认证 T2T |

恢复时保持原配置、环境和 run-dir，来源或参数变更由对应步骤重新处理。保留 config/、.snakemake/、清单和被链接的原始/旧交付文件。当前 layout_version=3；旧 01…10 布局保留为历史结果，不能在原目录直接接着计算，应使用新 run-dir。

### 独立脚本接口

三个后处理脚本均可独立调用 --help。固定输入为组装清单和质控后读段清单，参数由 Python 命令行指定：

```bash
python assembly/scripts/04_finish_chromosomes.py \
  --assembly-manifest /run/samples/rice01/03_assembly/delivery/manifest.json \
  --reads-manifest /run/samples/rice01/02_read_qc/manifest.json \
  --gap-closing --telomere-completion --telomere-motif TTTAGGG \
  --threads 8 --out-dir /run/samples/rice01/04_finishing

python assembly/scripts/05_polish_chromosomes.py \
  --assembly-manifest /run/samples/rice01/04_finishing/manifest.json \
  --reads-manifest /run/samples/rice01/02_read_qc/manifest.json \
  --kmer-sizes 21 31 --min-mapq 30 --threads 8 \
  --out-dir /run/samples/rice01/05_polishing

python assembly/scripts/06_assess_assembly.py \
  --assembly-manifest /run/samples/rice01/05_polishing/manifest.json \
  --reads-manifest /run/samples/rice01/02_read_qc/manifest.json \
  --busco-lineage /path/to/appropriate_lineage --telomere-motif TTTAGGG \
  --threads 8 --mem-mb 64000 --out-dir /run/samples/rice01/06_assessment
```

清单引用的数据必须保留且校验一致。完善参数包括锚定长度、一致性、支持读段数和局部窗口；纠错参数包括 k-mer、MAPQ、最小序列长度及缺 WGS 策略；评估参数包括 k-mer、MAPQ、覆盖窗口和 tidk 原生端粒参数。

## 5. 指标依据与结果边界

水稻组装中的连续性、BUSCO、QV、回贴和 LAI 组合可参见 [Y476 野生稻](https://doi.org/10.1038/s41467-024-48845-6)及[野生与栽培稻泛基因组](https://doi.org/10.1038/s41586-025-08883-6)。指标定义依据 [QUAST](https://doi.org/10.1093/bioinformatics/btt086)、[BUSCO](https://doi.org/10.1093/molbev/msab199)、[Merqury](https://doi.org/10.1186/s13059-020-02134-9)、[LAI](https://doi.org/10.1093/nar/gky730)和 [tidk](https://doi.org/10.1093/bioinformatics/btaf049)；多倍体路线依据 [GPhase](https://doi.org/10.1038/s41477-026-02422-5)。这是本流程选定的评估范围，不声称每篇论文执行同一清单。

评估执行完成、数值偏低、缺数据、不适用、工具失败和人工审阅分别记录。缺 WGS 的报告为 PARTIAL，LAI 不适用单列原因；工具失败返回非零并阻止生成新完成报告。高指标不自动转换成生物学验收结论。
