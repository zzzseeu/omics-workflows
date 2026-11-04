# Assembly：基因组组装

公共调用与项目边界见 [快速开始](../quickstart.md)，环境见 [环境管理](../environment.md)。本页集中维护组装的输入、配置、方法、产物和验收边界。

## 1. 能力与交付层级

| mode | representation | 必需数据 | 交付目标 |
|---|---|---|---|
| `hifi_only` | `primary_alternate` | HiFi | primary/alternate contigs、索引、统计报告 |
| `hifi_hic` | `phased` | HiFi、Hi-C | 两套分相组装、染色体候选 FASTA、AGP、接触图和报告 |
| `hifi_ul_hic` | `phased` | HiFi、ONT UL、Hi-C | 同上，并在组装中使用 UL 数据 |
| `hifi_ul_hic` | `unitigs`，backend=`gphase` | HiFi、ONT UL、Hi-C，明确 ploidy=4 | 原始 unitig 图、拷贝数估计、四份染色体候选、AGP、重新比对的接触图与报告 |

HiFi-only 可接受单/二倍体或未知倍性的初步组装；HapHiC 染色体路径必须显式提供 `ploidy: 2` 和每套染色体数。四倍体需要显式选择 GPhase/unitigs；其他已知多倍体或二/四倍体混合运行仍拒绝，不根据物种名推断核型。

03 基因组特征仍未接入。07 补洞/端粒完善、08 NextPolish2 和 09 扩展验收通过显式 `--target t2t` 接入，见本页“染色体组装后的完善与验收”。独立相位、污染、复杂重复拷贝和最终 Hi-C 人工审阅仍需外部证据；不会自动认证 T2T。其他 mode 或 representation 即使出现在 schema 中，也须经过运行预检确认，不能据此推断已开放。

## 2. 输入文件与样本表

`samples.tsv` 必须包含 `sample_id`、`species`、`profile` 三列。profile 留空表示使用默认及样本覆盖参数；非空时必须存在于 analysis_profiles。

`readsets.tsv` 每行代表一个 readset，列定义如下：

| 列 | 要求 |
|---|---|
| `sample_id` | 对应样本表 |
| `readset_id` | 同一样本内唯一 |
| `library_id` | 物理文库标识 |
| `technology` | `hifi`、`ont`、`hic`；`wgs` 在启用后处理纠错/k-mer 评估时参与 |
| `format` | 小写 `bam` 或 `fastq`；压缩 FASTQ 仍写 fastq |
| `file1` | 源文件路径 |
| `file2` | Hi-C/WGS 必填 R2；HiFi/ONT 留空 |
| `md5_1`、`md5_2` | 可选交付 MD5；提供时必须匹配 |
| `prepared_manifest` | 可选，仅 HiFi/ONT；指向已完成旧运行的 `samples/<样本>/01_inputs/<readset>.json`，启用准备结果复用 |

HiFi 支持未配对 BAM、普通/压缩 FASTQ；ONT UL 支持已完成 basecalling 的 FASTQ；Hi-C 为成对 FASTQ。FAST5 不属于该组装接口的必需输入。

同一样本多个 HiFi/ONT readset 可共同输入。一个物理 Hi-C 文库的多 lane 应先合并，再登记为一个 readset，使重复标记覆盖完整文库；不同文库分别处理。活动源路径不能重复登记。

ID 使用 ASCII 字母、数字、下划线、点或短横线，首字符为字母/数字。数据路径支持空格和中文。两张 TSV 的位置相对于 YAML，file1/file2 相对于 readsets.tsv，也可使用绝对路径。

缺少当前模式所需技术会报错，不自动降级。未提供 MD5 时仍计算源指纹，并记录 NOT_PROVIDED；这不代表已与交付方原件核对。

### 复用已准备的 HiFi/ONT 输入

在 readsets.tsv 保留原始 `file1`、`format` 和 MD5，另填 `prepared_manifest`。路径与 file1 一样相对于 TSV，或使用绝对路径。YAML 顶层必须指定 `prepared_input_cache: /project/results/assembly_inputs/prepared_reuse_cache`，缓存必须在新 run-dir 和流程仓库之外，不写入旧运行。

预检只读取小清单、QC 和文件签名，核对 sample/readset/library/technology、源路径及原始签名、交付 MD5 和无损转换策略。缺少准备 FASTQ、未完成 PASS QC、样本不符或来源不符时明确报错，不静默回退到重新转换。旧布局或其他软件的清单不能仅凭文件名复用。

首次使用旧版清单时，没有转换后 FASTQ 的历史 SHA256，因此在计算节点单次流式读取压缩文件：同时计算 SHA256、由 seqkit 解析 FASTQ 并核对原 QC 的 read/base 数量。该校验为现存准备文件建立来源缓存，不声称重新对 BAM 逐读段比对。新准备清单记录 FASTQ SHA256 和签名；若签名仍匹配，可直接使用生产者哈希，否则首次复用需要流式复核哈希。Snakemake 正常触碰输出时间戳也由该分支处理。

缓存按来源身份、旧清单/QC 哈希和准备文件路径寻址；文件锁避免并发消费者重复校验。命中后只核对小清单和 size/mtime/ctime/inode，不重复全文件读取。已缓存文件的签名变化会阻断，不能自动建立新的“正确”基线。删除缓存后可重新做只读校验，但仍不需要重新转换或复制。

新 `01_inputs/*.fastq.gz` 是指向旧准备文件的软链接，JSON 保留原始源哈希和准备 FASTQ 哈希，两者分开记录。新 02 复用原 PASS seqkit 统计并更新展示路径，不重新解析同一文件；T2T-only 目标本身也不会额外执行 02。没有 prepared_manifest 的普通输入继续使用原有准备/QC 路线。

必须保留被引用的旧 FASTQ、准备清单和 QC；原始源签名也必须可检查。移动或清理旧运行前应明确管理这些引用。软链接不会使旧数据变为可删除缓存。输入复用与版本变化后对新组装的 HiFi/ONT 比对是两件事，07–09 仍会运行针对其对应 FASTA 的必要比对。

## 3. 配置与模式选择

完整示例：[HiFi-only](../../assembly/examples/hifi_only/config.yaml)、[HiFi+UL+Hi-C](../../assembly/examples/hifi_ul_hic/config.yaml)。同目录包含匹配的样本表和 readset 表。将三份文件复制到项目后修改占位值；真实数据不沿用合成案例的小型资源配置。

配置顶层为 schema_version、samples_table、readsets_table，以及可选的 prepared_input_cache、defaults、analysis_profiles、samples。覆盖优先级为：模块默认值 → defaults → 样本所选 profile → samples 中该样本设置。

默认值以 [defaults.yaml](../../assembly/defaults.yaml) 为准，字段定义见 [schema.yaml](../../assembly/schema.yaml)。常用设置：

| 字段 | 作用 |
|---|---|
| `assembly.mode`、`representation` | 按上面的能力表选择 |
| `assembly.threads`、`mem_mb` | 单样本组装预算；不得超整份任务预算 |
| `assembly.bloom_filter_bits` | hifiasm 的 -f；默认 37，不能直接套用测试中的 0 |
| `assembly.purge_level` | null 使用后端默认；0—3 显式设置 -l |
| `ploidy`、`chromosomes_per_haplotype` | 倍性与每套目标染色体数，按材料填写 |
| `genome_size_bp` | 规模元数据，不据此强行截断或删减组装 |
| `scaffolding.enabled` | HiFi-only 为 false，两个 phased 路径为 true |
| `scaffolding.restriction_sites` | 按建库协议填写，默认 GATC |
| `scaffolding.enzyme_provenance` | 记录酶切位点的确认依据；未知默认不等于建库已确认 |
| `scaffolding.correct_nrounds` | 自动断点校正最大轮数，默认 0；按问题证据启用 |
| `scaffolding.mapq`、`nm_exclusive` | 成对过滤阈值，默认 MAPQ≥1、NM<3 |
| `scaffolding.hic_reference` | 默认 `combined`；二倍体已分相 contig 可显式选 `per_haplotype`，分别比对并分别运行 HapHiC，再汇总两套交付 |
| `scaffolding.min_retained_fraction` | 每个参考的 Hi-C 成对保留率下限，默认 0.001（0.1%），用于阻止信号几乎耗尽后继续聚类；不是通用的样本质量合格标准 |
| `scaffolding.gap_size`、`plot_bin_kb` | 连接 N 长度及接触图分箱 |
| `scaffolding.review_agp`、`review_manifest` | 成对提供已审阅修改，见第 7 节 |

未知参数会报错。染色体路径额外要求至少 4 核和 32,000 MB 总预算，每一步的声明预算也必须满足；这些下限不是实际基因组资源保证。Bloom bits=37 自身约占 17.2 GB，入口会检查明显不合理的内存组合。

## 4. 运行目标

公共选项为 `--config`、`--run-dir`、`--cores`、`--mem-mb`；可指定 `--env-prefix`、`--samples A B`、`--check`、`--dry-run`。完整参数通过 `run assembly --help` 查看。

| target | 当前执行范围 |
|---|---|
| `contigs` | 组装与索引，不生成评估/报告或开展挂载 |
| `chromosomes` | 仅适用于所有选中样本均启用挂载的任务；包括交付、接触图和样本报告 |
| `t2t` | 显式启用 07–09 完善/筛查和独立后处理报告；可接续外部已完成交付，不自动认证 T2T |
| `report`（默认） | 各样本对应路线的结果与报告，另加批次汇总 |

联合路线使用 hifiasm 的 --h1/--h2，UL 模式另加 --ul；此时不使用 --primary。输出 hap1/hap2 是 Hi-C 辅助分相标签，不代表父母本身份或已验证的全局相位。

以上 --h1/--h2 对应二倍体路线。GPhase 四倍体路线使用 hifiasm `--n-hap 4 -l 0 --ul` 生成原始 p_utg 图，不用二倍体的 Hi-C contig 分相输出；`contigs` 目标在该路线交付 unitig 图及索引。GPhase 及其附加依赖必须先在同一个公用 Python 前缀下准备好，见 [环境管理](../environment.md)。

### 四倍体配置与方法

```yaml
defaults:
  ploidy: 4
  chromosomes_per_haplotype: 12
  assembly:
    mode: hifi_ul_hic
    representation: unitigs
    purge_level: 0
    threads: 16
  scaffolding:
    enabled: true
    backend: gphase
    threads: 16
    processes: 4
    correct_nrounds: 0
```

这是方法片段，实际配置仍需完整的 schema_version、样本/数据表及与外层任务匹配的内存设置。ploidy=4、每份 12 条时，目标是四份共 48 条；1C 大小只是来源元数据，不据此截断每份组装。GPhase 的 copy 标签仅在染色体组内有效，不表示父母本，也不直接表示 C/D/K/L 亚基因组身份。

路线为：原始输入验证 → HiFi/UL 全量 QC 和跨 readset ID 唯一性检查 → 未经改名、过滤、purge 或 polishing 的 p_utg 图 → HiFi/minimap2/popCNV 拷贝数估计 → 按物理 Hi-C 文库 BWA/samblaster 比对并保留 MAPQ0 → GPhase 图辅助分相及挂载 → unitig 和 contig 两种表示的 AGP/FASTA 一致性校验 → 四份染色体候选 → Hi-C 重新比对到 GPhase 最终 contig components 后绘图 → 评估与报告。初次聚类使用 GPhase 默认 MAPQ≥1，scaffolding 使用 MAPQ≥0；接触图使用成对 MAPQ/NM 过滤。

GPhase 要求原始 p_utg 图，不能将二倍体 hap1/hap2 contig 简单扩成四份。其图闭合与 collapsed unitig 复制属于算法结果；导出会验证复制序列与原始 unitig 相同、AGP 坐标及同一副本区间不重叠，并核对原生 FASTA。图重叠允许被裁剪，不能套用 contig 路线的“每个原始碱基恰好一次”声明；未分配原始 unitig 保留在 unplaced.fa，可能包含图冗余。

每个 chr×copy 组必须恰有一条 scaffold；缺组、多条或未分类对象会标记 INCOMPLETE_SCAFFOLDING，保留产物和接触图并阻断最终成功，不挑最大 scaffold 凑 48 条。四份 FASTA 和索引位于同一 delivery/ 下的 hap1–hap4.chromosomes.fa；components.fa、unitig_components.fa、两类 AGP 及 contig_to_unitig.txt 保留图来源。

该路线目前要求无空白字符的 run-dir、gap_size=100、purge_level 为 0/null；使用 GPhase/YaHS 原生校正，不能设置 HapHiC 的 correct_nrounds。四倍体人工审阅结果回导尚未接入，提供 review_agp/review_manifest 会明确报错。后续可进入 07–09 完善和筛查；染色体组内四份标签仍不代表父母本或亚基因组身份，候选不自动成为 T2T。

算法依据：[GPhase 官方方法和输入要求](https://github.com/panlab-bioinfo/GPhase)。

## 5. 步骤与输出

源码按业务阶段组织，`Snakefile` 负责组装依赖和选择目标。同一步骤的多个规则放在同一个文件内：

| 规则文件 | 处理内容 | 调用的专用脚本（位于 scripts/） |
|---|---|---|
| `01_inputs.smk` | 长读长输入准备或身份校验后软链接复用、Hi-C 配对校验 | `prepare_long_reads.py`、`reuse_prepared_reads.py`、`validate_hic_pairs.py` |
| `02_read_qc.smk` | 运行 seqkit 或复用同一 FASTQ 的 PASS QC，解析统计结果 | `check_read_quality.py`、`summarize_read_qc.py` |
| `04_assembly.smk` | 两种组装路线、GFA/FASTA 规范化、直接建立 FASTA 索引 | `assemble_hifi_contigs.py`、`assemble_phased_contigs.py` |
| `05_linkage.smk` | 直接建立 BWA 索引、按物理文库比对过滤、保留配对顺序合并 | `map_filter_hic.py`、`merge_validate_hic_bams.py` |
| `06_scaffolding.smk` | HapHiC 挂载、AGP 校验与染色体导出、接触图 | `scaffold_haplotypes.py`、`export_chromosome_fasta.py`、`plot_hic_contacts.py` |
| `07_finishing.smk` | 唯一锚定读段、局部组装、补洞/端粒延伸及修改审计 | `finish_chromosomes.py` |
| `08_polishing.smk` | 普通 WGS 预处理、NextPolish2 及拷贝一致性筛查 | `polish_chromosomes.py` |
| `09_assessment.smk` | 基础统计；可选 tidk/Merqury/Inspector/BUSCO 及逐染色体检查 | `collect_contig_metrics.py`、`collect_chromosome_metrics.py`、`assess_t2t.py` |
| `10_report.smk` | 样本报告、染色体报告和批次汇总 | `write_sample_report.py`、`write_chromosome_report.py`、`write_batch_report.py` |

`scripts/` 保存需要自定义编程逻辑的处理。组装脚本除调用 hifiasm 外，还负责失败处理、ID/GFA 转换、来源绑定和完整产物发布；Hi-C 脚本负责检查每段管道退出状态、相邻 mate 和配对数量；HapHiC 适配器负责必需产物校验、目录发布和 Juicebox 脚本路径修复。接触图脚本根据交付 AGP 选择染色体对象，调用 HapHiC plot 并检查两类 PDF。

共享辅助代码仅供这些脚本导入，不作为独立分析步骤：

| 辅助文件 | 职责 |
|---|---|
| `io_utils.py` | JSON、文件指纹、子进程/管道、源文件变更检查、完整目录发布和失败保留 |
| `fasta_utils.py` | GFA→FASTA、FASTA 读写、序列统计和 ID 前缀处理 |
| `agp_utils.py` | AGP 坐标/方向/完整覆盖校验与序列重建 |
| `hic_bam_utils.py` | BAM 相邻配对和 mate 标记检查 |
| `report_utils.py` | Markdown 和转义后的离线 HTML 写出 |
| `prepared_input_utils.py` | 纯元数据的准备结果身份、来源、QC 和缓存描述校验 |

专用脚本的 `--help` 列出该任务的输入输出参数。用户运行分析仍通过仓库统一入口，通常无需单独调用脚本。未实现的 03 不创建空 `.smk` 文件。

输出布局为版本 2。样本产物位于 `samples/<sample_id>/`，阶段编号固定，跳过的阶段保留编号空位：

| 目录 | 当前用途 |
|---|---|
| `01_inputs/` | 全文件 MD5/SHA256、输入准备；Hi-C 全量配对格式、数量与 ID 检查 |
| `02_read_qc/` | HiFi/ONT 全量 FASTQ 解析和统计 |
| `03_genome_profile/` | 预留，未接入 |
| `04_assembly/hifiasm/` | 原始 GFA/缓存、标准 FASTA 和索引；phased 路线另有 combined.fa 和 ID 对照 |
| `05_linkage/` | 两套合并建 BWA 索引，逐 Hi-C 文库比对、去重、过滤和合并 BAM |
| `06_scaffolding/` | HapHiC 中间产物、contact_maps 和 delivery 交付包 |
| `07_finishing/` | 局部组装、quarTeT 补洞、端粒延伸、读段证据及坐标修改表；显式启用 |
| `08_polishing/` | NextPolish2 或明确的 DISABLED/SKIPPED_MISSING_WGS 状态、最终 genome.fa 和每份 hap FASTA |
| `09_assessment/` | 基础统计；`t2t/` 保存扩展验收、逐染色体表和待人工验收清单 |
| `10_report/` | 样本 manifest、Markdown、可离线打开的 HTML |

日志和资源统计位于样本的 logs/、benchmarks/，内部使用相同阶段编号。`report` 目标的批次汇总位于运行根目录 `10_reports/`。

HiFi-only 交付 primary.fa、alternate.fa、各自索引和标准 GFA；alternate 允许为空并显式记录。primary/alternate 不能直接解释成两套完整单倍型。

phased 路线的 `06_scaffolding/delivery/` 主要文件为：

| 文件 | 内容 |
|---|---|
| `chromosomes.fa` | 两套已挂载染色体候选合集 |
| `hap1.chromosomes.fa`、`hap2.chromosomes.fa` | 分开的两套候选 |
| `unplaced.fa` | 未挂载或未解决对象，保留源序列 |
| `*.fa.fai`、`*.chrom.sizes` | 索引与长度 |
| `assembly.agp`、`chromosomes.agp` | 全对象/仅染色体对象的原始 contig 坐标关系 |
| `chromosome_names.tsv` | 交付 ID 和上游对象对应 |
| `manifest.json`、`SHA256SUMS` | 状态、统计、来源和指纹 |
| `review.request.json` | 待审阅信息和版本绑定 |

## 6. 状态与质量边界

运行根目录的 config/plan.json 保存解析配置和工具/代码指纹；run_status.json 记录实际执行，dryrun_status.json 单独记录计划检查。`.snakemake` 用于依赖追踪和恢复。

COMPLETED 表示选定计算目标成功。判断交付还应查看样本 manifest：

| quality_level | 含义 |
|---|---|
| `CONTIG_ONLY_UNVALIDATED` | 基础 contig 结果，尚未独立认证质量 |
| `CHROMOSOME_CANDIDATE` | 数量达到配置目标，尚未人工审阅 |
| `CHROMOSOME_REVIEWED` | 应用了匹配该版本的有效审阅记录 |
| `INCOMPLETE_SCAFFOLDING` | 数量不符，保留部分交付及接触图并阻断最终成功 |
| `POLYPLOID_CHROMOSOME_CANDIDATE` | 四倍体各 chr×copy 组达到配置目标，尚未独立验证相位、人工审阅或认证 T2T |

交付校验检查 AGP 边界、方向、源序列每个碱基恰好覆盖一次，以及上游 FASTA 与 AGP 重建结果一致。混合 haplotype 的对象列入未解决集合。目标数量吻合、无 N 或高 N50 都不足以证明连接正确或 T2T。

染色体名按每套组内长度排序，如 hap1_Chr01，不表示已与参考 Chr01 同源对应；需要标准编号时须另做共线性命名。

## 7. 染色体挂载与审阅

Hi-C 使用 bwa mem -5SP、samblaster、-F3340 和 HapHiC 成对过滤，保留 query 相邻顺序，不做坐标排序。HapHiC 接收两套 GFA，phasing_weight=1，目标组数为两倍的 chromosomes_per_haplotype；所有未挂载序列保留。

当两个单倍型很相似时，合并参考可能使大量读段在两套上等分匹配、MAPQ=0。这时先检查比对得分、替代位置和组装覆盖信息，再显式设置 `scaffolding.hic_reference: per_haplotype`。流程会对每套保留原 MAPQ/NM 标准，以每套染色体数分别运行原生 HapHiC，保留两套原生输出并汇总 AGP/FASTA。BAM 汇总按序列名重建参考和 mate 字典，不直接拼接不同参考字典。统计单位为参考内的读对观察，同一物理分子可能贡献两套观察；汇总保留率以输入读对数乘参考数为分母。跨单倍型联系缺失是构造结果，不能用作分相正确的证据；结构变异区须结合原始图、读段和接触图审阅。该替代路线见 [HapHiC 作者的说明](https://github.com/zengxiaofei/HapHiC/issues/98)。映射或挂载失败时保留 `failed_attempts/` 中的原生工作目录和日志，成功后仅发布校验产物。

先查看 `06_scaffolding/contact_maps/` 的全局和逐染色体接触图，再结合原始 GFA/AGP 检查连接、方向和 hap 混合。若使用 Juicebox，可在项目计算任务中运行 HapHiC 保留的 `04.build/juicebox.sh` 生成交互材料；默认交付静态接触图。

编辑使用原始 contig 坐标的 scaffolds.raw.agp，保留所有拆分片段和未挂载对象。在运行目录外保存修改 AGP，并复制 review.request.json 填写实际 reviewer、reviewed_at、status=APPROVED、review_agp_sha256 和 chromosome_objects；输入 FASTA 和自动 AGP 的哈希必须匹配。

在原配置的 scaffolding 设置中同时提供 review_agp 和 review_manifest 的外部路径，再使用原 run-dir 恢复。入口拒绝过期哈希、未批准记录、序列丢失/重复、越界坐标和未知方向。审阅配置更新保持在下游，不重做未变化的联合组装。没有实际审阅时保持 NOT_REVIEWED。

### 染色体组装后的完善与验收

工具选择：07 使用 hifiasm 局部组装、quarTeT GapFiller 和 tidk；08 使用 NextPolish2；09 使用 tidk、Merqury/meryl、Inspector 和 BUSCO。工具选择与安装来源见[环境管理](../environment.md#t2t-完善工具)。不同时叠加多个补洞器，不将四倍体交给二倍体 purge 流程。

在已有项目配置中增加：

```yaml
defaults:
  finishing:
    enabled: true
    gap_closing: true
    telomere_completion: true
    telomere_motif: TTTAGGG  # 项目明确指定，不能对所有物种硬编码
    source_delivery: null   # null 表示在本次 DAG 的 06 交付之后继续
    threads: 16
    mem_mb: 64000
  polishing:
    enabled: true
    on_missing_wgs: skip    # 也可设 error，缺普通 WGS 时立即阻断
    threads: 16
    mem_mb: 64000
  assessment:
    read_mapping: true      # Inspector，使用最终序列重新比对 HiFi
    kmer: true              # Merqury；无普通 WGS 时明确 SKIPPED_MISSING_WGS
    busco: false            # 有适合物种的离线 lineage 后再设 true
    busco_lineage: null     # 离线目录，必须含 dataset.cfg
    threads: 16
    mem_mb: 64000
    min_qv: 60              # 项目筛查阈值，不是通用 T2T 定义
    min_busco_complete: 95
```

完整调用使用 `run assembly --config ... --run-dir ... --cores 16 --mem-mb 256000 --target t2t --env-prefix ...`。单个外层任务中的 local Snakemake 执行模型不变。`--target chromosomes` 仍止于 06，`--target report` 在新建完整分析中包含已启用的完善模块。

**接续冻结的旧组装**：复制原配置至项目新配置目录，开启上述模块，将 `finishing.source_delivery` 设为旧运行 `06_scaffolding/delivery/manifest.json` 的绝对路径；另用新 run-dir，先 `--target t2t --dry-run`。07 校验交付清单、每份 FASTA 的 SHA256、拷贝间序列一致性和目标数量，然后直接使用旧 06；不会调用旧 04–06。不要将正在写入的旧产物当成已完成输入，也不要直接修改旧快照。源样本和交付不一致、哈希改变或 INCOMPLETE_SCAFFOLDING 会阻断。

07 逐位点执行以下操作：

1. 对全部染色体拷贝和未挂载序列联合进行竞争比对，仅选取主比对、MAPQ≥30 的局部 HiFi/ONT 读段。
2. 用 hifiasm 对该位点局部组装；HiFi 种子不足默认 3 条时保留未解决状态。
3. 内部缺口要求两侧默认各 2 kb 锚定均唯一、方向一致、比对一致性≥99%。端粒延伸要求唯一末端锚定及延伸序列含近末端端粒阵列。
4. 将原始读段重比对至完整候选桥接序列，默认至少 3 个不同读段 ID 覆盖整个替换段及两侧锚定；大删除不算跨越证据。quarTeT 仅对通过证据筛查的单拷贝单缺口草图执行 fill，禁用 join。
5. 保存 accepted/rejected 原因、候选、局部图、比对和支持读段；不能解决的缺口不删除 N 来凑 gap-free。

这些阈值是保守筛查，不能独立证明多倍体相位。短于最小锚定、长于 `max_span`（默认 1 Mb）、无 HiFi 局部种子、重复路径歧义或缺少全跨越读段的位点保持未解决。hifiasm 局部路径无法覆盖的端粒不做任意重复延伸。`edits.tsv` 使用 0-based、右端不包含坐标，记录从 06 到 07 的替换；07 后原始 AGP 不再代表修改后的序列，不冒用旧 AGP。

08 仅接收 readsets 表中 `technology=wgs`、`format=fastq` 的普通成对基因组短读段，列 `file1/file2`。不使用 Hi-C 代替。fastp 后构建 21/31-mer yak 数据，NextPolish2 使用 MAPQ≥30 的 HiFi 比对，并以整序列比对检查输出 ID、拷贝数、方向和同拷贝覆盖。不满足条件时报错；缺少 WGS 时按配置记录 SKIPPED_MISSING_WGS 或阻断，不声称完成纠错。08 的 `genome.fa` 与 `hap1.fa` 等始终对应同一版本。

09 对最终 FASTA 重新运行 tidk，逐染色体列出 N/其他模糊碱基和两端近末端阵列。严格的末端筛查支持 motif 的循环移位及反向互补，内部重复或有模糊碱基不能计作“无缺口且两端端粒”。Merqury 分别报告总组装及每份 copy 的 QV；copy 的 k-mer 完整性不能当成全基因组完整性，无亲本 hap-mer 时不输出亲本相位结论。记录评估 WGS 是否也用于纠错。Inspector 只作为读段错误证据，不自动修改多倍体结构。BUSCO 离线逐 copy 执行，不将 duplicated 直接判作冗余。

结果在 `07_finishing/`、`08_polishing/`、`09_assessment/t2t/`、`10_report/t2t/`。状态为 INCOMPLETE、INCOMPLETE_QUALITY 或 CANDIDATE_NEEDS_REVIEW；最后一种仅说明计数、无缺口和末端阵列筛查符合要求。缺失 QV、禁用 BUSCO、独立相位、污染和复杂重复拷贝验收均单独记录，不能当作 PASS。`review.request.json` 绑定最终 FASTA 和 assessment 哈希，列出最终 Hi-C 审阅、相位、污染、着丝粒/rDNA 拷贝及修改位点审阅；程序不提供仅改状态即认证 T2T 的路径。

## 8. 恢复与排错

同一方案沿用配置、样本选择和 run-dir；入口自动启用 rerun-incomplete。新方案用新目录。变更参数、输入或相关代码可能使产物失效；保留原日志、manifest 和 `.snakemake`，检查对应编号阶段的日志。阶段参数分别记录相关源码指纹：仅修改接触图脚本不使输入准备、组装或比对失效；修改共享辅助函数会影响实际使用它的阶段。

旧源码使用综合任务脚本，改用当前源码后规则与指纹已变化，不能保证旧运行目录完全无重算。正在执行的项目继续使用其冻结快照；迁移已完成或运行中的真实项目需单独检查 DAG，不能直接替换源码后恢复。

| 问题 | 处理 |
|---|---|
| 缺输入、未知样本或 profile | 核对 TSV、相对路径和样本选择 |
| MD5 mismatch / FASTQ 解析失败 | 核对交付清单与源数据；修正原因后恢复 |
| 内存/Bloom 预算失败 | 协调实际申请和流程预算，按材料与 benchmark 调整 |
| 同目录已有启动器 | 检查旧任务及子进程，避免并发写入 |
| 强制终止后有 Snakemake 锁 | 确认旧任务全部退出后，使用同环境和计划执行 unlock，再恢复 |
| 交付数量不符 | 查看接触图及未解决对象，不能靠删小 contig 或强制连接凑数 |
| 历史迁移目录拒绝运行 | 读取已完成结果；真正的新分析另建 run-dir，见 [迁移说明](../development.md) |

旧报告可能留在磁盘，判断本次状态时核对 plan_id、assembly_id 和 run_status，不能仅以文件存在为依据。
