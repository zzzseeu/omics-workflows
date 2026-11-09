# Assembly 进阶操作

日常输入准备、步骤和运行命令见 [Assembly 手册](assembly.md)。本页用于四倍体组装、已准备读段的复用、染色体组装完善和 Hi-C 人工审阅。完整参数见 [defaults.yaml](../../assembly/defaults.yaml) 和 [schema.yaml](../../assembly/schema.yaml)，软件安装见[环境管理](../environment.md)。

## 1. 四倍体组装

本项目接入的是 hifiasm 原始 unitig 图结合 GPhase 的路线。hifiasm 本身已有多倍体组装研究，包括遗传图谱辅助的分相及实验性 Hi-C 分相；后者在四倍体马铃薯案例中仍有各份组装大小不均衡的问题，见[原始方法论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11214949/)。这些方法与本项目当前配置范围应分开理解。

必须同时有 HiFi、ONT UL 和 Hi-C。使用 `hifi_ul_hic` + `unitigs` + GPhase，不能把二倍体的两套 hap contig 直接扩成四份。沿用联合组装模板的两张 TSV，按实际材料编写以下 YAML；染色体数和资源为示例。

```yaml
schema_version: 1
samples_table: samples.tsv
readsets_table: readsets.tsv
defaults:
  assembly:
    mode: hifi_ul_hic
    representation: unitigs
    purge_level: 0
    threads: 16
    mem_mb: 64000
  scaffolding:
    enabled: true
    backend: gphase
    threads: 16
    mem_mb: 64000
    processes: 4
    correct_nrounds: 0
samples:
  sample_A:
    ploidy: 4
    chromosomes_per_haplotype: 12
```

GPhase 依赖须先安装到同一个 `omics_workflows` 环境。运行目录不能含空白字符；`gap_size` 须为 100，`purge_level` 只能为 0 或 null，`correct_nrounds` 须为 0，由 GPhase/YaHS 进行原生校正。四倍体人工修改 AGP 的导入尚不支持。

先检查跨 readset 的长读段 ID 唯一性，再用 hifiasm `--n-hap 4 -l 0 --ul` 生成原始 p_utg 图，随后用 HiFi/minimap2/popCNV 估计拷贝数。Hi-C 按物理文库 BWA/samblaster 比对并保留 MAPQ0；GPhase 分相挂载后，重新将 Hi-C 比对到最终 contig components 绘图。初次聚类默认 MAPQ≥1，scaffolding 使用 MAPQ≥0，接触图使用成对 MAPQ/NM 过滤。方法依据见 [GPhase 官方说明](https://github.com/panlab-bioinfo/GPhase)。

每个 chr×copy 组须恰有一条 scaffold。缺组、多条或未分类对象会保留产物和接触图，并以 `INCOMPLETE_SCAFFOLDING` 阻止分析报告成功。不能挑选最大 scaffold 凑齐目标染色体数，也不能按 `genome_size_bp` 截断序列。

四份序列为 `06_scaffolding/delivery/hap1.chromosomes.fa` 至 `hap4.chromosomes.fa`，同目录保留 `components.fa`、`unitig_components.fa`、两类 AGP 和 `contig_to_unitig.txt`。GPhase 允许图重叠裁剪和 collapsed unitig 复制；导出时核对复制序列、AGP 坐标与原生 FASTA，不能套用 contig 路线“每个原始碱基恰好一次”的校验。未分配 unitig 保留在 `unplaced.fa`，其中可能有图冗余。copy 标签仅在染色体组内有效，不表示父母本或亚基因组。

## 2. 复用已准备的 HiFi/ONT 数据

当同一批长读段已在旧运行中完成准备和 PASS QC，可复用旧 FASTQ，避免再次转换、复制和统计。原始 `file1`、`format` 和 MD5 保留，在 `readsets.tsv` 增加 `prepared_manifest` 列，填写旧运行的 `samples/<样本>/01_inputs/<readset>.json`；路径相对于 TSV，也可用绝对路径。YAML 顶层另设：

```yaml
prepared_input_cache: /project/results/assembly_inputs/prepared_reuse_cache
```

该目录须在代码仓库和新运行目录之外。运行前核对样本、readset、物理文库、测序技术、源文件、MD5 和转换方法；来源不符、FASTQ 缺失或 QC 未通过时会报错，不自动重新转换。旧布局或其他软件的清单不能仅凭文件名使用。

旧清单若没有转换后 FASTQ 的 SHA256，首次使用时在计算节点读取一次压缩文件，同时计算 SHA256、用 seqkit 解析并核对原 QC 的读段数和碱基数；这没有与 BAM 逐读段比对。新清单已有 SHA256 且文件属性未变时可直接采用，否则流式复核。并发使用同一来源时只校验一次。

后续使用会核对来源清单、QC 以及文件大小、时间戳和 inode。已校验文件发生变化会阻断分析；删除校验缓存后可重新检查，无需转换或复制。新 `01_inputs/*.fastq.gz` 通过软链接引用旧 FASTQ，02 沿用原 PASS QC；仅运行 `t2t` 不额外执行 02。07 至 09 仍需针对对应版本的组装序列重新比对长读段。

旧 FASTQ、准备清单和 QC 必须保留，原始源文件的属性也须可核对。清理或移动旧运行前，应检查这些引用关系；被引用的数据不能作为临时文件删除。

## 3. 补洞、纠错与扩展评估

只适用于染色体级组装。将以下设置合并到已有 YAML 的 `defaults`，不要重复添加第二个 `defaults`。`telomere_motif` 须按物种填写，资源按实际需要调整；BUSCO 启用后必须提供离线 lineage。

```yaml
defaults:
  finishing:
    enabled: true
    gap_closing: true
    telomere_completion: true
    telomere_motif: TTTAGGG
    source_delivery: null
    threads: 16
    mem_mb: 64000
  polishing:
    enabled: true
    on_missing_wgs: skip
    threads: 16
    mem_mb: 64000
  assessment:
    read_mapping: true
    kmer: true
    busco: false
    busco_lineage: null
    threads: 16
    mem_mb: 64000
    min_qv: 60
    min_busco_complete: 95
```

在[主手册](assembly.md#4-运行示例)的路径变量基础上，填写已完成组装的配置和新的运行目录：

```bash
"$WORKFLOW_ROOT/run" assembly \
  --env-prefix "$OMICS_ENV_PREFIX" \
  --config "$PROJECT_ROOT/config/finishing.yaml" \
  --run-dir "$PROJECT_ROOT/runs/finishing_001" \
  --cores 16 --mem-mb 256000 --target t2t --dry-run
```

检查后去掉 `--dry-run` 执行。运行目录不能含空白字符。`source_delivery: null` 表示先完成本次染色体挂载，再继续后处理；从已完成的旧组装接续时，将其改为旧运行 `06_scaffolding/delivery/manifest.json` 的绝对路径，使用新运行目录。流程核对样本、每份 FASTA 的 SHA256、拷贝间序列一致性及染色体数后，直接读取旧 06 结果，不重跑旧组装和挂载。未完成、仍在写入、序列改变或 `INCOMPLETE_SCAFFOLDING` 的结果不能使用，原有组装与代码副本应保留。

### 07 补洞与端粒延伸

长读段与全部染色体拷贝和未挂载序列联合竞争比对，仅采用主比对且 MAPQ≥30 的局部 HiFi/ONT 读段。hifiasm 局部组装默认至少需要 3 条 HiFi 种子；补洞要求两侧各至少 2 kb 的唯一锚定、方向一致且序列一致性≥99%，端粒延伸还要求近末端端粒阵列。

候选桥接序列须有至少 3 个不同读段 ID 跨越完整替换区及两侧锚定，大删除不算跨越证据。quarTeT 只对通过筛查的单拷贝单缺口执行 fill，禁用 join。短于最小锚定、长于 `max_span`（默认 1 Mb）、重复路径不明确或缺少支持读段的位点保留未解决状态；不删除 N 来制造无缺口序列，也不任意延伸端粒重复。

`07_finishing/` 保留候选、局部组装图、比对和支持读段，以及接受或拒绝修改的原因。`edits.tsv` 用 0-based、右端不包含的坐标记录修改。修改后须使用新 FASTA，旧 AGP 不再表示新序列坐标；这些筛查也不能独立证明多倍体相位。

### 08 碱基纠错

普通双端 WGS 在读段表中填写 `technology=wgs`、`format=fastq` 和 `file1/file2`。fastp 处理后构建 21/31-mer yak 数据，NextPolish2 结合 MAPQ≥30 的 HiFi 比对纠错，再检查序列 ID、拷贝数、方向和同拷贝覆盖。不会使用 Hi-C 代替 WGS，也不将四倍体按二倍体 purge 处理。

缺少 WGS 时，`polishing.on_missing_wgs: skip` 记录 `SKIPPED_MISSING_WGS`，设为 `error` 则报错；关闭纠错记录 `DISABLED`。`08_polishing/genome.fa` 与各份 `hap*.fa` 始终对应同一版本，跳过不表示已完成纠错。

### 09 扩展质量评估

| 评估                | 解释结果时注意什么                                                                                             |
| ------------------- | -------------------------------------------------------------------------------------------------------------- |
| tidk 与逐染色体检查 | 统计 N、其他模糊碱基和两端近末端阵列；支持 motif 循环移位及反向互补，内部端粒重复不算染色体末端                |
| Merqury/meryl       | 分别报告总组装和各 copy 的 QV、k-mer 完整性；单份 copy 的完整性不能当作全基因组完整性，记录 WGS 是否也用于纠错 |
| Inspector           | 对最终序列重新比对 HiFi，提供读段错误证据，不自动修改多倍体结构                                                |
| BUSCO               | 使用指定的离线 lineage 逐 copy 评估保守基因完整性，duplicated 不能直接判为冗余                                 |

结果位于 `09_assessment/t2t/`，报告位于 `10_report/t2t/`。`INCOMPLETE`、`INCOMPLETE_QUALITY` 或 `CANDIDATE_NEEDS_REVIEW` 表示不同筛查结果；最后一种只说明计数、无缺口和末端阵列符合要求。QV 60 和 BUSCO 95% 是示例筛查阈值，不是通用 T2T 定义。缺少 QV、禁用 BUSCO、独立相位、污染和复杂重复拷贝均须单独记录，不能视为通过；没有亲本 hap-mer 时不作亲本相位判断。

`review.request.json` 对应最终 FASTA 和评估结果，列出最终 Hi-C、相位、污染、着丝粒/rDNA 拷贝及修改位点的审阅要求。通过筛查后仍需人工验收，不能只修改状态就认定 T2T。

## 4. Hi-C 挂载与人工审阅

二倍体默认将两套 hap contig 合并比对：bwa mem -5SP、samblaster、-F3340 和 HapHiC 成对过滤，保留 query 相邻顺序，不做坐标排序。默认 MAPQ≥1、NM<3；HapHiC 使用两套 GFA、`phasing_weight=1`，目标组数为每套染色体数的两倍，未挂载序列保留。`scaffolding.correct_nrounds` 默认 0，按断点证据决定是否启用自动校正。

若两个单倍型过于相似，大量 Hi-C 读段可能等分匹配且 MAPQ=0。先检查比对得分、替代位置和组装覆盖，再考虑 `scaffolding.hic_reference: per_haplotype`：分别比对两套 hap contig、按原过滤标准各自运行 HapHiC，最后汇总 AGP/FASTA。方法依据见 [HapHiC 作者说明](https://github.com/zengxiaofei/HapHiC/issues/98)。

这种方法按每个参考内的读对统计，同一物理分子可能计入两套观察，汇总保留率分母为输入读对数乘参考数；跨单倍型联系缺失不能证明分相正确。每个参考默认至少保留 0.1% 的读对（`scaffolding.min_retained_fraction: 0.001`），只是防止信号耗尽后继续聚类，不是通用样本合格标准。结构变异区仍需结合原始图、读段和接触图检查。比对或挂载失败的原始工作目录和日志保留在 `failed_attempts/`。

### 审阅和导入二倍体 AGP

先看 `06_scaffolding/contact_maps/` 的全局及逐染色体接触图，结合 GFA/AGP 检查连接、方向和 hap 混合。需要 Juicebox 交互材料时，可在项目计算任务中运行 HapHiC 保留的 `04.build/juicebox.sh`。

1. 编辑使用原始 contig 坐标的 `scaffolds.raw.agp`，保留所有拆分片段和未挂载对象，另存到运行目录外。
2. 复制 `review.request.json`，填写真实的 `reviewer`、`reviewed_at`、`status=APPROVED`、`review_agp_sha256` 和 `chromosome_objects`；输入 FASTA 及自动 AGP 的哈希须与待审阅版本一致。
3. 在原配置的 `scaffolding` 同时填写外部 `review_agp` 和 `review_manifest` 路径，以原运行目录恢复。过期、未批准、序列丢失或重复、坐标越界、方向未知的记录会被拒绝。

有效审阅结果标为 `CHROMOSOME_REVIEWED`，没有实际审阅仍为 `NOT_REVIEWED`。导入修改仅更新相关下游结果，不重做未变化的联合组装；更换代码或共享分析脚本后，相关步骤仍可能需要重跑。二倍体导出核对原始序列每个碱基的覆盖及 FASTA/AGP 一致性，混合 haplotype 的对象保留为未解决序列。
