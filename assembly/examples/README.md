# Assembly 使用案例

运行方式见 [快速开始](../../docs/quickstart.md)，输入、参数与输出见 [Assembly 手册](../../docs/workflows/assembly.md)。

本目录只保存示例代码和配置，reads 与运行结果由项目保存。每个案例目录包含匹配的 `config.yaml`、`samples.tsv` 和 `readsets.tsv`，模块顶层无需另放一套 `*.example.*`。

- `hifi_only/`：通用两样本配置。sample_A 包含两个虚拟文库（BAM + FASTQ），sample_B 包含一个 FASTQ 文库。复制到项目后，需提供对应数据或修改输入路径，并按真实材料评估资源；默认资源仅作示意。
- `hifi_ul_hic/`：二倍体联合组装配置模板；按实际材料填写 HiFi、ONT UL、成对 Hi-C、倍性和染色体数。
- `prepare_demo.py`：在仓库外的新目录创建小型合成数据项目，同时复制并调整配置、填入 MD5。通过统一环境的 Python 执行，不联网、不安装包、不提交任务。

合成示例用于检查流程能否运行，不能用于设定真实材料参数或评价 T2T。首次生成时，项目目录必须尚不存在；后续恢复分析使用已有运行目录，无需重新生成数据。

迁移后的生成器、统一接口执行和无重算恢复已由作业 `77958207`（DONE）验证，记录见 [开发与维护](../../docs/development.md)。历史验收记录保存在维护档案。
