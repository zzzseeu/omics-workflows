# Assembly 使用案例

运行方式见 [快速开始](../../docs/quickstart.md)，输入、参数与输出见 [Assembly 手册](../../docs/workflows/assembly.md)。

本目录是配置模板与案例的唯一入口，只包含示例代码和配置，不保存 reads 或运行结果。每个案例目录包含匹配的 `config.yaml`、`samples.tsv` 和 `readsets.tsv`；模块顶层不再另放三份 `*.example.*`。

- `hifi_only/`：通用两样本配置。sample_A 包含两个虚拟文库（BAM + FASTQ），sample_B 包含一个 FASTQ 文库。复制到项目后需提供对应数据或修改输入路径；默认资源是示意，需按真实材料评估。
- `hifi_ul_hic/`：二倍体联合组装配置模板；按实际材料填写 HiFi、ONT UL、成对 Hi-C、倍性和染色体数。
- `prepare_demo.py`：在仓库外的新目录创建小型合成数据项目，同时复制并调整配置、填入 MD5。通过统一环境的 Python 执行，不联网、不安装包、不提交任务。

合成示例只检验流程可运行，不用于设定真实材料参数或评价 T2T。首次生成要求项目目录不存在；后续分析恢复使用已有运行目录，不重新生成数据。

迁移后生成器、统一接口执行和无重算恢复已通过作业 `77958207`（DONE）验证；记录入口见 [开发与维护](../../docs/development.md)。历史验收保存在维护档案。
