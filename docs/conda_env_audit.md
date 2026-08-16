# Conda 环境依赖审计报告|Conda Environment Dependency Audit Report

> 生成日期|Generated: 2026-08-16
> 生成脚本|Generator: `scripts/audit_conda_envs.py`（只读扫描，不改动任何代码|Read-only scan, no code changes）
> 完整明细|Full details: `docs/conda_env_audit.csv`（871 行引用明细）

## 一、总览|Overview

| 指标|Metric | 数值|Value |
|---|---|---|
| 引用点总数|Total reference points | 871 |
| 环境名字符串（原始）|Raw env name strings | 134 |
| 占位符/无效名|Placeholders / invalid names | 11 |
| 真实环境名|Real env names | 123 |
| **归一化后独立环境|Canonical envs** | **112** |
| ├─ 被代码引用|Cited in code | 102 |
| └─ 仅文档提及|Docs-only | 10 |
| 涉及文件|Files involved | 340 |
| 涉及模块|Modules involved | 156 |
| 同环境多写法（别名组）|Alias groups | 10 |
| 版本号钉在环境名里|Version-pinned env names | 74（占 112 的 66%） |
| 同工具分布在多个环境|Tools split across envs | 13 组 |

**结论|Conclusion:** 仓库代码实际依赖 **约 102 个 conda 环境**（112 个独立环境名减去 10 个仅文档提及），
平均每 1.5 个模块就有一个专属环境；其中 2/3 的环境把软件版本号写进了名字。

## 二、引用方式分布|Reference Style Breakdown

| 引用方式|Style | 数量|Count | 说明|Note |
|---|---|---|---|
| `envs/<NAME>/bin/<tool>` 路径默认值 | 690 | 绝大多数，散落在各模块 config/main 默认参数里 |
| `conda run -n <NAME>` 字面量 | 80 | 硬编码环境名的命令包装 |
| `conda activate <NAME>` | 33 | 主要在文档/脚本示例 |
| `xxx_env: str = "NAME"` 配置字段 | 16 | 环境名作为配置字段，配合动态包装 |
| `conda run -n {var}` 动态包装 | 29 | 环境名运行时解析（`get_conda_env()` 或配置字段传入），不引入新环境 |

> 动态包装的 29 处（rnaseq/mixrace/kmc/braker/allhic/rmvp/pixy 等 20+ 个文件）是**好模式**：
> 环境名由工具路径自动检测（§13.4），其实际环境已被上表 690 处路径默认值覆盖，未重复计数。

## 三、环境清单（按代码引用数 Top 30）|Top 30 Envs by Code Citations

| # | 环境|Env | 代码引用|Code | 模块数|Mods | 主要工具|Tools |
|---|---|---|---|---|---|
| 1 | Population_genetics | 25 | 15 | admixture, bedtools, bwa, plink, vcftools |
| 2 | GATK_v.4.6.2.0 | 22 | 16 | gatk, samtools, bcftools, wgsim |
| 3 | RNA_Seq | 21 | 12 | hisat2, stringtie, fastp, gffread |
| 4 | singularity_v.3.8.7 | 19 | 10 | singularity |
| 5 | BUSCO_v.6.0.0 | 10 | 7 | busco |
| 6 | psvcp_v.1.0.1 | 10 | 3 | Rscript, python3, samtools, bedtools, nucmer |
| 7 | yahs_v.1.2.2 | 10 | 2 | yahs, juicer |
| 8 | signalp6 | 9 | 6 | signalp6 |
| 9 | EDTA_v.2.2.2 | 9 | 4 | EDTA |
| 10 | resistify_v.1.3.0 | 9 | 4 | resistify, hmmsearch |
| 11 | repeatmodeler_v.2.0.7 | 9 | 5 | RepeatModeler, BuildDatabase |
| 12 | selective_sweep | 8 | 4 | RAiSD, xpclr, bcftools, vcftools |
| 13 | haphic | 8 | 6 | haphic, samblaster |
| 14 | jcvi_v.1.5.7（别名 JCVI_v.1.5.6） | 8 | 3 | python(jcvi) |
| 15 | Augustus_v.3.5.0 | 8 | 4 | augustus, etraining, bam2hints |
| 16 | samplot_v.1.3.0 | 8 | 2 | samplot |
| 17 | sv_calling | 7 | 5 | bcftools, bedtools, samtools |
| 18 | telocomp | 7 | 5 | minimap2, samtools |
| 19 | purge_dups_v.1.2.6 | 7 | 3 | (purge_dups 家族) |
| 20 | Rqtl | 7 | 2 | mstmap |
| 21 | repeat_identiy | 7 | 5 | RepeatMasker |
| 22 | needle_v.1.0.3 | 7 | 2 | needle(EMBOSS) |
| 23 | transdecoder_v.5.5.0 | 7 | 2 | TransDecoder |
| 24 | fanc_v.0.9.23b | 6 | 4 | fanc |
| 25 | eggnog-mapper_v.2.1.15 | 6 | 5 | emapper.py, mmseqs |
| 26 | braker_v.3.0.8 | 6 | 7 | hmmscan, miniprot |
| 27 | cphasing（别名 cphasing_v.0.2.10） | 6 | 6 | bwa-mem2, minimap2 |
| 28 | ltr_retriever_v.3.0.1（别名 ltr_retriever / LTR_retriever_v.3.0.4） | 6 | 5 | LTR_retriever |
| 29 | canu_v.2.3（别名 canu） | 6 | 3 | canu |
| 30 | pairtools_v.1.1.3 | 6 | 4 | pairtools |

> 完整 112 行清单见 `docs/conda_env_audit.csv`（canonical_env 列）。

## 四、冲突风险信号|Conflict Risk Signals

> 合并环境时最大的坑（即「有一类软件会冲突」），从静态扫描可以提前锁定以下几类：

### 4.1 同工具分布在多个环境（合并时的直接冲突源）|Same Tool Split Across Envs

同一个工具存在多份、版本可能不一致，合并时必须二选一或统一版本：

| 工具|Tool | 分布环境|Envs |
|---|---|---|
| bcftools | GATK_v.4.6.2.0 / bcftools_v.1.22 / selective_sweep / sv_calling |
| samtools | GATK_v.4.6.2.0 / psvcp_v.1.0.1 / sv_calling / telocomp |
| minimap2 | Genome_dedup / cphasing / pan-blocks / telocomp |
| bwa | Population_genetics / bwa_env / bwa_v0.7.17 |
| bedtools | Population_genetics / psvcp_v.1.0.1 / sv_calling |
| vcftools | Population_genetics / selective_sweep |
| python | biopytools / deeptmhmm_v.1.0 / jcvi_v.1.5.7 |
| miniprot | braker_v.3.0.8 / miniprot_v.0.18 |
| kmtricks | biopytools / kmtricks_v.1.5.1 |
| nucmer | pan-blocks / psvcp_v.1.0.1 |
| Rscript | WGCNA_v.1.73 / psvcp_v.1.0.1 |
| gt | genometools_v.1.6.5 / ltr_harvest_parallel_v.1.2 |
| biopytools | biopytools / primer3_v.2.6.1 |

### 4.2 版本钉死环境名（升级即改代码）|Version-Pinned Env Names

112 个环境里 **74 个**把版本号写进环境名（GATK_v.4.6.2.0、BUSCO_v.6.0.0、yahs_v.1.2.2…），
软件一升级就要改代码默认值，且同环境新旧版本并存时会产生隐性冲突。
完整名单见 CSV；典型高风险组合：
- 旧版本/legacy 栈：`juicer_v.1.6`（内含 **python2** + java + matlock）、`signalp_v.3.0b`（老 SignalP 3）、`Rqtl`（老 R）
- 版本分裂：`miniprot_v.0.18` vs `miniprot_v.1.0.0`、`cphasing` vs `cphasing_v.0.2.10`、`LTR_retriever_v.3.0.1` vs `v.3.0.4`

### 4.3 混合栈环境（合并时的冲突高发区）|Mixed-Stack Envs

这些环境装了一整个工具链（R + Python + C 工具混装），最容易被依赖冲突卡死：

- `Population_genetics`：plink + vcftools + admixture + bwa + bedtools（多模块共享，动它影响面最大）
- `psvcp_v.1.0.1`：Rscript + python3 + samtools + nucmer + assemblytics
- `yahs_v.1.2.2`：yahs + juicer（Java）
- `RNA_Seq`：hisat2 + stringtie + fastp + gffread + 自带脚本

### 4.4 占位符/无效引用（代码质量欠账）|Placeholder / Invalid Refs

- `biopytools/nlr_annotator/config.py:27`、`main.py:150`、`cli/commands/nlr_annotator.py:73`：Java 路径注释写死 `~/miniforge3/envs/xxx/bin/java`（占位符 `xxx` 未落实）
- 其余占位符均为文档示例（`tool_env`/`your_env`/`custom_panman_env`）和测试断言（`ENV`/`None`/`e`/`BCFTOOLS_ENV`），无需处理。

### 4.5 仅文档提及、代码未用（可暂缓/确认是否遗留）|Docs-Only Envs

smudgeplot / swave_v.1.2 / busco_env / my_ltr_harvest / minigraph / centier / pbbam_v.2.4.0 / jcvi_custom / bwa_v0.7.17 / sra-tools —— 共 10 个，代码无引用，合并方案里不需要为它们留位置（需人工确认无历史遗留）。

## 五、模块依赖排行|Modules by Env Count（代码范围 Top 20）

| 模块|Module | 依赖环境数|Envs |
|---|---|---|
| cli（各命令包装器） | 35 |
| tests | 27 |
| common（配置模板/路径指南） | 15 |
| oomycete_anno | 10 |
| assembly_qc | 9 |
| braker | 7 |
| hic_heatmap | 6 |
| haphic / phyto_effector / gene_rnaseq_check | 各 5 |
| insert2locus / yahs / hifi_hic / mixrace / kmertools | 各 4 |

## 六、下一步|Next Steps

1. **超算对照检查**（本审计无法覆盖的部分）：把本报告的环境清单与实际环境对比，找出「代码引用了但环境不存在」的死依赖：

   ```bash
   conda env list | awk '{print $1}' | sort > ~/tmp/conda_envs_actual.txt
   # 与 docs/conda_env_audit.csv 的 canonical_env 列 diff
   ```

2. **冲突预检**：对 4.1/4.3 的合并候选组，在超算上用 `mamba create -n test_merge --dry-run` 让求解器先验冲突，再决定合并粒度（不建环境、只求解）。

3. **合并优先级建议**（Phase 1，需等第 1 步结果）：
   - 先合并「同工具多环境」的 13 组（消除版本分裂），再合并混合栈；
   - `juicer_v.1.6`（python2）、`signalp_v.3.0b`、`Rqtl` 等 legacy 栈大概率保留为「例外环境」；
   - 环境名去版本号，版本信息统一进 `00_pipeline_info/software_versions.yml`（§12.5）。

## 七、审计方式|Audit Method

`scripts/audit_conda_envs.py` 只读扫描 `.py/.md/.yml/.yaml/.sh/.R/.txt`，按 6 类模式提取引用
（env-path / conda-run-literal / conda-activate / env-field / env-const / conda-run-var），
去重后按「去版本号 + 小写」归一化合并别名。可随时重跑：
```bash
python3 scripts/audit_conda_envs.py .
```
