# 基因组组装质量评估 | Genome Assembly Quality Control

一句话理解：**给一份基因组组装做一次「全面体检」**——只靠基因组本身看完整度（BUSCO）和重复序列组装质量（LAI），再用测序 reads 看准确度（QV）和比对率（mapping），最后汇总成一份 HTML 报告和一张发表用表格。

## 功能概述 | Overview { #overview }

- 七步流水线：基因组统计 → BUSCO 完整度 → LAI 指数 → QV 质量值 → NGS 比对 → 三代比对 → 报告
- 核心评估（只靠基因组本身）：BUSCO（默认开）、LAI（默认关，EDTA 耗时长）
- 可选评估（靠外部 reads）：QV（k-mer 光谱）、NGS mapping、三代 mapping（提供 reads 后自动启用）
- 输出 HTML 综合报告 + TSV/Excel 发表用表格
- 支持断点续传：各步骤输出已存在时自动跳过，方便失败后重跑

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools assembly-qc --genome genome.fa --lineage embryophyta_odb10 --ngs-reads ./illumina_reads --long-reads ./hifi_reads --long-read-type hifi -o qc_results
```

最小运行只需 `--genome` 和 `--lineage`（跑基因组统计 + BUSCO）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| BUSCO | 用一套「几乎所有物种都该有的核心基因」当探针，看组装里找回了多少；比例越高越完整 |
| LAI 指数 | 衡量「重复序列（尤其 LTR 转座子）装得对不对」的分数；越高说明重复区也装得好 |
| QV 值 | 组装「准确度」打分，QV40 约等于每 1 万碱基错 1 个，越大越准 |
| 比对率（mapping rate） | 测序 reads 能成功贴回基因组的比例，贴不上说明组装缺东西或错了 |
| 覆盖度（coverage） | 基因组上有多少比例位点被 reads 覆盖到 |
| N50 | 序列从长到短累加到总长 50% 时那条序列的长度，越大组装越连贯 |
| GC 含量 | 基因组里 G+C 占的比例，物种特有、常用于检查是否混样 |
| LTR 转座子 | 基因组里一类会「自我复制搬家」的重复元件，植物基因组尤其多 |
| EDTA | 一个识别重复序列的软件，LAI 评估靠它先找 LTR |

## 输入 | Input { #input }

- `--genome`：基因组 FASTA（必需）
- `--lineage`：BUSCO 谱系名（如 `embryophyta_odb10`）或完整路径（必需，哪怕 `--skip-busco`）
- `--ngs-reads`：Illumina 短读文件或目录（可选，用于 QV + NGS 比对；文件名默认匹配 `_1.clean.fq.gz`）
- `--long-reads`：三代 reads 文件或目录（可选，用于三代 QV + 比对；类型由 `--long-read-type` 指定）

目录下按命名规则自动发现样本：NGS 找 `*_1.clean.fq.gz` 并推断对应的 `_2` 文件；三代按 `*.fq.gz` 收集。

## 参数说明 | Parameters { #parameters }

### 核心评估开关 | Core evaluation switches

**通俗理解|In plain words:** BUSCO 默认开；LAI 默认关（因为它要跑 EDTA，对植物大基因组可能几小时到几天）。想加 LAI 就 `--enable-lai`；确认不想要 BUSCO 才 `--skip-busco`。**LAI 只在重复序列含量足够时才有意义，含量太低会返回「不适用」。**

相关参数：`--skip-busco`、`--enable-lai`、`--lai-full-mode`、`--lai-repeatmasker-species`。

### 可选评估与 reads | Optional evaluations and reads

**通俗理解|In plain words:** QV、NGS 比对、三代比对默认都开启，但「没米下锅」——只有给了对应的 reads 才会真正跑。给了 `--ngs-reads` 就自动跑 NGS 的 QV 和比对，给了 `--long-reads` 就自动跑三代的。**提供哪种 reads 就评估哪种，缺数据时对应步骤会被跳过而不是报错。**

相关参数：`--enable-qv`、`--qv-kmer-size`、`--enable-mapping`、`--enable-long-read-mapping`、`--ngs-reads`、`--long-reads`、`--long-read-type`、`--mapping-pattern`。

### 报告输出 | Report output

**通俗理解|In plain words:** HTML 报告和发表用表格默认都生成；`--no-html` / `--no-table` 关掉对应产物；`--table-format` 选表格是 TSV、Excel 还是两者都要。**默认 `both` 最省心，只要纯文本就 `--table-format tsv`。**

相关参数：`--no-html`、`--no-table`、`--table-format`。

### 线程 | Threads

**通俗理解|In plain words:** `-t` 是全局线程数，各步骤串行执行、每个步骤都用满全部线程，调大一般更快。注意 LAI 用的 EDTA 有「内存随线程近线性增长」的坑，默认把 EDTA 线程封顶在 24，大基因组加高线程容易内存爆掉（88 线程实测 1.2 TB）。**一般不用动；要压 EDTA 内存可显式给 `--edta-threads`。**

相关参数：`-t, --threads`、`--edta-threads`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 FASTA (+ 可选 reads)
    │
    ├─ 步骤1 基因组统计（大小/N50/GC）
    ├─ 步骤2 BUSCO 完整度评估（默认开）
    ├─ 步骤3 LAI 指数评估（默认关，--enable-lai 开启）
    ├─ 步骤4 QV 质量值评估（有 reads 时）
    ├─ 步骤5 NGS mapping 评估（有 --ngs-reads 时）
    ├─ 步骤6 三代 mapping 评估（有 --long-reads 时）
    │
    ▼
步骤7 生成 HTML 报告 + TSV/Excel 发表表格
```

## 输出 | Output { #output }

```text
qc_results/
├── 01_busco_evaluation/               # BUSCO 结果（busco_output/ 内 short_summary.*.json/txt）
├── 02_lai_evaluation/                 # LAI 结果（EDTA 输出 + *.LAI）
├── 03_qv_evaluation/                  # QV 结果（reads.meryl + qv_result_*.qv）
├── 04_mapping_evaluation/             # NGS 比对（bam_files/ 内各样本 sorted.bam/flagstat/coverage）
├── 05_long_read_mapping_evaluation/   # 三代比对（bam_files/ 内各样本）
├── assembly_qc_report.html            # 综合 HTML 报告
├── assembly_qc_table.tsv              # 发表用表格（TSV）
├── assembly_qc_table.xlsx             # 发表用表格（Excel，table-format 含 xlsx 时）
└── assembly_qc.log                    # 运行日志
```

发表表格列：sample、size、n50、gc、busco、lai、qv_ngs、qv_long、qv_ngs_error、qv_long_error、ngs_mapping_rate、ngs_mapping_cov、long_mapping_rate、long_mapping_cov。

## 结果解读 | Interpreting Results { #interpreting-results }

- **BUSCO 完整度**：≥90% 良好，≥95% 优秀；`duplicated` 偏高提示组装可能有冗余（如未去单倍型冗余）
- **LAI**：≥15 优秀，10–15 良好，<10 重复区组装较差；返回「不适用」说明 LTR 含量低于 5%，不能据此判质量
- **QV**：≥40 优秀，≥30 良好，<20 较差（错误率约 >1%）
- **NGS / 三代比对率**：通常 >95% 算正常，明显偏低提示组装缺失或 reads 与组装不匹配；三代 HiFi 比对率略低于短读属正常
- **覆盖度（coverage_fraction）**：指「有多少比例位点被 reads 盖到」，越高说明组装被 reads 支持得越充分

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 快速自检：只给 `--genome --lineage`，跑统计 + BUSCO
- 完整植物基因组评估：加 `--enable-lai` + `--ngs-reads` + `--long-reads --long-read-type hifi`
- 只需发表表格字段、不要 HTML：`--no-html --table-format tsv`
- 大基因组 / 内存吃紧：给 `--edta-threads 8` 压低 EDTA 内存
- 断点续传默认开启，重跑会自动跳过已完成步骤；换参数重跑见 FAQ

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--lineage, -l` | 必填 |  | BUSCO数据集名称（如embryophyta_odb10）或完整路径｜BUSCO dataset name (e.g., embryophyta_odb10) or full path |
| `--output-dir, -o` | `./assembly_qc_output` | Path | 输出目录｜Output directory |
| `--sample-name, -s` | `genome_sample` |  | 样品名称｜Sample name |
| `--skip-busco` | — |  | 跳过BUSCO评估｜Skip BUSCO evaluation |
| `--enable-lai` | — |  | 启用LAI评估（默认禁用，EDTA流程耗时长）｜Enable LAI evaluation (disabled by default, EDTA is time-consuming) |
| `--lai-full-mode` | — |  | LAI完整模式（不使用-qq，运行blastn计算，用于种间比较）｜LAI full mode (no -qq, run blastn for interspecies comparison) |
| `--lai-repeatmasker-species` | `Viridiplantae` |  | RepeatMasker物种参数（EDTA失败回退时使用）｜RepeatMasker species for EDTA fallback |
| `--enable-qv` | — |  | 启用QV评估（默认启用）｜Enable QV evaluation (default: enabled) |
| `--qv-kmer-size` | — | int | k-mer大小（None表示自动选择）｜K-mer size (None for auto) |
| `--enable-mapping` | — |  | 启用NGS Mapping评估（默认启用）｜Enable NGS mapping evaluation (default: enabled) |
| `--enable-long-read-mapping` | — |  | 启用三代数据Mapping评估（默认启用）｜Enable long-read mapping evaluation (default: enabled) |
| `--mapping-pattern` | `_1.clean.fq.gz` |  | FASTQ文件匹配模式｜FASTQ file pattern |
| `--ngs-reads` | — | Path | NGS reads文件或目录（用于QV和mapping）｜NGS reads file or directory (for QV and mapping) |
| `--long-reads` | — | Path | Long-reads文件或目录（用于QV和mapping）｜Long-reads file or directory (for QV and mapping) |
| `--long-read-type` | `hifi` | ont/pacbio/hifi | Long-read数据类型｜Long-read data type |
| `--no-html` | — |  | 不生成HTML报告｜Do not generate HTML report |
| `--no-table` | — |  | 不生成表格｜Do not generate table |
| `--table-format` | `both` | tsv/xlsx/both | 表格格式｜Table format |
| `--threads, -t` | `12` | int | 线程数（自动分配给各子模块）｜Threads (automatically distributed to sub-modules) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--lineage, -l` | 必填 |  | BUSCO数据集名称（如embryophyta_odb10）或完整路径｜BUSCO dataset name (e.g., embryophyta_odb10) or full path |
| `--output-dir, -o` | `./assembly_qc_output` |  | 输出目录｜Output directory |
| `--sample-name, -s` | `genome_sample` |  | 样品名称｜Sample name |
| `--skip-busco` | — | store_true | 跳过BUSCO评估｜Skip BUSCO evaluation |
| `--enable-lai` | — | store_true | 启用LAI评估（默认禁用，EDTA流程耗时长）｜Enable LAI evaluation (disabled by default, EDTA is time-consuming) |
| `--lai-full-mode` | — | store_true | LAI完整模式（不使用-qq，运行blastn计算，用于种间比较）｜LAI full mode (no -qq, run blastn for interspecies comparison) |
| `--lai-repeatmasker-species` | `Viridiplantae` |  | RepeatMasker物种参数（EDTA失败回退时使用）｜RepeatMasker species for EDTA fallback (default: Viridiplantae) |
| `--enable-qv` | — | store_true | 启用QV评估（默认启用）｜Enable QV evaluation (default: enabled) |
| `--qv-kmer-size` | — | int | k-mer大小（None表示自动选择）｜K-mer size (None for auto) |
| `--enable-mapping` | — | store_true | 启用NGS Mapping评估（默认启用）｜Enable NGS mapping evaluation (default: enabled) |
| `--enable-long-read-mapping` | — | store_true | 启用三代数据Mapping评估（默认启用）｜Enable long-read mapping evaluation (default: enabled) |
| `--mapping-pattern` | `_1.clean.fq.gz` |  | FASTQ文件匹配模式｜FASTQ file pattern |
| `--ngs-reads` | — |  | NGS reads文件或目录（用于QV和mapping）｜NGS reads file or directory (for QV and mapping) |
| `--long-reads` | — |  | Long-reads文件或目录（用于QV和mapping）｜Long-reads file or directory (for QV and mapping) |
| `--long-read-type` | `hifi` | ont/pacbio/hifi | Long-read数据类型｜Long-read data type |
| `--no-html` | — | store_true | 不生成HTML报告｜Do not generate HTML report |
| `--no-table` | — | store_true | 不生成表格｜Do not generate table |
| `--table-format` | `both` | tsv/xlsx/both | 表格格式｜Table format |
| `--threads, -t` | `12` | int | 线程数（自动分配给各子模块）｜Threads (automatically distributed to sub-modules) |
| `--edta-threads` | — | int | EDTA专用线程数（默认自动封顶24, 因EDTA内存随线程近线性增长, 88线程实测1.2TB易OOM）｜EDTA-only threads (default auto-cap 24, since EDTA memory scales ~linearly with threads; observed 1.2TB at 88 threads) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- BUSCO v6（conda 环境 `BUSCO_v.6.0.0`，数据集目录默认 `~/database/busco`）
- EDTA（`EDTA_v.2.2.2`）、LTR_retriever（`LTR_retriever_v.3.0.4`）、LTR_harvest（`ltr_harvest_parallel_v.1.2`）、LTR_FINDER_parallel（`ltr_finder_parallel_v.1.3`）——仅 LAI
- Merqury / meryl（`merqury_v.1.3`）——仅 QV
- BWA（`Population_genetics`）、samtools（`GATK_v.4.6.2.0`）、minimap2（`Genome_dedup`）、bedtools——比对评估；比对用工具走系统 PATH 直接调用
- Python 3（pandas、openpyxl 用于表格）

## 常见问题 | FAQ { #faq }

**Q1：报「BUSCO 数据集路径不存在」？**
默认在 `~/database/busco` 找数据集。把 `--lineage` 对应的数据集目录放到该路径下，或用 `--lineage` 传完整路径（程序会取目录名作谱系名）。

**Q2：NGS mapping 没发现样本？**
确认 reads 文件名匹配默认模式 `_1.clean.fq.gz`，且同一目录有对应的 `_2.clean.fq.gz`；命名不同就用 `--mapping-pattern` 指定。

**Q3：LAI 返回「不适用（LTR 含量低于 5%）」？**
这是正常结果，不是错误——重复序列太少时 LAI 指标无意义，应改用其他指标评估该基因组。

**Q4：换参数重跑为什么结果没变？**
断点续传按输出文件存在性判断。换 BUSCO 谱系、QV k 值、mapping 模式等参数后，需先删对应步骤目录（如 `03_qv_evaluation`）里的旧产物，否则会复用旧结果。

**Q5：Excel 表格生成失败？**
缺 `pandas` 或 `openpyxl`。装 `pip install pandas openpyxl`，或改用 `--table-format tsv`。
