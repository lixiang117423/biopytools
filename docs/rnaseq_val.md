# 转录组验证注释 | Transcriptome Validation for Genome Annotation

一句话理解：**用二代/三代转录组数据检验基因组注释的转录本结构对不对，自动给每个转录本分级(哪些可信、哪些要改/要删、哪些要人工看)，并输出校正后的注释 GTF。**

## 功能概述 | Overview { #overview }

- 两条数据路线：二代(HISAT2 比对 + StringTie 组装)、三代(minimap2 比对 + FLAIR 组装)，可只用其一或两者都用
- GFFcompare 把组装结果与参考注释做结构比较，产出 class code
- 自动分级决策：高置信 / 三代支持 / 二代支持 / 带标记保留 / 需人工复核
- 输出校正后注释 GTF + 各类别 TSV + 汇总报告
- 断点续传：按输出文件存在性自动跳过已完成步骤

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools rnaseq-val -g genome.fa -a anno.gtf --sr-dir ./sr/ -o out
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 转录本(transcript) | 一个基因转录出的 RNA 全长，有固定的外显子排列 |
| 转录组 | 一个样本所有转录本的集合 |
| 剪接 / 外显子结构 | mRNA 由外显子拼接而成，中间的内含子被去掉 |
| 二代 / 三代测序 | 二代读长短(几百 bp)、三代读长长(几 kb 到几十 kb) |
| 链特异性(strandness) | 记录 RNA 来自哪条 DNA 链，RF/FR 是两种常见建库方向 |
| GTF | 转录本结构的标准注释格式 |
| TPM / 覆盖度 | 转录本的表达量 / 被测序支持的程度 |
| class code | GFFcompare 给「组装转录本与参考转录本的异同」打的分类(如 = 完全一致、j 新亚型、u 新基因) |
| 校正(correction) | 根据证据决定保留、修改还是删掉注释里的转录本 |

## 输入 | Input { #input }

- `--genome`：基因组 FASTA
- `--annotation`：参考注释 GTF
- `--sr-dir`(可选)：二代 clean reads 目录，自动检测配对 fastq(R1/R2 等命名)
- `--lr-dir`(可选)：三代 clean reads 目录，自动检测单端 fastq/fq

至少提供 `--sr-dir` 或 `--lr-dir` 之一，否则报错。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required { #parameters-required }

**通俗理解|In plain words:** `-g` 基因组、`-a` 注释 GTF、`-o` 输出目录，三者必填。

### 数据输入 | Data input { #parameters-data }

**通俗理解|In plain words:** `--sr-dir` 放二代配对 fastq，`--lr-dir` 放三代单端 fastq。`--sr-pattern` 在命名不规范时指定配对模式(如 `*_1.clean.fq.gz`)。`--lr-platform` 告诉程序三代数据是 PacBio HiFi 还是 ONT，影响比对和组装参数。

### 运行参数 | Run { #parameters-run }

**通俗理解|In plain words:** `--strandness` 是建库链方向(RF/FR/unstranded)，错了会直接影响比对和组装，务必和建库方式一致。`--steps` 可只跑其中几步(逗号分隔)，默认 all 全跑。`--max-intron` 是三代比对允许的最大内含子长度，一般不用动。`--sample-timeout` 是单样本超时，大样本可调大。

### StringTie 组装 | StringTie { #parameters-stringtie }

**通俗理解|In plain words:** 控制二代组装的灵敏度。`--st-min-cov` 越低越敏感(能抓到低表达转录本但也更容易出噪声)；`--st-min-junction` 和 `--st-min-isoform` 同理。默认值比较稳妥，一般不用动。

### 校正阈值 | Correction thresholds { #parameters-threshold }

**通俗理解|In plain words:** 决定「一个转录本算不算被数据支持」。`--min-cov` / `--min-tpm` 越高越严格(更多转录本被判为需人工复核)；`--min-tpm-sr-only` 是「只有二代支持」时更严的 TPM 门槛；`--min-junction-ont` 是 ONT 数据的最低 junction 数。想保留更多转录本就调低，想更严格就调高。

### 运行控制 | Run control { #parameters-control }

**通俗理解|In plain words:** `--force` 强制重跑(忽略断点续传)；`--dry-run` 演练；`-v` / `--quiet` 调日志。

## 分析流程 | Pipeline { #pipeline }

```text
输入 genome.fa + annotation.gtf + (sr_dir 和/或 lr_dir)
    │
    ├─ align_2nd：HISAT2 建索引 + 比对二代 reads → 02_align_2nd/*.sorted.bam
    ├─ align_3rd：minimap2 比对三代 reads → 03_align_3rd/*.sorted.bam
    ├─ assemble_2nd：StringTie 组装 + 合并 → 04_assemble_2nd/merged.gtf
    ├─ assemble_3rd：FLAIR correct + collapse → 05_assemble_3rd/.../collapse.isoforms.gtf
    ├─ compare：GFFcompare 比较组装 vs 参考 → 06_compare/merged.annotated.gtf
    ├─ correct：分级决策 + 导出分类 TSV + 校正后 GTF → 07_correction/
    └─ report：汇总报告 → 07_correction/summary.tsv
```

## 输出 | Output { #output }

```text
out/
├── 00_pipeline_info/software_versions.yml   # 软件版本记录
├── 01_index/                  # HISAT2 索引
├── 02_align_2nd/              # 二代 BAM + flagstat
├── 03_align_3rd/              # 三代 BAM + flagstat
├── 04_assemble_2nd/           # StringTie 单样本 GTF + merged.gtf
├── 05_assemble_3rd/           # FLAIR 组装(flair/<样本>/collapse.isoforms.gtf)
├── 06_compare/                # GFFcompare：merged.annotated.gtf / .stats / .tracking
├── 07_correction/             # 校正结果(核心)
│   ├── keep_high_confidence.tsv
│   ├── keep_lr_supported.tsv
│   ├── keep_sr_supported.tsv
│   ├── keep_with_flag.tsv
│   ├── manual_review.tsv
│   ├── corrected_annotation.gtf   # 校正后注释(核心产物)
│   └── summary.tsv               # 汇总报告
└── 99_logs/                   # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- `07_correction/corrected_annotation.gtf`：按证据保留/修改后的注释，可直接用于后续分析
- 各类别 TSV 表示转录本的置信度分级：`keep_high_confidence`(二代+三代都支持，最可信)、`keep_lr_supported`(三代全长支持)、`keep_sr_supported`(仅二代且 TPM 充足)、`keep_with_flag`(class j/k 有支持)、`manual_review`(证据不足，需人工看)
- `summary.tsv` 里 GFFcompare 的 Sensitivity / Precision 反映「组装结果与参考注释的重合程度」，越高越好
- 每类 TSV 的占比可直接看出注释质量：manual_review 占比过高说明注释与转录组证据差异大，需重点检查
- `06_compare/*.stats` 是 GFFcompare 原始统计，class code 分布可辅助判断差异来源(新基因 vs 新亚型)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 有二代又有三代：默认全跑，两个路线互相验证，`keep_high_confidence` 最可信
- 只有二代：自动跳过三代相关步骤，靠 StringTie + `keep_sr_supported` 判断
- 只想快速看差异不做完整流程：`--steps compare,correct,report`(但需要已有组装产物)
- 想更激进地保留转录本：调低 `--min-cov` / `--min-tpm`；想更严格则调高
- ONT 数据建议配二代数据(会用二代 junction 辅助 FLAIR 校正)，`--lr-platform ont`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-a, --annotation` | 必填 |  | 参考注释GTF文件｜Reference annotation GTF file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `--sr-dir` | — |  | 二代clean reads目录｜SR clean reads directory |
| `--lr-dir` | — |  | 三代clean reads目录｜LR clean reads directory |
| `--lr-platform` | `pacbio_hifi` | pacbio_hifi/ont | 三代测序平台｜LR platform |
| `-t, --threads` | `16` |  | 线程数｜Threads |
| `--strandness` | `RF` | RF/FR/unstranded | 文库链特异性｜Library strandness |
| `--steps` | `all` |  | 执行步骤(逗号分隔或all)｜Steps (comma-separated or all) |
| `--force` | — |  | 强制重新运行｜Force re-run |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA 文件｜Genome FASTA file |
| `-a, --annotation` | 必填 |  | 参考注释 GTF 文件｜Reference annotation GTF file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `--sr-dir` | — |  | 二代 clean reads 目录（自动检测配对 fastq）｜SR clean reads directory (auto-detect paired fastq) |
| `--sr-pattern` | — |  | 二代 fastq 自定义命名模式｜SR fastq custom pattern (e.g. *_1.clean.fq.gz) |
| `--lr-dir` | — |  | 三代 clean reads 目录（自动检测单端 fastq）｜LR clean reads directory (auto-detect single fastq) |
| `--lr-platform` | `pacbio` | pacbio/ont | 三代测序平台｜LR sequencing platform |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--strandness` | `RF` | RF/FR/unstranded | 文库链特异性｜Library strandness |
| `--max-intron` | `500000` | int | 最大内含子长度｜Max intron length |
| `--steps` | `all` |  | 执行步骤（逗号分隔或 all）｜Steps to run (comma-separated or all) |
| `--sample-timeout` | `21600` | int | 单样本超时时间(秒)｜Sample timeout in seconds |
| `--st-min-cov` | `5.0` | float | StringTie 最低覆盖度｜StringTie min coverage |
| `--st-min-junction` | `3` | int | StringTie 最低 junction reads｜StringTie min junction reads |
| `--st-min-isoform` | `0.1` | float | StringTie 最低亚型丰度比｜StringTie min isoform fraction |
| `--min-cov` | `5.0` | float | 最低覆盖度阈值｜Min coverage threshold |
| `--min-tpm` | `0.5` | float | 最低 TPM 阈值｜Min TPM threshold |
| `--min-tpm-sr-only` | `1.0` | float | 仅二代支持时最低 TPM｜Min TPM for SR-only support |
| `--min-junction-ont` | `5` | int | ONT 最低 junction reads｜Min junction reads for ONT |
| `-v, --verbose` | `0` | count | 详细模式 (-v, -vv)｜Verbose mode |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--force` | — | store_true | 强制重新运行｜Force re-run |
| `--dry-run` | — | store_true | 试运行模式｜Dry run mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

统一使用 conda 环境 `rnaseq_val`，内含：

- hisat2 / hisat2-build / hisat2_extract_splice_sites.py / hisat2_extract_exons.py
- samtools、minimap2、stringtie、flair、gffcompare
- bedtools / regtools(ONT 时提取二代 junction)、isoquant.py、multiqc

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持。索引、各样本 BAM、组装 GTF、GFFcompare 结果都按「输出文件存在性」跳过。想强制重跑加 `--force`。

**Q2：`--steps` 怎么用？**
传逗号分隔的步骤名(align_2nd / align_3rd / assemble_2nd / assemble_3rd / compare / correct / report)或 all。中间某步产物已存在时，只跑后续步骤即可。

**Q3：CLI 的线程默认值到底是多少？**
命令行帮助显示 `-t` 默认 16，但内部 main.py 默认是 12；两者不一致，实际「不传 -t」时按 12 线程跑。要精确控制就显式传 `-t`。

**Q4：`--lr-platform` 为什么写 pacbio_hifi？**
CLI 提供 `pacbio_hifi` / `ont` 两个选项，内部把 pacbio_hifi 映射为 pacbio(PacBio HiFi 的比对/组装预设)。选默认或 pacbio_hifi 都按 HiFi 处理。

**Q5：报「必须提供至少一种数据源」？**
要传 `--sr-dir` 或 `--lr-dir` 之一，不能只有基因组和注释。

**Q6：链特异性选错了会怎样？**
会导致比对/组装方向判断错误、下游结果偏差。务必和建库方式一致(RF 是 dUTP 建库最常见)。