# 互作转录组分析 | Dual RNA-seq Analysis

一句话理解：**把混在一个样品里的「宿主 + 病原体」两类测序读段自动拆开，分别算出每个物种每个基因的表达量（FPKM/TPM），输出两套表达矩阵**，解决「一个样品里两个物种的 RNA 混在一起、没法分开定量」的问题。

## 功能概述 | Overview

- 双物种（如宿主 host 与病原体 pathogen）转录组同时分析，一次运行产出两套独立结果
- 核心是「读段物种分类」：每条读段同时比对到两个参考基因组，按 MAPQ 高低判定归属
- 对分类后的读段分别做 StringTie 定量，提取 FPKM 与 TPM
- 生成两个物种各自的表达矩阵与比对统计报告
- 可选从分类结果 BAM 反向提取 FASTQ，供下游二次分析
- 断点续传：索引、BAM、定量结果已存在时自动跳过

## 快速开始 | Quick Start

```bash
biopytools dual-rnaseq --species1-name host --species1-genome host.fa --species1-gtf host.gtf --species2-name pathogen --species2-genome pathogen.fa --species2-gtf pathogen.gtf -i ./fastq_data -o results
```

最小输入：两个物种的基因组 FASTA + GTF 注释各一份，加一个双端 FASTQ 目录（或样本信息文件）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 读段 read | 测序仪吐出来的一小段一小段序列，像拆散的书页碎片 |
| 双端测序 | 每条 DNA 片段从头尾各读一次，得 R1/R2 两份，像一张纸正反面各拍一张 |
| 参考基因组 | 该物种「标准答案」的完整 DNA 序列，用来对照 |
| GTF 注释 | 标注「基因组上哪些位置是基因、外显子」的坐标说明书 |
| 比对 alignment | 把每条读段放到参考基因组上找它最像的位置 |
| MAPQ | 比对质量分，越高表示「这条读段放在这个位置越有把握」 |
| 物种分类 | 判断一条读段到底来自宿主还是病原体 |
| FPKM / TPM | 两种「基因表达量」的计量单位，都做了深度和基因长度校正，数值越大表达越强 |
| StringTie | 负责「数每个基因有多少读段」并算出表达量的软件 |
| 表达矩阵 | 一张表：行是基因，列是样本，格子是表达量 |

## 输入 | Input

### 两个物种的参考文件

- 物种 1（如宿主）与物种 2（如病原体）各需要：基因组 FASTA 文件 + GTF 注释文件
- 基因组 FASTA：标准格式（`>` 开头的序列名 + 序列行）
- GTF 文件：供提取剪接位点与外显子，用于构建 HISAT2 索引

### 测序数据（二选一）

1. **FASTQ 目录**：目录内放成对的双端 FASTQ 文件，按 `-p/--pattern` 命名（默认 `*_1_clean.fq.gz`，`*` 是样本名，R2 会自动替换为 `_2_clean.fq.gz`）
2. **样本信息文件**：制表符分隔三列 `样本名\tR1路径\tR2路径`，一行一个样本

```text
sampleA    /data/sampleA_1_clean.fq.gz    /data/sampleA_2_clean.fq.gz
sampleB    /data/sampleB_1_clean.fq.gz    /data/sampleB_2_clean.fq.gz
```

## 参数说明 | Parameters

### 两个物种的参考文件 | Species references

**通俗理解|In plain words:** 这是「标准答案」——程序要知道宿主和病原体各自的基因组长什么样、基因在哪。**名字（name）会直接用作输出文件名前缀**，建议用简单英文如 host、pathogen，别用空格或特殊字符。

- `--species1-name` / `--species1-genome` / `--species1-gtf`：物种 1 的名称、基因组 FASTA、GTF 注释
- `--species2-name` / `--species2-genome` / `--species2-gtf`：物种 2 的名称、基因组 FASTA、GTF 注释

### 输入与输出 | Input & output

**通俗理解|In plain words:** `-i` 告诉程序去哪找测序数据（目录或样本清单文件），`-o` 是结果放哪。`-p` 是「文件名长什么样」的模板——你的文件命名跟默认不一致时才需要改。**默认命名 `*_1_clean.fq.gz` 是标准质控产物命名，一般不用动。**

- `-i/--input`：输入 FASTQ 目录或样本信息文件
- `-o/--output-dir`：输出目录（默认 `./dual_rnaseq_output`）
- `-p/--pattern`：FASTQ 命名模式，`*` 是样本名占位符

### 分类参数 | Classification

**通俗理解|In plain words:** 这三项决定「读段怎么判归属」。`--min-mapq` 是质量门槛——调高更保守（只信最确定的比对，但可能丢一些），调低更宽松，默认 20 是常用值。`--no-unique-only` 会关掉「只保留高质量比对」的过滤（即不再要求 MAPQ 不低于 min-mapq），只在特殊诊断场景用。`--no-extract-fastq` 跳过「从分类结果反提取 FASTQ」这一步，想省时间/磁盘时可加。**三项默认值都经实践验证，一般不用动。**

- `--min-mapq`：最小比对质量阈值（默认 20，范围 0-60）
- `--no-unique-only`：关闭「仅保留高质量唯一比对」过滤（默认开启过滤）
- `--no-extract-fastq`：不提取 FASTQ（默认会提取）

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t/--threads` 是并行线程数，机器核多可以调大加速比对和定量，默认 12 对多数场景够用。

## 分析流程 | Pipeline

```text
两个物种基因组 + GTF
    |
    v
步骤1: 分别构建两个物种的 HISAT2 索引(含剪接位点/外显子)
    |
    v
步骤2: 解析输入样本(成对 FASTQ)
    |
    v
步骤3: 读段物种分类
   +- 每条读段分别比对到两个基因组
   +- 只在物种1 -> 归物种1
   +- 只在物种2 -> 归物种2
   +- 两边都有 -> 取 MAPQ 更高的那侧;相同 -> ambiguous
   +- 两边都没有 -> unassigned
    |
    v
步骤3b: 比对统计(每个样本各物种读段数与比例)
    |
    v
步骤4: 对每个物种的分类读段分别 StringTie 定量,提取 FPKM/TPM
    |
    v
步骤5: 合并成两个物种各自的表达矩阵 + 总结报告
```

### 分类判定规则 | Classification rules

| 读段比对情况 | 归属 |
|---|---|
| 仅比对到物种 1 | 物种 1 |
| 仅比对到物种 2 | 物种 2 |
| 两边都比对上，物种 1 的 MAPQ 更高 | 物种 1 |
| 两边都比对上，物种 2 的 MAPQ 更高 | 物种 2 |
| 两边都比对上且 MAPQ 相同 | ambiguous（模棱两可） |
| 两边都没比对上 | unassigned（未归属） |

## 输出 | Output

```text
output_dir/
+-- 01_index/                              # 两个物种的 HISAT2 索引
|   +-- {species1}_hisat2_index.*.ht2      # 物种1索引(多份分片)
|   +-- {species1}_normalized.fa           # 大写标准化后的基因组
|   +-- {species1}.ss / {species1}.exon    # 剪接位点/外显子
|   +-- {species2}...(同上)
+-- 02_classification/                     # 读段分类结果(按样本分目录)
|   +-- {sample}/
|   |   +-- {sample}_{species1}.bam        # 归物种1的读段
|   |   +-- {sample}_{species2}.bam        # 归物种2的读段
|   |   +-- {sample}_ambiguous.bam         # 模棱两可
|   |   +-- {sample}_unassigned.bam        # 未归属
|   +-- .temp/                             # 中间原始比对(可删)
+-- 03_alignment_statistics/
|   +-- mapping_statistics.tsv             # 比对统计(中英列表头)
|   +-- mapping_statistics_en.tsv          # 纯英文版(供下游脚本)
+-- 04_quantification/                     # 定量结果
|   +-- {species1}/{sample}.gtf            # StringTie 定量 GTF
|   +-- {species1}/{sample}_fpkm.txt       # 该样本 FPKM/TPM 表
|   +-- {species2}/...(同上)
+-- 05_expression_matrix/                  # 表达矩阵
|   +-- {species1}_matrix.txt              # 物种1合并表达表
|   +-- {species2}_matrix.txt              # 物种2合并表达表
+-- 06_extracted_fastq/                    # (默认开启)反提取的 FASTQ
|   +-- {sample}_{species1}_1_clean.fq.gz  # 等
+-- dual_rnaseq_summary.txt                # 总结报告
+-- dual_rnaseq_processing_*.log           # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 比对统计（`03_alignment_statistics/mapping_statistics.tsv`）

**通俗理解|In plain words:** 这是「一份样品里到底有多少读段属于宿主、多少属于病原体」的台账。表头为中文，另有一份纯英文版方便下游脚本。

关键列：总读段数、仅物种1、仅物种2、同时比对（both）、未比对（neither），以及各自的百分比。**「同时比对」比例高，说明两个基因组有较多相似序列（如保守基因），属正常**；「未比对」比例高则提示测序质量问题或还混有其他物种。

### 2. 表达量文件（`04_quantification/{species}/{sample}_fpkm.txt`）

**通俗理解|In plain words:** 每个基因/转录本一行，记录它在这个样本里的表达量。

```text
gene_id    transcript_id    cov    FPKM    TPM    sample
gene1      gene1.1          12.3   3.45    5.67   sampleA
```

- `cov`：覆盖度；`FPKM`、`TPM`：两种表达量单位，值越大表达越强
- 用途差异：比较同一样本内不同基因用 TPM，比较同一基因在不同样本间可用 FPKM（如需严格跨样本比较，建议下游再做归一化）

### 3. 表达矩阵（`05_expression_matrix/{species}_matrix.txt`）

**通俗理解|In plain words:** 把该物种所有样本的表达量表拼成一张长格式大表（一行 = 一个转录本在一个样本里的记录，含 `sample` 列）。可直接喂给 Excel/R 做差异分析。

### 4. 总结报告（`dual_rnaseq_summary.txt`）

记录本次运行的输入文件、参数、样本清单与输出位置，是「运行留档」，出问题先看它。

## 参数选择建议 | Parameter Guidance

- **`--min-mapq`**：默认 20。若宿主与病原体亲缘较近（基因组相似度高）、想让分类更保守，可调到 30 甚至 40；反之想保留更多读段可降到 10~15
- **`--no-unique-only`**：仅当你想关掉 MAPQ 门槛、观察「全量读段」的归属分布时用；正常分析不要加
- **`--no-extract-fastq`**：不需要下游对分类读段二次组装/比对时加，省磁盘和时间
- **`-p/--pattern`**：你的文件不是 `*_1_clean.fq.gz` 命名时才改，`*` 只允许一个，且必须能识别 R1/R2 标识

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--species1-name` | 必填 |  | 物种1名称｜Species 1 name (e.g., host) |
| `--species1-genome` | 必填 |  | 物种1基因组FASTA文件｜Species 1 genome FASTA file |
| `--species1-gtf` | 必填 |  | 物种1GTF注释文件｜Species 1 GTF annotation file |
| `--species2-name` | 必填 |  | 物种2名称｜Species 2 name (e.g., pathogen) |
| `--species2-genome` | 必填 |  | 物种2基因组FASTA文件｜Species 2 genome FASTA file |
| `--species2-gtf` | 必填 |  | 物种2GTF注释文件｜Species 2 GTF annotation file |
| `--input, -i` | 必填 |  | 输入FASTQ文件目录或样本信息文件｜Input FASTQ file directory or sample information file |
| `--output-dir, -o` | `./dual_rnaseq_output` | Path | 输出目录｜Output directory |
| `--pattern, -p` | `*_1_clean.fq.gz` |  | FASTQ文件命名模式｜FASTQ file naming pattern |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-mapq` | `20` | int | 最小mapping quality值｜Minimum mapping quality value |
| `--no-unique-only` | — |  | 不禁用非唯一比对｜Do not disable non-unique mappings |
| `--no-extract-fastq` | — |  | 不提取FASTQ文件｜Do not extract FASTQ files from BAM (default: extract) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--species1-name` | 必填 |  | 物种1名称｜Species 1 name (e.g., host) |
| `--species1-genome` | 必填 |  | 物种1基因组FASTA文件｜Species 1 genome FASTA file |
| `--species1-gtf` | 必填 |  | 物种1GTF注释文件｜Species 1 GTF annotation file |
| `--species2-name` | 必填 |  | 物种2名称｜Species 2 name (e.g., pathogen) |
| `--species2-genome` | 必填 |  | 物种2基因组FASTA文件｜Species 2 genome FASTA file |
| `--species2-gtf` | 必填 |  | 物种2GTF注释文件｜Species 2 GTF annotation file |
| `-i, --input` | 必填 |  | 输入FASTQ文件目录或样本信息文件｜Input FASTQ file directory or sample information file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-p, --pattern` | `*_1_clean.fq.gz` |  | FASTQ文件命名模式｜FASTQ file naming pattern (e.g., "*.R1.fastq.gz") |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--min-mapq` | `20` | int | 最小mapping quality值｜Minimum mapping quality value |
| `--no-unique-only` | — | store_false | 不禁用非唯一比对｜Do not disable non-unique mappings |
| `--no-extract-fastq` | — | store_false | 不提取FASTQ文件｜Do not extract FASTQ files from BAM (default: extract) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（`pandas`、`pysam` 用于统计与分类）
- HISAT2（`hisat2`、`hisat2-build`，含 `extract_splice_sites.py`、`extract_exons.py`）
- SAMtools（`samtools`）
- StringTie（`stringtie`）
- 软件通过 conda 环境自动检测并包装调用，无需手动指定环境

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
断点续传按「输出文件是否存在」判断。改过滤参数（如 `--min-mapq`）重跑前，先删除旧输出（尤其是 `02_classification/`、`04_quantification/`），否则会复用旧结果。

**Q2：两个物种的 name 可以随便起吗？**
name 会直接作为输出文件名前缀（如 `{species}_matrix.txt`、`{sample}_{species}.bam`）。建议用 host/pathogen 这类简单英文，避免空格、特殊字符，否则文件名难看且易出问题。

**Q3：为什么会有「同时比对到两个物种」的读段？**
两个物种存在同源基因（序列相似），读段在两边都能高质量比对。程序会把两边 MAPQ 相同、都达到阈值的归为 ambiguous。这类读段比例不高属正常。

**Q4：输入是目录还是样本信息文件？**
两者都行。目录模式按 `-p` 模式自动配对 R1/R2；若文件命名不规则，建议用「样本名、R1路径、R2路径」三列的制表符文件，最省心。

**Q5：运行日志在哪？**
在输出目录根下的 `dual_rnaseq_processing_时间戳.log`（不在 `99_logs/` 子目录）。它同时打印到 stdout（INFO 级）与 stderr（WARNING 及以上）。
