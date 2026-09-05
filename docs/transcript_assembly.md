# 转录本组装 | Transcript Assembly (FASTQ/BAM to GFF3)

一句话理解：**把 RNA-seq 比对结果（或直接给 BAM/FASTQ）拼成一个个「基因结构」，输出标准 GFF3 文件，并可顺带提取转录本序列、预测编码区（CDS）**——是做无参考转录组、基因结构注释的基础工具。

## 功能概述 | Overview

- 基于 HISAT2 + StringTie 的转录本组装流程
- 输入二选一：FASTQ 目录（走 HISAT2 比对）或 BAM 文件（直入，跳过比对）
- 支持短读段与长读段（长读段走 StringTie `-L` 模式）
- 主输出 GFF3 基因结构；可选输出转录本 cDNA（`--transcripts`）、预测 CDS（`--predict-cds`）
- 可选参考注释引导组装（`--guide-gff`）
- 断点续传：各步骤按输出文件存在性跳过，`--force` 强制重跑

## 快速开始 | Quick Start

```bash
biopytools transcript-assembly -b sample.bam -o ./out
```

最小输入：一个 BAM 文件 + 输出目录（BAM 模式无需基因组）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 转录本 transcript | 一个基因转录出的 RNA 序列，同一基因可能有多种剪接形式 |
| 基因结构 | 一个转录本由哪些外显子（片段）按什么顺序拼成 |
| GFF3 / GTF | 描述「基因组上基因、外显子、CDS 坐标」的标准格式 |
| 组装 assembly | 从比对结果推断出转录本的完整结构 |
| CDS | 编码蛋白的那段序列 |
| TransDecoder | 从转录本里找最可能编码蛋白的开放阅读框（ORF）的软件 |
| 引导组装 guided | 用已有注释当「骨架」，让组装更贴合已知基因 |

## 输入 | Input

### 输入（二选一，互斥）

1. **FASTQ 目录**（`-i/--input`）：需要同时给 `-g/--genome`（走 HISAT2 建索引 + 比对）
2. **BAM 文件**（`-b/--bam`，可重复指定多个）：直入组装，无需基因组（除非要 `--transcripts`/`--predict-cds`）

### 其他输入

- `--guide-gff`：参考注释 GTF/GFF3（引导组装，可选）
- BAM 输入要求：坐标排序（`SO:coordinate`）、有索引（缺 `.bai` 会自动补建）

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i`（FASTQ）和 `-b`（BAM）二选一。用 FASTQ 必须给 `-g` 基因组；用 BAM 通常不用基因组，除非要额外输出转录本序列或预测 CDS。`--guide-gff` 有已知注释时提供，组装会更贴近已知基因。

- `-o/--output`：输出目录（必需）
- `-g/--genome`：基因组 FASTA（FASTQ 模式或 `--transcripts`/`--predict-cds` 时必需）
- `-i/--input`：FASTQ 目录（与 `-b` 互斥）
- `-b/--bam`：BAM 文件（可多个，与 `-i` 互斥）
- `--guide-gff`：参考注释（引导组装）

### 读长与附加输出 | Read type & extra output

**通俗理解|In plain words:** `--read-type` 告诉程序是短读还是长读；默认 `auto` 会按读长中位数自动判断（约 500 bp 为界），**一般选 auto 即可**。`--transcripts` 额外输出转录本序列 FASTA，`--predict-cds` 用 TransDecoder 预测编码区——两者都需要 `-g` 基因组。

- `--read-type`：`auto`/`short`/`long`（默认 `auto`）
- `--transcripts`：额外输出 transcripts.fa cDNA（需 `-g`）
- `--predict-cds`：TransDecoder 预测 CDS（需 `-g`）

### 样本与运行 | Samples & runtime

**通俗理解|In plain words:** `-p` 是 FASTQ 命名模式（`*` 为样本名）。`-t` 线程数。`--sample-timeout` 是单样本超时（默认 43200 秒 = 12 小时）。**默认值一般不用动。**

- `-p/--pattern`：FASTQ 命名模式（默认 `*_1_clean.fq.gz`）
- `-t/--threads`：线程数（默认 12）
- `--sample-timeout`：单样本超时秒数（默认 43200）

### 步骤控制与日志 | Step control & logging

**通俗理解|In plain words:** `-s/--step` 只跑某一步（1-7，调试/续跑用）；`--force` 强制重跑已完成步骤。`-v`/`--quiet` 控制输出多少。**正常跑全流程这些都无需指定。**

- `-s/--step`：运行指定步骤（1-7）
- `--force`：强制重跑已完成步骤
- `-v/--verbose` / `--quiet`：详细程度 / 静默

## 分析流程 | Pipeline

```text
FASTQ 模式:                       BAM 模式:
步骤1 hisat2-build 建索引          (跳过 1-3,校验+建索引+测读长)
步骤2 hisat2 --dta 比对            |
步骤3 samtools sort 排序+索引      |
    +-------------------------------+
    v
步骤4 stringtie 逐样本组装(长读 -L / 引导 -G)
    |
    v
步骤5 stringtie --merge 合并(仅多样本)
    |
    v
步骤6 gffread GTF -> GFF3(主输出) + 可选提取 cDNA
    |
    v
步骤7 (可选 --predict-cds) TransDecoder 预测 CDS
```

### 七步说明 | The 7 steps

| 步骤 | 内容 | 说明 |
|---|---|---|
| 1 | 构建 HISAT2 索引 | 仅 FASTQ 模式 |
| 2 | HISAT2 比对（`--dta`） | 仅 FASTQ 模式 |
| 3 | SAM 转排序 BAM + 索引 | 仅 FASTQ 模式 |
| 4 | StringTie 逐样本组装 | FASTQ/BAM 共用 |
| 5 | StringTie 合并 GTF | 仅样本数 > 1 |
| 6 | GFF3 输出（+ 可选 cDNA） | 主输出 |
| 7 | TransDecoder 预测 CDS | 可选（`--predict-cds`） |

## 输出 | Output

```text
output_dir/
+-- 00_pipeline_info/software_versions.yml  # 软件版本记录
+-- 01_hisat2_index/                        # HISAT2 索引(仅 FASTQ 模式)
+-- 02_hisat2_align/{sample}.sam            # 比对 SAM(仅 FASTQ 模式)
+-- 03_bam_sort/{sample}_sorted.bam(.bai)  # 排序 BAM(仅 FASTQ 模式)
+-- 04_stringtie/
|   +-- {sample}.gtf                       # 逐样本组装 GTF
|   +-- gtf_list.txt                       # GTF 列表(供合并)
+-- 05_merge/merged.gtf                    # 合并后 GTF(仅多样本)
+-- 06_gff3/merged.gff3                    # 主输出:基因结构 GFF3
+-- 06_transcripts/transcripts.fa          # (可选)cDNA 序列
+-- 07_transdecoder/                       # (可选)CDS 预测结果
|   +-- transcripts.fa.transdecoder.gff3
|   +-- transcripts_fa_transdecoder_genome.gff3
|   +-- transcripts.fa.transdecoder.pep / .cds
+-- 99_logs/pipeline.log                    # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 主输出 GFF3（`06_gff3/merged.gff3`）

**通俗理解|In plain words:** 最终成果——标准 GFF3 基因结构文件，记录组装出的转录本、外显子等坐标。可直接用于基因组注释、可视化（IGV）或下游分析。单 BAM 且未合并时，用该样本的 GTF 转换而成。

### 2. 转录本序列（`06_transcripts/transcripts.fa`，可选）

每条转录本的 cDNA 序列（用 gffread 从基因组提取），供同源性比对、功能注释、CDS 预测等使用。

### 3. CDS 预测（`07_transdecoder/`，可选）

- `transcripts_fa_transdecoder_genome.gff3`：映射回基因组坐标的 CDS 结构（gene/mRNA/CDS）
- `.pep`：预测的蛋白序列；`.cds`：编码区核酸序列
- 注意：本流程的 TransDecoder 为基础版，未使用 blastp/hmmscan 同源数据库，CDS 预测偏保守

### 4. 软件版本（`00_pipeline_info/software_versions.yml`）

记录 stringtie/gffread/hisat2/samtools 等工具的路径与版本，供论文方法复现。

## 参数选择建议 | Parameter Guidance

- **有 BAM 优先用 `-b` 直入**：跳过建索引和比对，省大量时间
- **长读段**：`--read-type long` 或用 `auto` 自动识别（读长中位数 ≥ 500 bp 判为长读，走 StringTie `-L`）
- **有参考注释**：加 `--guide-gff` 引导组装，结果更贴近已知基因结构
- **要蛋白序列**：加 `--predict-cds`（需 `-g`），会额外跑 TransDecoder
- **单样本**：样本数 = 1 时自动跳过合并，直接用该样本 GTF 出 GFF3

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--genome, -g` | — |  | 基因组FASTA(FASTQ模式或--transcripts时必需)｜Genome FASTA (required for FASTQ mode or --transcripts) |
| `--input, -i` | — | Path | 输入FASTQ文件目录(与-b互斥)｜Input FASTQ dir, mutually exclusive with -b |
| `--bam, -b` | — | Path | 输入BAM文件(可多次-b,与-i互斥)｜Input BAM file(s), repeatable, mutually exclusive with -i |
| `--guide-gff` | — | Path | 参考注释GTF/GFF3(-G guided)｜Reference annotation for guided assembly |
| `--read-type` | `auto` | auto/short/long | 读长类型｜Read type |
| `--transcripts` | — |  | 额外输出transcripts.fa(需-g)｜Also output cDNA (needs -g) |
| `--predict-cds` | — |  | TransDecoder预测CDS(需-g,输出gene/mRNA/CDS)｜TransDecoder CDS prediction (needs -g) |
| `--pattern, -p` | `*_1_clean.fq.gz` | str | FASTQ文件命名模式｜FASTQ file naming pattern (* is sample name placeholder) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--sample-timeout` | `43200` | int | 单个样本处理超时时间（秒）｜Sample processing timeout in seconds |
| `--step, -s` | — | 1/2/3/4/5/6/7 | 运行指定步骤｜Run only specified step |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--force` | — |  | 强制重新处理已完成的步骤｜Force re-process completed steps |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-g, --genome` | — |  | 基因组FASTA(FASTQ模式或--transcripts时必需)｜Genome FASTA (required for FASTQ mode or --transcripts) |
| `-i, --input` | — |  | 输入FASTQ文件目录｜Input FASTQ file directory |
| `-b, --bam` | — |  | 输入BAM文件(一个或多个,与-i互斥)｜Input BAM file(s), mutually exclusive with -i |
| `--guide-gff` | — |  | 参考注释GTF/GFF3(-G guided)｜Reference annotation for guided assembly |
| `--read-type` | `auto` | auto/short/long | 读长类型(auto=自动检测)｜Read type (auto=auto-detect) |
| `--transcripts` | — | store_true | 额外输出transcripts.fa cDNA(需-g)｜Also output cDNA transcripts.fa (needs -g) |
| `--predict-cds` | — | store_true | TransDecoder预测CDS(需-g,输出gene/mRNA/CDS)｜TransDecoder CDS prediction (needs -g) |
| `-p, --pattern` | `*_1_clean.fq.gz` |  | FASTQ文件命名模式（*为样本名占位符）｜FASTQ file naming pattern (* is sample name placeholder) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--sample-timeout` | `43200` | int | 单个样本处理超时时间（秒）｜Sample processing timeout in seconds |
| `-s, --step` | — | 1/2/3/4/5/6/7 | 运行指定步骤｜Run only specified step (1: 索引｜index, 2: 比对｜alignment, 3: 排序｜sort, 4: 组装｜assembly, 5: 合并｜merge, 6: GFF3输出｜GFF3 output, 7: TransDecoder CDS｜TransDecoder CDS) |
| `-v, --verbose` | `0` | count | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--force` | — | store_true | 强制重新处理已完成的步骤｜Force re-process completed steps |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- HISAT2（`hisat2`、`hisat2-build`，仅 FASTQ 模式）
- SAMtools（`samtools`）
- StringTie（`stringtie`）
- gffread（GTF→GFF3、提取 cDNA）
- TransDecoder（`TransDecoder.LongOrfs`/`TransDecoder.Predict` + 两个坐标转换脚本，仅 `--predict-cds`）
- 工具路径默认按 conda 环境解析（stringtie/gffread/hisat2 在 rna 域、TransDecoder 在 annot/transdecoder 环境），可用环境变量覆盖

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。每一步按输出文件是否存在判断（索引 `.1.ht2`、SAM 非空、BAM+索引、GTF 非空、合并 GTF、GFF3、transcripts.fa、CDS 结果等）。改参数重跑请加 `--force`，否则复用旧产物。

**Q2：BAM 输入有什么要求？**
必须坐标排序（`SO:coordinate`），否则报错「BAM 非坐标排序」；缺 `.bai` 会自动补建。BAM 模式不支持步骤 1-3（`-s 1/2/3` 会报错），只支持 4-7。

**Q3：什么时候必须给 `-g` 基因组？**
FASTQ 模式、`--transcripts`、`--predict-cds`、或 `-s 7` 都需要基因组。纯 BAM 直入且只要 GFF3 结构时可不给。

**Q4：`--read-type auto` 怎么判断长读短读？**
程序采样 BAM 中约 1000 条读段，按读长中位数判断：≥ 500 bp 判为长读（走 StringTie `-L`），否则短读。FASTQ 模式恒按短读处理（HISAT2 短读比对）。

**Q5：TransDecoder 预测的 CDS 准确吗？**
本流程用基础版 TransDecoder（LongOrfs + Predict），未接入 blastp/hmmscan 同源证据，结果偏保守、可能漏掉部分真实 CDS。需要更完整注释建议用 BRAKER 等完整流程。

**Q6：单个 BAM 会合并 GTF 吗？**
不会。样本数 = 1 时跳过 StringTie merge，直接用该样本 GTF 转 GFF3，避免空合并报错。
