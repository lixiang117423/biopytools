# RNA-seq 表达定量流程 | RNA-seq Expression Quantification Pipeline

一句话理解：**用 HISAT2 把 RNA 测序读段比到参考基因组，再用 StringTie 数每个基因有多少读段，算出 FPKM/TPM 表达量，最后合并成一张「基因 × 样本」的表达矩阵**——这是绝大多数 RNA-seq 差异表达分析的第一步。

## 功能概述 | Overview

- 标准 RNA-seq 定量流程：HISAT2 比对 + StringTie 定量
- 自动构建 HISAT2 索引（基因组 + GTF 的剪接位点/外显子）
- 逐样本计算 FPKM 与 TPM 表达量
- 合并所有样本为一个表达矩阵文件 `all_fpkm_tpm.txt`
- 支持样本目录或样本信息文件两种输入
- 断点续传：已完成的样本自动跳过；`--force` 强制重跑

## 快速开始 | Quick Start

```bash
biopytools rnaseq -g genome.fa -f genes.gtf -i ./fastq_data -o rnaseq_results
```

最小输入：参考基因组 FASTA + 基因注释 GTF + 双端 FASTQ 目录。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 比对 alignment | 把每条读段放到参考基因组上找最像的位置 |
| HISAT2 | 一个剪接感知的 RNA 比对软件，会处理内含子 |
| StringTie | 数「每个基因有多少读段」并算表达量的软件 |
| GTF 注释 | 标注基因组上基因/外显子位置的说明书 |
| FPKM / TPM | 两种表达量单位，都做了深度和长度校正，值越大表达越强 |
| 表达矩阵 | 一张表：行是基因，列是样本，格子是表达量 |
| 剪接位点 | 内含子两端的边界序列，帮助比对软件正确跨内含子 |

## 输入 | Input

### 参考文件

- `-g/--genome`：参考基因组 FASTA
- `-f/--gtf`：基因注释 GTF（用于提取剪接位点与外显子、构建索引、StringTie 定量）

### 测序数据（二选一）

1. **FASTQ 目录**：目录内放成对双端 FASTQ，按 `-p/--pattern` 命名（默认 `*_1.clean.fq.gz`，`*` 是样本名，R2 自动换为 `_2.clean.fq.gz`）
2. **样本信息文件**：制表符三列 `样本名\tR1路径\tR2路径`

```text
sampleA    /data/sampleA_1.clean.fq.gz    /data/sampleA_2.clean.fq.gz
sampleB    /data/sampleB_1.clean.fq.gz    /data/sampleB_2.clean.fq.gz
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 四项都必填：参考基因组、基因注释、测序数据、结果目录。缺一不可。

- `-g/--genome`：基因组 FASTA
- `-f/--gtf`：基因注释 GTF
- `-i/--input`：FASTQ 目录或样本信息文件
- `-o/--output`：输出目录

### 样本与文件 | Samples & files

**通俗理解|In plain words:** `-p` 是「文件名长什么样」的模板，文件命名与默认 `*_1.clean.fq.gz` 不一致时才改。`-r` 决定「跑完后删不删 BAM」——BAM 很占磁盘，只想要表达量的话可设 `yes` 删掉省空间。

- `-p/--pattern`：FASTQ 命名模式（默认 `*_1.clean.fq.gz`）
- `-r/--remove`：处理后是否删除 BAM（`yes`/`y`/`no`/`n`，默认 `no`）

### 运行与超时 | Runtime & timeout

**通俗理解|In plain words:** `-t` 是线程数。`--sample-timeout` 是「单个样本最多跑多久」——超时就跳过该样本继续下一个，防止一个坏样本卡住整个流程；默认 21600 秒（6 小时），一般够用，样本特别大时可调大。

- `-t/--threads`：线程数（默认 12）
- `--sample-timeout`：单样本超时秒数（默认 21600）

### 日志与运行模式 | Logging & mode

**通俗理解|In plain words:** `-v`/`--quiet` 控制输出多少；`--dry-run` 只打印命令不真跑（检查参数用）；`--force` 强制重跑已完成样本（换参数后想重算时用）。**这些是运维/调试选项，正常跑一次都不用动。**

- `-v/--verbose`：增加详细程度（可 `-vv`/`-vvv`）
- `--quiet`：静默模式，只输出错误
- `--log-file` / `--log-level`：日志文件路径 / 级别
- `--dry-run`：试运行，不实际执行
- `--force`：强制重新处理已完成样本

## 分析流程 | Pipeline

```text
基因组 + GTF
    |
    v
步骤1: 提取剪接位点/外显子 -> 构建 HISAT2 索引(索引建在基因组同目录)
    |
    v
步骤2: 解析样本(成对 FASTQ)
    |
    v
步骤3: 逐样本处理
   +- HISAT2 比对(小文件走管道,大文件自动拆分)
   +- StringTie -G -e 定量
   +- 提取 FPKM/TPM
   +- (可选)删除 BAM
    |
    v
步骤4: 合并所有样本 -> all_fpkm_tpm.txt
    |
    v
生成总结报告
```

## 输出 | Output

```text
output_dir/
+-- 01_bam/{sample}.sorted.bam        # 比对结果(旧版目录名 01.bam)
+-- 02_stringtie/{sample}.gtf         # StringTie 定量 GTF(旧版 02.stringtie)
+-- 03_fpkm_tpm/{sample}.fpkm.txt     # 单样本 FPKM/TPM(旧版 03.fpkm_tpm)
+-- all_fpkm_tpm.txt                  # 所有样本合并表达矩阵(主要产物)
+-- rnaseq_summary.txt                # 总结报告
+-- rnaseq_processing_*.log           # 运行日志
```

HISAT2 索引不建在输出目录内，而建在**基因组 FASTA 同目录**下：`{基因组名}.hisat2.index.*.ht2`、`{基因组名}.ss`、`{基因组名}.exon`。

## 结果解读 | Interpreting Results

### 1. 表达矩阵（`all_fpkm_tpm.txt`）

**通俗理解|In plain words:** 最核心的结果——一张长格式大表，一行 = 一个转录本在一个样本里的表达量。

```text
gene_id    transcript_id    cov    FPKM    TPM    sample
gene1      gene1.1          12.3   3.45    5.67   sampleA
gene1      gene1.1          14.1   4.02    6.10   sampleB
```

- 含 `sample` 列，是多样本拼接的长格式，下游可用 R（reshape2/tidyr）转成「基因 × 样本」宽表后做差异分析
- `FPKM`、`TPM` 越大表达越强；`cov` 是覆盖度

### 2. 单样本表达量（`03_fpkm_tpm/{sample}.fpkm.txt`）

与矩阵中该样本的内容一致，表头为 `gene_id	transcript_id	cov	FPKM	TPM	sample`。

### 3. StringTie GTF（`02_stringtie/{sample}.gtf`）

含每个转录本的结构与 `cov`/`FPKM`/`TPM` 属性，供下游组装或进一步处理。

### 4. 总结报告（`rnaseq_summary.txt`）

记录输入、参数、样本清单与输出位置，是运行留档。

## 参数选择建议 | Parameter Guidance

- **`-r/--remove yes`**：磁盘紧张、只做差异表达时建议开启，删除中间 BAM
- **`--force`**：换了 `-p` 或想重新定量时加，否则已完成样本会被跳过
- **`-t/--threads`**：样本多、核多时可调大，HISAT2 与 samtools 都能利用多线程
- **`--sample-timeout`**：个别样本特别大（如单样本超 6 小时）才需调大

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `--gtf, -f` | 必填 |  | 基因注释GTF文件路径｜Gene annotation GTF file path |
| `--input, -i` | 必填 | Path | 输入FASTQ文件目录或样本信息文件｜Input FASTQ file directory or sample information file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--pattern, -p` | `*_1.clean.fq.gz` | str | FASTQ文件命名模式｜Fastq file naming pattern (e.g., "*.R1.fastq.gz" or "*_1.fq.gz"), * represents sample name |
| `--remove, -r` | `no` | yes/y/no/n | 处理后删除BAM文件｜Remove BAM files after processing |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity (-v, -vv, -vvv) |
| `--quiet` | — |  | 静默模式｜Quiet mode, only output errors |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — |  | 试运行模式｜Dry run mode, no actual execution |
| `--force` | — |  | 强制重新处理已完成的样本｜Force re-process completed samples |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--sample-timeout` | `21600` | int | 单个样本处理超时时间（秒）｜Sample processing timeout in seconds |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组fasta文件路径｜Genome fasta file path |
| `-f, --gtf` | 必填 |  | 基因注释GTF文件路径｜Gene annotation GTF file path |
| `-i, --input` | 必填 |  | 输入fastq文件目录或样本信息文件｜Input fastq file directory or sample information file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-p, --pattern` | `*_1.clean.fq.gz` |  | Fastq文件命名模式｜Fastq file naming pattern (e.g., "*.R1.fastq.gz" or "*_1.fq.gz"), * represents sample name |
| `-r, --remove` | `no` | yes/y/no/n | 处理后删除BAM文件｜Remove BAM files after processing |
| `-t, --threads` | `8` | int | 线程数｜Number of threads |
| `--sample-timeout` | — | int | 单个样本处理超时时间（秒），默认不限制｜Sample processing timeout in seconds, default no limit |
| `-v, --verbose` | `0` | count | 增加输出详细程度 (-v, -vv, -vvv)｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — | store_true | 试运行模式，不实际执行｜Dry run mode, no actual execution |
| `--force` | — | store_true | 强制重新处理已完成的样本｜Force re-process completed samples |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- HISAT2（`hisat2`、`hisat2-build`，含 `extract_splice_sites.py`、`extract_exons.py`）
- SAMtools（`samtools`）
- StringTie（`stringtie`）
- Python 3（`pandas` 用于合并矩阵）
- 软件通过 conda 环境自动检测并包装调用

## 常见问题 | FAQ

**Q1：换参数重跑为什么结果没变？**
断点续传按「样本输出 `{sample}.fpkm.txt` 是否已存在」判断，已完成的样本会跳过。要重算请加 `--force`。

**Q2：HISAT2 索引建在哪？**
建在基因组 FASTA 的**同目录**下（不在输出目录），前缀为 `{基因组名}.hisat2.index`。同一基因组多次运行会复用该索引，不会重复构建。

**Q3：一个样本跑很久会不会卡住？**
`--sample-timeout`（默认 21600 秒）到期会跳过该样本继续下一个，流程不会整体卡死。超时的样本会在日志中标记。

**Q4：大文件和小文件的比对方式有何区别？**
程序按 FASTQ 文件大小自动切换：超过约 10 GB 走「先出 SAM 再排序」的拆分管道，更省内存；小的走「HISAT2 | samtools sort」管道。用户无需干预。

**Q5：为什么有的基因在结果里缺失？**
StringTie `-e` 模式只输出 GTF 里已注释的基因/转录本；若某基因在 GTF 中不存在或该样本无读段覆盖，就不会出现在表达矩阵里。
