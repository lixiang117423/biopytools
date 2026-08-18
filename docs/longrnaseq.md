# 三代转录组比对 | Long RNA-seq Alignment

一句话理解：**用 minimap2 把三代长读段（ONT/PacBio 的 RNA 数据）比到参考基因组上，产出排序好的 BAM 文件与比对统计**，是长读段转录组组装、定量、可变剪接分析之前的公共第一步。

## 功能概述 | Overview

- 面向三代长读段（Oxford Nanopore ONT / PacBio）的剪接感知比对（splice-aware）
- 输入灵活：BAM 或 FASTQ（单端/双端）均可，也支持整个目录批量处理
- 使用 minimap2 的 `-x splice` 模式，自动处理内含子（intron）
- 自动建参考基因组索引（缺 `.fai` 时用 samtools faidx 补）、比对、排序、建索引、统计一气呵成
- 输出排序 BAM + 比对率/覆盖度统计 + 汇总报告
- 无断点续传：重复运行会从头重新比对

## 快速开始 | Quick Start

```bash
biopytools longrnaseq -i input.bam -r genome.fa -o output_dir
```

最小输入：一个 BAM 或 FASTQ 文件 + 一个参考基因组 FASTA。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 长读段 long read | 一次能读几千到几万碱基的测序数据（ONT/PacBio），像「长文章」而非碎片 |
| 剪接感知比对 | 比对时允许读段跨过内含子，把两端分别对到外显子上 |
| 内含子 intron | 基因里会被剪掉、不进成熟 RNA 的那段 |
| 参考基因组 | 该物种「标准答案」的完整 DNA 序列 |
| BAM | 比对结果的标准二进制格式，记录每条读段落在基因组哪里 |
| MAPQ | 比对质量分，越高越可信 |
| 次优比对 secondary | 一条读段除最佳位置外，可能还有的备选位置 |
| 比对率 mapping rate | 比对上的读段占总读段的比例，衡量数据与参考是否匹配 |

## 输入 | Input

### 输入文件（三选一）

1. **BAM 文件**（`.bam`）：先内部转成 FASTQ 再比对（注意：BAM 转 FASTQ 会丢失双端配对信息）
2. **FASTQ 文件**：单端（`.fq`/`.fastq`）或双端（`.fq.gz`/`.fastq.gz`，按 `_R1/_R2`、`_1/_2` 等命名自动配对）
3. **目录**：批量处理目录内所有 BAM/FASTQ 文件（自动排除 R2 文件，只保留 R1）

### 参考基因组

- FASTA 格式；若同目录缺 `.fai` 索引会自动用 samtools 生成

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 是测序数据，`-r` 是参考基因组，`-o` 是结果目录。`-s` 给样本起名——不填就自动从文件名里猜（去掉 `.clean`/`.trimmed`/`R1` 等标记）。**一般不用手动指定 `-s`。**

- `-i/--input-file`：输入 BAM/FASTQ 文件或目录
- `-r/--ref-genome`：参考基因组 FASTA
- `-o/--output-dir`：输出目录
- `-s/--sample-name`：样本名（可选，默认从文件名提取）

### 比对参数 | Alignment

**通俗理解|In plain words:** `--max-intron` 是「允许跨过的内含子最长多长」——真核基因有的内含子很长，默认 100000 bp 覆盖绝大多数情况，一般不用动。`--no-secondary` 关掉次优比对输出（只保留最佳位置），结果更干净但会丢掉多位置比对的读段。

- `--max-intron`：最大内含子长度（默认 100000 bp）
- `--no-secondary`：不输出次优比对（默认会输出）
- `--min-mapq`：最小比对质量（默认 20，范围 0-60；见 FAQ 说明）

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t/--threads` 是线程数，长读段比对较吃 CPU，机器核多可调大。`--minimap2-path`、`--samtools-path` 是软件路径，默认自动解析功能域环境，缺失时回退 PATH，**装好并进 PATH 就不用动**。

- `-t/--threads`：线程数（默认 12）
- `--minimap2-path`：minimap2 路径（默认自动解析 align 功能域环境，缺失时回退 PATH 查找）
- `--samtools-path`：samtools 路径（默认 `samtools`）

## 分析流程 | Pipeline

```text
输入(BAM 或 FASTQ 或目录)
    |
    v
步骤1: 准备 FASTQ(BAM 则先 samtools fastq 转换)
    |
    v
步骤2: minimap2 剪接比对(-x splice) -> SAM
    |
    v
步骤3: SAM -> BAM -> 排序 -> 建索引
    |
    v
步骤4: 生成统计(flagstat/stats/depth)
    |
    v
步骤5: 质量检查(quickcheck + 比对率提示)
    |
    v
输出: 排序 BAM + 统计 + 汇总报告
```

## 输出 | Output

```text
output_dir/
+-- {sample}.sorted.bam            # 排序后的比对结果
+-- {sample}.sorted.bam.bai        # BAM 索引
+-- logs/alignment_{sample}.log    # 运行日志
+-- stats/{sample}_stats.txt       # flagstat + 平均深度
+-- stats/{sample}_detail_stats.txt# samtools stats 详细统计
+-- stats/{sample}_summary.txt     # 样本汇总(比对率等)
+-- alignment_summary.txt          # 整体汇总报告
+-- tmp/                           # 中间文件(运行后可清理)
```

## 结果解读 | Interpreting Results

### 1. 排序 BAM（`{sample}.sorted.bam`）

**通俗理解|In plain words:** 这是最核心产物——每条读段「落在基因组哪个位置」的完整记录，可直接喂给下游 StringTie 组装、featureCounts 定量、rMATS 剪接分析等。

### 2. 比对统计（`stats/{sample}_stats.txt`）

**通俗理解|In plain words:** 一份体检报告，回答「有多少读段比上去了」。

- 关键指标 `mapped` 一行：比对率。**长读段 RNA 比对率一般在 70% 以上算正常**；若低于 70%，程序会打印警告，需检查数据质量或参考基因组是否选对
- 平均深度（`Average depth`）：覆盖基因组的平均层数

### 3. 汇总报告（`alignment_summary.txt`）

记录本次运行的参考基因组、样本名、输入类型、比对参数、输出文件位置，并附下游分析建议。

## 参数选择建议 | Parameter Guidance

- **`--max-intron`**：默认 100000 覆盖绝大多数真核内含子；分析植物/动物复杂基因组一般不用改，若已知物种有超长内含子可适当调大
- **`--no-secondary`**：做定量/组装前建议加，让结果只保留最佳比对；研究多拷贝基因或转座子时才需要保留次优比对
- **`-t/--threads`**：minimap2 对多线程利用好，核多（如 32/64）可调大显著加速
- **`-s/--sample-name`**：批量目录模式下会自动用文件夹/文件名，一般无需指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-file` | 必填 |  | 输入文件或文件夹（BAM/FASTQ）｜Input file or directory (BAM/FASTQ) |
| `-r, --ref-genome` | 必填 |  | 参考基因组文件｜Reference genome file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-s, --sample-name` | — |  | 样本名称(可选，默认从输入文件名提取)｜Sample name (optional, auto-extracted from input filename) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--max-intron` | `100000` | int | 最大intron长度｜Maximum intron length |
| `--min-mapq` | `20` | int | 最小mapping quality｜Minimum mapping quality |
| `--no-secondary` | `False` |  | 不输出次优比对｜Do not output secondary alignments |
| `--minimap2-path` | `minimap2` |  | minimap2可执行文件路径｜minimap2 executable path |
| `--samtools-path` | `samtools` |  | samtools可执行文件路径｜samtools executable path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-file` | 必填 |  | 输入文件或文件夹（BAM/FASTQ）｜Input file or directory (BAM/FASTQ) |
| `-r, --ref-genome` | 必填 |  | 参考基因组文件｜Reference genome file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-s, --sample-name` | — |  | 样本名称 (可选，默认从输入文件名提取)｜Sample name (optional, auto-extracted from input filename) |
| `-t, --threads` | `64` | int | 线程数 (默认: 64)｜Number of threads (default: 64) |
| `--max-intron` | `100000` | int | 最大intron长度，默认: 100000｜Maximum intron length, default: 100000 |
| `--min-mapq` | `20` | int | 最小mapping quality，默认: 20｜Minimum mapping quality, default: 20 |
| `--no-secondary` | — | store_true | 不输出次优比对｜Do not output secondary alignments |
| `--minimap2-path` | `minimap2` |  | minimap2可执行文件路径 (默认: minimap2)｜minimap2 executable path (default: minimap2) |
| `--samtools-path` | `samtools` |  | samtools可执行文件路径 (默认: samtools)｜samtools executable path (default: samtools) |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- minimap2（自动解析 align 域环境并经 conda run 调用，可用 `--minimap2-path` 或环境变量 MINIMAP2_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- samtools（自动解析 align 域环境并经 conda run 调用，可用 `--samtools-path` 或环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。每次运行都会删除旧日志、从头重新比对。中断后重跑会重新执行全部步骤，不跳过已完成步骤。

**Q2：`--min-mapq` 设了但好像没起作用？**
本模块的 minimap2 比对命令实际用的是分数阈值（`-M 0.9` 最小比对分数、`-p 0.1` 次优分数阈值），`--min-mapq` 目前只做范围校验（0-60）并记录到日志/报告，并未作为硬过滤应用到输出 BAM。如需按 MAPQ 过滤，请比对后用 `samtools view -q` 自行过滤。

**Q3：输入 BAM 会怎样处理？**
内部先用 `samtools fastq` 转成 FASTQ 再比对，此过程会丢失双端配对信息（转成单端）。若你有原始 FASTQ，建议直接喂 FASTQ 以保留配对信息。

**Q4：目录模式怎么配对双端？**
目录模式下只挑 R1 文件（自动排除 `_R2.`/`_2.`/`-R2.`/`-2.` 命名），每个 R1 再自动找配对 R2。单个文件模式下也会按常见命名自动找 R2。

**Q5：比对率偏低怎么办？**
先确认参考基因组是正确物种/版本；再检查数据是否被污染、接头是否已去除。长读段 RNA 数据本身错误率较高，比对率 70%-90% 属常见范围。
