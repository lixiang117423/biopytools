# trimal - 多序列比对修剪 | trimal MSA Trimming

一句话理解：**对多序列比对做「瘦身」——把 gap 太多、太保守的列去掉，留下真正有信息量的列，让后续建树更准更快。**

## 功能概述 | Overview

- 6 种修剪方法：`automated1`(默认,自动)、`gappyout`、`strict`、`strictplus`、`gt`(按 gap 阈值)、`cons`(按保守度阈值)
- 输出格式可沿用输入(keep)或转成 fasta / phylip / nexus / clustal 等 7 种
- 支持输出「新旧列号映射」、被删列的「互补比对」、蛋白到核酸的「反向翻译」
- 主修剪步骤支持断点续传(已有结果则跳过)

## 快速开始 | Quick Start

```bash
biopytools trimal -i aln.fasta -o out/ --method automated1
```

最小输入：一个多序列比对文件(`-i`)。结果写到 `out/01_trimal/`。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| 多序列比对(MSA) | 把多条同源序列「对齐排好」的表格，每列是同一「位点」 |
| gap(空位) | 某条序列在这个位置「缺了一块」，比对里用 `-` 表示 |
| 列(column) | 比对里竖着的一列，即所有序列在同一个位点的状态 |
| 修剪(trimming) | 删掉「质量差」的列(gap 多、太保守)，只留「有信息量」的列 |
| 保守度 | 一列里所有序列「长得有多像」；全一样=完全保守(无信息量)，有差异=有信息量 |
| 反向翻译(backtranslation) | 拿蛋白质比对去对应回核酸(CDS)序列的「倒推」操作 |

## 输入 | Input

### 多序列比对文件

trimAl 支持 FASTA、PHYLIP、NEXUS 等多种比对格式。FASTA 示例：

```text
>seq1
ACGT--ACGTACGT
>seq2
ACGT--ACGTACGT
>seq3
ACGTAAACGTACGT
```

比对必须「对齐」(各序列等长、含 gap)，这是所有多序列比对的共同前提。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 只有 `-i` 是必填(输入比对)。`-o/--output-dir` 默认 `./trimal_output`。

### 修剪方法 | Trimming method

**通俗理解|In plain words:** `--method` 决定「按什么规则删列」，是核心参数：
- `automated1`(默认)：自动平衡「删 gap」和「留信息」，**最适合做最大似然(ML)建树前的预处理**
- `gappyout`：只按 gap 多少自动删，更保守
- `strict` / `strictplus`：更激进的删除，删得更狠
- `gt`：手动指定「一列里 gap 最多允许占多少比例」(`--gt-threshold`，默认 0.9)
- `cons`：手动指定「一列至少要多保守才保留」(`--cons-threshold`，默认 80)

`--gt-threshold` 和 `--cons-threshold` **只在对应方法下才生效**，用其它方法时忽略。

### 输出格式 | Output format

**通俗理解|In plain words:** `--format` 决定输出文件格式，默认 `keep` 即「输入是什么格式、输出就是什么格式」，**一般不用改**。需要下游工具吃特定格式时再换(如 `phylip` 给 RAxML)。

### 附加输出 | Extra outputs

**通俗理解|In plain words:** 这些是可选的「加分项」，默认都关：
- `--colnumbering`：输出一张「新旧列号对照表」，告诉你保留的列在原比对里是第几列
- `--complementary`：额外输出「被删掉的列」拼成的互补比对(反着留)
- `--backtrans`：输入是蛋白质比对时，配合 CDS 文件把结果「倒推」回核酸
- `--keep-header`：保留完整 FASTA 头信息(默认可能截断)

### 输出与日志 | Output and logging

**通俗理解|In plain words:** `--sample-name` 是输出文件名前缀(默认用输入文件名)。`--log-file` 自定义日志路径(默认写到 `99_logs/trimal.log`)。`-v/--verbose` 打印更详细日志。

## 输出 | Output

```text
out/
├── 01_trimal/
│   ├── <sample>.trimmed.fasta        # 修剪后的比对(主结果)
│   ├── <sample>.complementary.fasta  # 互补比对(仅 --complementary)
│   ├── <sample>.backtrans.fasta      # 反向翻译比对(仅 --backtrans)
│   └── <sample>.colnumbering.tsv     # 新旧列号映射(仅 --colnumbering)
├── 00_pipeline_info/
│   └── software_versions.yml         # 软件版本与参数记录
└── 99_logs/
    └── trimal.log                    # 运行日志
```

其中 `<sample>` 是 `--sample-name`(默认输入文件 basename)，扩展名随 `--format` 变化。

## 结果解读 | Interpreting Results

### 1. 修剪后比对(trimmed)

**通俗理解|In plain words:** 主结果，删掉「坏列」后的比对。它可以直接喂给 RAxML、IQ-TREE 等建树软件。通常修剪后序列会变短，但留下的列「更有信息量」。

### 2. 新旧列号映射(colnumbering.tsv)

**通俗理解|In plain words:** 一张对照表：修剪后比对的第 1 列对应原比对第几列。当你需要「把结果里的位点位置映射回原始坐标」时用它。

### 3. 互补比对(complementary)

**通俗理解|In plain words:** 和主结果「互补」——主结果留的列它删掉，主结果删的列它留下。用来检查「删掉的到底是不是垃圾」很有用。

## 参数选择建议 | Parameter Guidance

- **建树前常规预处理**：默认 `--method automated1`，这是专为 ML 建树调好的自动方案
- **只想去掉明显的 gap 区**：用 `gappyout`(更保守，删得少)
- **想删得更狠、更快**：用 `strict` 或 `strictplus`
- **要精确控制**：`gt` 配 `--gt-threshold`(gap 容忍度)或 `cons` 配 `--cons-threshold`(保守度)
- **蛋白比对要转回核酸**：加 `--backtrans <CDS文件>`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入比对文件｜Input alignment file |
| `--output-dir, -o` | `./trimal_output` | Path | 输出目录｜Output directory |
| `--method` | `automated1` | automated1/gappyout/strict/strictplus/gt/cons | 修剪方法｜Trimming method (automated1 适合 ML 建树｜suited for ML trees) |
| `--gt-threshold` | `0.9` | float | gap 阈值[0,1],仅 method=gt｜gap threshold [0,1], only method=gt |
| `--cons-threshold` | `80` | int | 保守度阈值[0,100],仅 method=cons｜conservation [0,100], only method=cons |
| `--format` | `keep` | keep/fasta/phylip/phylip_paml/clustal/nexus/nbrf/mega | 输出格式(keep=沿用输入)｜Output format (keep=input format) |
| `--colnumbering` | — |  | 输出新旧列号映射｜Output old/new column mapping |
| `--backtrans` | — | Path | CDS 文件,AA→NT 反向翻译｜CDS file for AA→NT backtranslation |
| `--complementary` | — |  | 输出被修剪列的互补比对｜Output complementary alignment of trimmed columns |
| `--keep-header` | — |  | 保留完整 FASTA 头｜Keep full FASTA headers |
| `--sample-name` | — |  | 输出文件名前缀(默认输入 basename)｜Output filename prefix (default: input basename) |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--verbose, -v` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入比对文件｜Input alignment file |
| `-o, --output-dir` | `./trimal_output` |  | 输出目录｜Output directory (default: ./trimal_output) |
| `--method` | `automated1` |  | 修剪方法｜Trimming method (default: automated1) |
| `--gt-threshold` | `0.9` | float | gap 阈值[0,1],仅 method=gt｜gap threshold [0,1], only method=gt (default: 0.9) |
| `--cons-threshold` | `80` | int | 保守度阈值[0,100],仅 method=cons｜conservation [0,100], only method=cons (default: 80) |
| `--format` | `keep` |  | 输出格式｜Output format (default: keep=沿用输入｜input format) |
| `--colnumbering` | — | store_true | 输出新旧列号映射｜Output old/new column mapping |
| `--backtrans` | — |  | CDS 文件,AA→NT 反向翻译｜CDS file for AA→NT backtranslation |
| `--complementary` | — | store_true | 输出被修剪列的互补比对｜Output complementary alignment of trimmed columns |
| `--keep-header` | — | store_true | 保留完整 FASTA 头｜Keep full FASTA headers |
| `--sample-name` | — |  | 输出文件名前缀｜Output filename prefix (default: 输入 basename｜input basename) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 软件<br>Software | 说明<br>Description |
|------|------|
| trimal | trimAl 独立编译的 C++ 二进制，默认路径 `~/miniforge3/envs/phylo/bin/trimal`，可用环境变量 `TRIMAL_PATH` 覆盖 |

trimAl 是独立二进制，直接按绝对路径调用(无需 conda 环境激活，也不经过 `conda run` 包装)。

## 常见问题 | FAQ

### 1. 会断点续传吗？

**主修剪会**：如果目标文件 `<sample>.trimmed.<ext>` 已存在，会跳过主修剪步骤。但 `--complementary` 和 `--backtrans` 两个附加输出**不参与断点续传**，每次都会重新生成。

### 2. gt 方法的 gap 阈值参数名是什么？

是 `-gt <值>`(本工具已按 trimAl v1.5 适配)。注意不是旧版文档里常见的 `-gapthreshold`。用 `--method gt --gt-threshold 0.8` 即可。

### 3. 反向翻译(--backtrans)是干什么的？

当你手上是**蛋白质比对**，但最终想要**核酸比对**建树时用：提供对应的 CDS 核酸序列文件，trimAl 会把蛋白比对「倒推」回核酸比对。

### 4. 输入格式必须是 FASTA 吗？

不是。trimAl 支持 FASTA、PHYLIP、NEXUS、CLUSTAL 等多种格式，默认 `--format keep` 会沿用输入格式输出。
