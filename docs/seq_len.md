# seq-len 序列长度统计 | FASTA Sequence Length Statistics

一句话理解：**给 FASTA 文件（或整个文件夹）里的每条序列量一量长度，输出一张「序列名 + 长度」的表格，并顺手算出 N50、总长、最长/最短、平均长度等汇总指标**，用来快速评估基因组组装或转录本的连续性。

## 功能概述 | Overview

- 输入可以是单个 FASTA 文件，也可以是一个文件夹（自动识别其中的 FASTA，批量合并统计）
- 文件夹模式在结果里多一列 `source_file`，标明每条序列来自哪个文件
- 主表输出每条序列的长度，汇总表输出序列数 / 总长 / N50 / 最短 / 最长 / 平均
- 支持 `.gz` 压缩文件，流式读取，不把整条序列载入内存，大文件也能跑
- 可选按最小长度过滤、按长度降序排序
- 纯 Python 实现，无第三方软件依赖，开箱即用

## 快速开始 | Quick Start

```bash
biopytools seq-len -i genome.fa -o out.tsv
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| FASTA | 存放序列的标准文本格式，每条序列由 `>` 开头的一行（名称行）加若干行序列字母组成 |
| 序列 ID | 名称行 `>` 后面第一段文字（到第一个空格/制表符为止），是本工具统计的「名字」 |
| N50 | 把序列按长度从长到短排好，一条条累加，加到总长一半时手里那条的长度；N50 越大说明「大块头」越多，组装越连续 |
| 汇总表 | 一张只放几行统计数字的小表，回答「一共几条、总共多长、最长最短多少」 |

## 输入 | Input

- **单个文件**：标准 FASTA（扩展名 `.fa` / `.fasta` / `.fna` / `.faa` / `.ffn` / `.frn`），也支持对应的 `.gz` 压缩版。
- **文件夹**：自动挑出该目录下（仅当前层，不递归）名字像 FASTA 的文件，按文件名排序后逐个统计。
- 序列名取名称行 `>` 之后、第一个空白字符（空格/制表符）之前的那段；空序列（只有名称行没有字母）记为长度 0。

格式示例：

```text
>seq1 description here
ACGTACGTACGT
>seq2
TTTTGGGG
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这两个参数告诉工具「读什么、写到哪」，每次必填。

- `-i, --input`：输入 FASTA 文件或文件夹。
- `-o, --output`：输出位置。工具会「智能判断」——如果给的是目录（以 `/` 结尾、或已存在目录、或不像文件名），就在里面写 `<前缀>.seq_len.tsv` 和 `<前缀>.seq_len.summary.tsv`；如果给的是文件路径，就直接当主表写，汇总表自动变成 `xxx.summary.tsv`。

### 过滤与排序 | Filter & sort

**通俗理解|In plain words:** 控制「哪些序列进结果、结果按什么顺序排」。最小长度过滤像筛子，把太短的碎片挡在门外；排序只是换个展示顺序，不改变统计数字。**一般都用默认值即可。**

- `--min-length`：只保留长度不小于此值的序列（默认 0，即不过滤）。调大它可以滤掉测序/组装产生的短碎片。
- `--sort`：加上后按长度从长到短排（默认保持输入顺序）。

### 输出控制 | Output control

**通俗理解|In plain words:** 控制结果文件叫什么名字、要不要那份汇总表。**一般不用动。**

- `--prefix`：输出文件名前缀，默认取输入名（文件模式取文件名去掉扩展名，文件夹模式取目录名）。
- `--no-summary`：加上后只写主表，不写 N50 等汇总表。

### 日志 | Logging

**通俗理解|In plain words:** 控制运行日志写到哪、记多详细，只在排查问题时才需要动。

- `--log-file`：把完整日志额外写到这个文件（默认只打到屏幕）。
- `--log-level`：日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL，默认 INFO）。
- `-v, --verbose`：详细日志（等价于 DEBUG 级别）。

## 输出 | Output

假设 `-o out/`（目录模式，输入文件 `genome.fa`）：

```text
out/
├── genome.seq_len.tsv            # 主表：每条序列的长度
└── genome.seq_len.summary.tsv    # 汇总表：N50/总长等统计
```

- 主表（文件模式）两列：`sequence_id`、`length`。
- 主表（文件夹模式）三列：`source_file`、`sequence_id`、`length`。
- 汇总表列：`source_file`、`num_seqs`、`total_length`、`n50`、`min_length`、`max_length`、`mean_length`；每个文件一行，文件夹模式（超过 1 个文件）末尾再加一行 `ALL` 全局汇总。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 主表是「流水账」，汇总表是「体检报告」。看结果主要看汇总表。

- `num_seqs`：序列条数。条数多不一定好，可能是碎成了很多小段。
- `total_length`：总长度，即所有序列首尾相接的总碱基数（约等于基因组大小或转录本总规模）。
- `n50`：连续性指标，越大说明长序列占比越高。比较两个组装时，N50 大的通常更连续（但要结合总长一起看）。
- `min_length` / `max_length`：最短/最长序列，看有没有异常的超短或超长记录。
- `mean_length`：平均长度，等于总长除以条数。
- 若主表里某些序列长度为 0，说明该序列只有名称行、没有字母内容。

## 参数选择建议 | Parameter Guidance

- **`--min-length`**：评估组装连续性时可设成 500 或 1000 滤掉短碎片，看「有效序列」的 N50；纯粹统计全貌时保持默认 0。
- **`--sort`**：想直接肉眼扫「最长的是哪几条」时加上，不影响任何统计值。
- **`--prefix`**：批处理多个文件、又想统一命名时才需要，单个文件一般不用。
- **`-o` 写文件还是写目录**：想直接拿到一张表就用文件路径（如 `out.tsv`）；想同时要主表 + 汇总表且命名规范，就用目录（如 `out/`）。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | FASTA 文件或文件夹｜FASTA file or folder |
| `--output, -o` | 必填 |  | 输出 TSV 路径或目录｜Output TSV path or directory |
| `--prefix` | — |  | 输出前缀(默认取输入名)｜Output prefix (default: input name) |
| `--min-length` | `0` | int | 最小长度过滤(0=不过滤)｜Min length filter (0=no filter) |
| `--sort` | — |  | 按长度降序(默认保持输入顺序)｜Sort by length descending |
| `--no-summary` | — |  | 不输出汇总表｜Skip summary table |
| `--log-file` | — |  | 日志文件｜Log file |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--verbose, -v` | — |  | 详细日志｜Verbose |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | FASTA 文件或文件夹｜FASTA file or folder |
| `-o, --output` | 必填 |  | 输出 TSV 路径或目录｜Output TSV path or directory |
| `--prefix` | — |  | 输出前缀(默认取输入名)｜Output prefix (default: input name) |
| `--min-length` | `0` | int | 最小长度过滤(0=不过滤)｜Min length filter (0=no filter) |
| `--sort` | — | store_true | 按长度降序(默认保持输入顺序)｜Sort by length descending |
| `--no-summary` | — | store_true | 不输出汇总表｜Skip summary table |
| `--log-file` | — |  | 日志文件｜Log file |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `-v, --verbose` | — | store_true | 详细日志｜Verbose |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3 标准库（`gzip`、`logging`），无第三方软件、无 conda 环境依赖。

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。本工具是单遍流式统计，速度很快，重跑会直接覆盖同名输出，无需也不存在断点续传。

**Q2：`-o` 到底会写成文件还是目录？**
按「看起来像不像文件」判断：路径以 `/` 结尾、已存在为目录、或末尾不含 `.` 都按目录处理，否则按文件处理。拿不准就显式写成目录（末尾加 `/`）。

**Q3：序列名怎么取的？**
取 `>` 后到第一个空格/制表符之前的那段。所以 `>chr1 some comment` 只认 `chr1`；两个序列若这段相同，会在结果里重复出现（本工具不去重）。

**Q4：能处理 `.gz` 压缩文件吗？**
能，`.fa.gz`、`.fna.gz` 等 FASTA 压缩文件都会透明解压读取。

**Q5：文件夹里混了非 FASTA 文件会怎样？**
会被忽略，只统计名字以 `.fa/.fasta/.fna/.faa/.ffn/.frn`（可带 `.gz`）结尾的文件。
