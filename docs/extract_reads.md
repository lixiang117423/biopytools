# 提取 Reads | Extract Reads

一句话理解：**给你一张「contig 对应哪些 reads」的对应关系表，从一个大 FASTQ 文件里把这些指定的 reads 挑出来**，打包成一个新 FASTQ，常用于只保留组装结果里某几条序列用到的原始读段。

## 功能概述 | Overview

- 依据 contig-reads 对应关系表（TSV）从 FASTQ 中提取指定 reads
- 输入 FASTQ 支持 gzip 压缩，输出默认压缩（可关闭）
- 纯 Python 实现，无外部软件依赖，逐行扫描、内存占用低
- 命中即写、全部命中可提前退出，大批量 reads 也能跑

## 快速开始 | Quick Start

```bash
biopytools extract-reads -m contig_reads.tsv -i input.fq.gz -o output.fq.gz
```

最小输入：一个对应关系表 + 一个 FASTQ + 一个输出文件名。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| contig | 组装后得到的一段连续序列，像「拼好的长片段」 |
| read | 测序仪读出的短片段，通常几十到几百个碱基，FASTQ 里一行条记录就是一条 read |
| 对应关系表 mapping | 记录「哪个 contig 由哪些 read 拼成」的两列表格 |
| gzip | 一种压缩格式，`.gz` 后缀文件是压缩过的文本，本工具读写都支持 |
| FASTQ 四行格式 | 每条 read 固定四行：名字、序列、加号、质量值 |

## 输入 | Input

### 对应关系表（TSV）

两列、制表符分隔：第一列 contig 名，第二列 read 名；可选一行 `#` 开头的表头（有则跳过）。示例：

```text
#contig_id	read_name
contig1	read_001
contig1	read_002
contig2	read_003
```

### 输入 FASTQ

标准 FASTQ 文件，扩展名须为 `.fq` / `.fq.gz` / `.fastq` / `.fastq.gz` 之一。read 名按「`@` 后、第一个空格前」的部分去和对应关系表匹配（如 `@read_001 1:N:0:...` 取 `read_001`）。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 三样缺一不可：对应关系表告诉程序「要挑哪些 read」，输入 FASTQ 是「从哪挑」，输出文件是「挑出来放哪」。输出若是 `.gz` 结尾就压缩，不是则默认也会自动补上 `.gz`。

### 是否压缩 | Compression

**通俗理解|In plain words:** `--no-compress` 关闭压缩输出。默认压缩（省空间、下游更常用）；只有你的下游工具非要纯文本 FASTQ 时才加这个开关。

## 分析流程 | Pipeline

```text
对应关系表 + FASTQ
  |
  v
1. 解析对应关系表，收集所有要提取的 read 名
  |
  v
2. 逐条扫描 FASTQ，命中目标 read 就暂存
  |
  v
3. 写入输出 FASTQ(默认 gzip 压缩) -> output.fq.gz
```

## 输出 | Output

```text
output.fq.gz          # 提取出的 reads(默认压缩;--no-compress 则为纯文本 .fq/.fastq)
```

- 单一 FASTQ 文件，内容为对应关系表中列出的、且在输入 FASTQ 中实际找到的 reads
- 文件名若未以 `.gz` 结尾，压缩模式下会自动补 `.gz`
- 日志打印到终端（stdout），不单独写日志文件

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 跑完后看日志里的「找到目标reads / 目标reads总数」比例：接近 100% 说明对应关系表基本都命中；明显偏低说明很多 read 名对不上（表里的名字和 FASTQ 里的名字格式不一致）。

- **命中比例**：日志中「找到目标reads X / 目标reads总数 Y」。X / Y 越接近 1 越好
- **完全没命中**：程序会报「未找到任何目标reads」并失败，通常是 read 名格式不一致（比如表里带 `/1` 后缀、FASTQ 里不带）
- **好坏判据**：提取出的 reads 数应约等于对应关系表中去重后的 read 数；对不上先核对两边 read 名格式

## 参数选择建议 | Parameter Guidance

- **默认全默认**：直接给三个文件路径即可，压缩交给默认行为
- **下游要纯文本**：加 `--no-compress`
- **只想要某条 contig 的 reads**：把对应关系表预先过滤成只含那条 contig 的行再喂进来

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--mapping, -m` | 必填 |  | contig-reads对应关系文件(TSV格式)｜contig-reads mapping file (TSV format) |
| `--input, -i` | 必填 |  | 输入FASTQ文件｜Input FASTQ file |
| `--output, -o` | 必填 |  | 输出文件｜Output file |
| `--no-compress` | — |  | 不压缩输出文件｜Do not compress output files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-m, --mapping` | 必填 |  | contig-reads对应关系文件(TSV格式)｜contig-reads mapping file (TSV format) |
| `-i, --input` | 必填 |  | 输入FASTQ文件(支持gzip压缩)｜Input FASTQ file (gzip supported) |
| `-o, --output` | 必填 |  | 输出文件｜Output file |
| `--no-compress` | — | store_true | 不压缩输出文件｜Do not compress output files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- 仅 Python 标准库（gzip），无任何外部软件、无 conda 环境名

## 常见问题 | FAQ

**Q1：为什么输出文件名多了一个 .gz？**
压缩模式（默认）下，若你给的输出名不是 `.gz` 结尾，程序会自动补 `.gz`。想保持原样用 `--no-compress`。

**Q2：对应关系表第一行要不要写表头？**
可选。若第一行以 `#` 开头会被当作表头跳过；否则从第一行就按数据解析。别把真正的数据行写成 `#` 开头。

**Q3：read 名怎么才算匹配上？**
程序取 FASTQ 头 `@` 后、第一个空格前的部分（如 `read_001`），与对应关系表第二列逐字比对。两边命名规则必须完全一致，否则匹配不上。

**Q4：输入 FASTQ 扩展名报错？**
输入只认 `.fq` / `.fq.gz` / `.fastq` / `.fastq.gz`。其他后缀（如 `.txt`）会被拒，改个名或转格式即可。
