# FASTQ GC 与长度过滤 | FASTQ GC & Length Filter

一句话理解：**按「GC 含量」和「序列长度」两个条件，把 FASTQ 里符合条件的 reads 留下来、不符合的丢掉**，常用于去除异常 GC 或过短的读段。

## 功能概述 | Overview

- 按 GC 含量范围（最小 / 最大）和序列长度范围（最短 / 最长）过滤 reads
- 基于 seqkit + awk 的流式实现，速度快、适合大文件
- 输入输出都支持 gzip 压缩（`.gz`）
- 输出一个过滤后的 FASTQ，并打印通过率统计

## 快速开始 | Quick Start

```bash
biopytools fastq-gc-filter -i input.fastq.gz -o filtered.fastq.gz --min-gc 25 --min-length 50
```

最小输入：一个输入 FASTQ + 一个输出文件名，其余用默认值（GC 25-100%，长度 >=50）。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| GC 含量 | 序列里 G 和 C 这两种碱基占的比例（百分比），像「一锅汤里两种主要配料的浓度」 |
| 读段长度 | 一条 read 有多少个碱基，太短的 read 信息量少、容易被误比对 |
| FASTQ | 测序读段的四行格式文本（名字 / 序列 / 加号 / 质量值） |
| seqkit | 一个功能强大的序列处理工具，本工具靠它提取信息和抽取序列 |
| awk | 一个按行处理文本的小工具，本工具靠它判断「这条 read 该不该留」 |

## 输入 | Input

### 输入 FASTQ

标准 FASTQ 文件，支持 `.gz` 压缩（`.fastq` / `.fq` / `.fastq.gz` / `.fq.gz`）。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 一个输入、一个输出。输出文件若以 `.gz` 结尾就输出压缩格式，否则输出纯文本。

### GC 含量过滤 | GC content filtering

**通俗理解|In plain words:** 决定「GC 在什么范围内的 read 才留下」。`--min-gc`（默认 25）和 `--max-gc`（默认 100）围成一个区间，GC 落在区间内的保留。调大 `--min-gc` 会丢掉更多低 GC 的 read，调小 `--max-gc` 会丢掉更多高 GC 的 read。默认 25-100 相当宽松，一般只有怀疑存在异常 GC 污染（如接头、外源序列）时才需要收紧。

### 序列长度过滤 | Sequence length filtering

**通俗理解|In plain words:** 决定「多长的 read 才留下」。`--min-length`（默认 50）是最短长度，短于它的丢掉；`--max-length`（默认不设 = 不限）是最长长度，长于它的丢掉。默认只卡最短 50，基本够用；只有想剔除异常超长 read 时才设 `--max-length`。

## 分析流程 | Pipeline

```text
输入 FASTQ
  |
  v
1. seqkit fx2tab 提取每条 read 的名字、长度、GC 含量
  |
  v
2. awk 按 GC / 长度条件筛选出符合条件的 read 名
  |
  v
3. seqkit grep 按名单抽取序列 -> 输出 FASTQ
  |
  v
4. 统计总 reads 数、通过数、通过率
```

## 输出 | Output

```text
filtered.fastq.gz     # 过滤后保留的 reads(格式随输出后缀)
```

- 单一 FASTQ 文件，只含通过 GC 和长度双重条件的 reads
- 日志打印到终端，末尾给出「总 reads 数 / 通过数 / 过滤掉数 / 通过率」

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看日志末尾的通过率：接近 100% 说明数据本身就很干净；通过率骤降说明有大量 read 落在你设的阈值之外，要么阈值太严，要么数据确实异常。

- **通过率**：日志里「通过率 Pass rate」百分比。正常数据用默认阈值通常 90% 以上
- **过滤掉数**：若异常高，先确认阈值是否设得过严（比如把 `--min-length` 设成几百），或数据是否真有大量超短 / 极端 GC 的 read
- **好坏判据**：通过率显著低于预期时，别急着改参数，先想想「这批数据是不是本来就该过滤这么多」

## 参数选择建议 | Parameter Guidance

- **常规去短读段**：默认参数即可（GC 25-100，长度 >=50）
- **去接头 / 外源污染**：把 `--min-gc` 调高（如 30-35）或 `--max-gc` 调低（如 60-70）试探
- **只要足够长的 read**：把 `--min-length` 调到 100 或更高
- **怀疑有异常超长 read**：设 `--max-length` 一个上限

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTQ文件｜Input FASTQ file (支持.gz压缩｜supports .gz compression) |
| `--output, -o` | 必填 |  | 输出FASTQ文件｜Output FASTQ file (支持.gz压缩｜supports .gz compression) |
| `--min-gc` | `25.0` | float | 最小GC含量百分比｜Minimum GC content percentage |
| `--max-gc` | `100.0` | float | 最大GC含量百分比｜Maximum GC content percentage |
| `--min-length` | `50` | int | 最短序列长度｜Minimum sequence length |
| `--max-length` | — | int | 最长序列长度｜Maximum sequence length |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTQ文件路径｜Input FASTQ file path (支持.gz压缩｜supports .gz compression) |
| `-o, --output` | 必填 |  | 输出FASTQ文件路径｜Output FASTQ file path (支持.gz压缩｜supports .gz compression) |
| `--min-gc` | `25.0` | float | 最小GC含量百分比｜Minimum GC content percentage |
| `--max-gc` | `100.0` | float | 最大GC含量百分比｜Maximum GC content percentage |
| `--min-length` | `50` | int | 最短序列长度｜Minimum sequence length |
| `--max-length` | — | int | 最长序列长度｜Maximum sequence length |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- seqkit（`seqkit fx2tab`、`seqkit grep`、`seqkit stats`）
- awk（Linux 自带）

无固定 conda 环境名；seqkit 直接通过 PATH 查找。

## 常见问题 | FAQ

**Q1：为什么 seqkit grep 用的是 24 线程，我没法改？**
本工具的 `seqkit grep` 步骤在代码里写死了 `-j 24`，没有对应参数可调。如果机器核少，可能抢资源，但目前版本无法改，需注意运行环境。

**Q2：输出想压缩还是不想压缩怎么控制？**
由输出文件名的后缀决定：以 `.gz` 结尾就压缩，否则纯文本。程序不自动加 `.gz`（与 extract-reads 不同）。

**Q3：min-gc 设成 100 以上或负数会怎样？**
程序会校验 GC 必须在 0-100、且 min 不能大于 max，否则直接报错退出。

**Q4：过滤完文件是空的？**
说明没有 read 同时满足 GC 和长度条件。检查阈值是否过于极端（比如 min-gc 比 max-gc 大、min-length 超过所有 read 的长度）。
