# bam_cov - BAM 覆盖度统计 | BAM Coverage Statistics

一句话理解：**算出一个或多个 BAM 在指定区域（单段区间或一批 BED 区间）里每个碱基被测了多少次**（测序深度），输出一张"位置 × 样本"的深度表，并附均值/中位数等汇总。

## 功能概述 | Overview

- 基于 `samtools depth` 计算单碱基覆盖度，输出"位置 × 样本"的合并深度矩阵
- 支持两种模式：**手动单区间**（`-c` 染色体 + `-s` 起始 + 可选 `-e` 终止）或 **BED 批量区间**（`-b`）
- 支持单个 BAM 或目录（自动扫描所有 `*.bam`），多个 BAM 自动合并成一张表
- 缺索引自动用 `samtools index` 补建
- 自动生成统计摘要（均值/中位/最大/最小深度、覆盖比例等），可用 `--no-summary` 关掉
- 无断点续传（每次重算）

## 快速开始 | Quick Start

~~~bash
biopytools bam-cov -i sample.bam -c chr1 -s 1000 -e 2000
~~~

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 覆盖度(深度) | 某个碱基位置上"压着多少条读段"，深度越高说明这个位置被测得越充分 |
| BAM 索引 | BAM 的"目录"，`samtools depth` 靠它快速跳到指定区域 |
| BED | 记录一批基因组区间的文本格式（3 列：染色体、起点、终点），像一张"要查哪些位置"的清单 |
| 0-based / 1-based | 坐标起点从 0 还是 1 开始数；BED 用 0 起，samtools 用 1 起，本工具内部自动转换 |
| MAPQ / baseQ | 比对质量 / 碱基质量，低于阈值的读段/碱基会被排除 |

## 输入 | Input

- BAM 文件或目录（目录自动收集所有 `*.bam`），须为 `.bam` 格式
- 手动模式：必须给 `-c 染色体` + `-s 起始`（1-based），`-e` 不填则默认到该染色体末端
- BED 模式：3 列制表符分隔，格式如下

~~~text
chr1    1000    2000
chr1    5000    6000
chr2    100     900
~~~

BED 坐标是 0-based 半开区间（起点从 0 数、终点不含），工具内部会转成 1-based 闭区间再交给 samtools，无需手动换算。

## 参数说明 | Parameters

### 区域指定 | Region

**通俗理解|In plain words:** 二选一。要么用 `-b` 给一个 BED 批量查很多区间，要么用 `-c` + `-s` 手动给一个区间（`-e` 可选，不填查到染色体末尾）。**两者互斥**：BED 模式下不能再给 `-c`/`-s`；非 BED 模式则 `-c` 和 `-s` 都必填。

### 质量过滤 | Quality filters

**通俗理解|In plain words:** `--min-mapq`（默认 0）和 `--min-baseq`（默认 0）分别是最低比对质量和最低碱基质量。调大=只要"更可靠"的读段/碱基，深度会变低但更可信；默认 0 即不过滤。做严谨的深度评估时常见设 `--min-mapq 20`，普通看深度不用动。

### 输出与合并 | Output & merge

**通俗理解|In plain words:** `-o` 指定输出；`--no-merge` 关闭多样本合并、`--no-summary` 关闭统计摘要。多样本时默认合并成一张宽表，方便直接对比；单个样本不需要合并。这两个开关一般不用动。

## 分析流程 | Pipeline

~~~text
输入 BAM(文件/目录) + 区间(BED 或 -c/-s)
    │
    ▼
1. 检查 samtools 可用
2. 收集 BAM 列表; 缺索引则 samtools index 补建
3. 逐个 BAM: samtools depth 提取区间深度 → temp/{sample}.depth
4. 合并所有样本 → coverage.txt(位置 × 样本矩阵)
5. 生成统计摘要 → coverage_summary.txt(可关)
~~~

## 输出 | Output

~~~text
coverage.txt              # 合并后的深度矩阵(核心结果)
coverage_summary.txt      # 统计摘要(均值/中位/最大/最小/覆盖比例)
temp/                     # 中间临时文件
bam_coverage_YYYYMMDD_HHMMSS.log   # 运行日志
~~~

- **手动单区间模式** 的 `coverage.txt` 列：`Chrom  Pos  sample1  sample2 ...`（每行一个碱基位置，后接各样本深度）
- **BED 模式** 的列：`Chrom  Start  End  Pos  sample1  sample2 ...`（多了区间起止，方便对应回 BED 行）
- `coverage_summary.txt`：每个样本的 `Total_Positions / Mean_Coverage / Median_Coverage / Max / Min / Positions_>0X / >10X / >30X / Coverage_%_>0X / >10X` 等
- 若 `-o` 给的是一个已存在目录或无扩展名的路径，工具会自动在目录内生成 `coverage.txt`

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开 `coverage.txt`，找到关心的位置那一行，看对应样本列的数字就是深度；再配合 `coverage_summary.txt` 看整体深浅和覆盖是否均匀。

- 深度 = 0：该位置没有读段覆盖（可能是低测序区、重复区或过滤掉了）
- 各样本深度差异大：提示测序量或文库质量不齐
- `Coverage_%_>10X` 很低：整体测得很浅，需要更多数据
- `Mean_Coverage` 与 `Median_Coverage` 差距大：深度分布不均（有少数极深或极浅区）

## 参数选择建议 | Parameter Guidance

- **只看一个基因/区间**：手动模式 `-c chr1 -s 1000 -e 2000`
- **批量查一批候选区间**：BED 模式 `-b regions.bed`
- **要可靠深度**：加 `--min-mapq 20`；只要粗略覆盖看 0 即可
- **多 BAM 一起看**：输入目录，默认自动合并；单样本不想合并可用 `--no-merge`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入路径（BAM文件或包含BAM的目录）｜Input path (BAM file or directory containing BAM files) |
| `--bed, -b` | — | Path | BED文件路径(3列: chrom, start, end)｜BED file path (3 columns: chrom, start, end) |
| `--chromosome, -c` | — |  | 染色体名称｜Chromosome name (e.g., chr1, Chr12) |
| `--start, -s` | — | int | 起始位置｜Start position (1-based) |
| `--end, -e` | — | int | 终止位置｜End position |
| `--output, -o` | `./coverage.txt` | str | 输出文件路径｜Output file path |
| `--min-mapq` | `0` | int | 最小mapping质量｜Minimum mapping quality |
| `--min-baseq` | `0` | int | 最小碱基质量｜Minimum base quality |
| `--no-merge` | — |  | 不合并样本输出｜Do not merge sample outputs |
| `--no-summary` | — |  | 不生成统计摘要｜Do not generate summary statistics |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — |  | 试运行模式｜Dry run mode |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入路径（BAM文件或包含BAM的目录）｜Input path (BAM file or directory containing BAM files) |
| `--bed, -b` | — | str | BED文件路径(3列: chrom, start, end)｜BED file path (3 columns: chrom, start, end) |
| `-c, --chromosome` | — |  | 染色体名称｜Chromosome name (e.g., chr1, Chr12) |
| `-s, --start` | — | int | 起始位置｜Start position (1-based) |
| `-e, --end` | — | int | 终止位置｜End position |
| `-o, --output` | `./coverage.txt` |  | 输出路径：目录(已存在或无扩展名)则合并文件/temp/摘要都写入该目录；有扩展名则作为输出文件｜Output path: a directory (existing or no extension) → merged/temp/summary written inside it; with extension → output file |
| `--min-mapq` | `0` | int | 最小mapping质量｜Minimum mapping quality |
| `--min-baseq` | `0` | int | 最小碱基质量｜Minimum base quality |
| `--no-merge` | — | store_true | 不合并样本输出｜Do not merge sample outputs |
| `--no-summary` | — | store_true | 不生成统计摘要｜Do not generate summary statistics |
| `-v, --verbose` | `0` | count | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--dry-run` | — | store_true | 试运行模式｜Dry run mode |
| `-t, --threads` | `64` | int | 线程数｜Number of threads |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `samtools`（自动解析 align 域环境并经 conda run 调用；可用环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用，要求支持 `samtools depth` 子命令）

## 常见问题 | FAQ

**Q1：会不会断点续传？**
不会。每次运行都重新提取深度。BAM 索引是复用/补建的（存在则跳过），但深度计算本身不缓存。

**Q2：为什么 BED 模式和 -c/-s 不能一起用？**
两种指定区域的方式互斥，一起给会直接报参数错误。BED 是"批量清单"，`-c/-s` 是"单个区间"。

**Q3：BED 坐标要不要 +1？**
不用。BED 是 0-based 半开区间，工具内部已自动转换为 samtools 的 1-based 闭区间，你直接按 BED 惯例写即可。

**Q4：`-o` 传了一个目录会发生什么？**
若该路径是已存在目录或无扩展名，工具会在里面生成 `coverage.txt`、`coverage_summary.txt` 和 `temp/`；有扩展名则当作输出文件名（向后兼容）。

**Q5：报"未找到染色体"或深度文件为空？**
先确认 BAM 头里确实有该染色体（名字完全一致，区分大小写如 chr1 vs Chr1），再确认 BAM 索引已生成。
