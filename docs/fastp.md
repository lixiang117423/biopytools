# Fastp FASTQ 质控 | Fastp Quality Control

一句话理解：**用 fastp 对一批测序 reads 做质量控制和过滤（去接头、去低质量碱基、去短读段），并可自动配对修复双末端文件**，把原始 FASTQ 变成干净、可直接下游分析的 clean 数据。

## 功能概述 | Overview

- 批量处理一个目录（或多个样本）的 FASTQ 文件，支持双末端和单末端
- 用 fastp 做质控：质量阈值、最小长度、不合格碱基比例、N 碱基上限
- 默认只跑 fastp（可选 `--enable-pair` 追加 seqkit pair 配对修复；fastp 双端输出本身已严格成对）
- 自动检测后缀（`_1.fq.gz`/`_2.fq.gz`、`_1.fastq.gz`/`_2.fastq.gz`）、自动识别单末端
- 内置模拟数据检测（质量值全同自动把质量阈值降为 0，避免全丢）
- 支持断点续传：已完成的样本自动跳过（`--force` 强制重跑）

## 快速开始 | Quick Start

```bash
biopytools fastp -i raw_data/ -o clean_data/
```

最小输入：一个原始 FASTQ 目录（或单个文件）+ 一个输出目录。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| 质控 QC | 把 reads 里质量差、带接头、太短的部分裁掉或整条丢掉，让数据更干净 |
| 质量值 Phred / Q | 每个碱基的「可信度打分」，Q30 表示测错概率约千分之一 |
| 双末端 paired-end | 一条片段两头各测一次，R1、R2 两个文件需严格一一对应 |
| 配对修复 pair | 质控可能让某条 read 被丢、它的搭档还在，pair 步骤把这些「落单」的 read 挑出来，保证两端数量一致 |
| 接头 adapter | 测序时加在片段两端的已知短序列，太长读穿时会跑到 read 里，需裁掉 |
| 单末端 single-end | 只测片段一端，只有一份 FASTQ（如 PacBio HiFi 长读） |

## 输入 | Input

### 输入目录或文件

- 目录：扫描其中符合后缀规则的 FASTQ 文件，自动配对
- 单个文件：直接处理这一个文件；若找不到配对文件会自动按单末端处理
- 双末端默认按 `_1.fq.gz`/`_2.fq.gz` 或 `_1.fastq.gz`/`_2.fastq.gz` 自动识别，也可用 `--read1-suffix`/`--read2-suffix` 手动指定
- 目录里没有 `_1`/`_2` 配对命名时，自动扫描普通 FASTQ 并切换单末端模式

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是「原始数据在哪」（目录或单文件），`-o` 是「干净数据放哪」。两者都必填。

### 软件路径 | Software paths

**通俗理解|In plain words:** 指定 fastp 和 seqkit 的可执行文件位置。默认 `fastp` / `seqkit`（即从 PATH 找），并会自动检测它们所在的 conda 环境做包装。**一般不用动**，只有装在非标准位置时才需要给完整路径。

### 线程数 | Threads

**通俗理解|In plain words:** `-t` 控制 fastp 的 worker 线程数，默认 12。**fastp 的 worker 上限是 16**：实测给到 16 以上（如 `-t 88`）fastp 会进入"假死"状态——线程互相卡住、几乎不占 CPU、永远不出结果。所以模块会把超过 16 的值自动钳制到 16 并打印 WARNING；调 `-t` 在 1-16 之间有意义，核多也**不要超过 16**。

### 质控阈值 | QC thresholds

**通俗理解|In plain words:** 这四个数决定「什么样的 read 该被裁、被丢」。
- `-q`（质量阈值，默认 30）：低于此质量值的碱基算「不合格」
- `-u`（不合格比例，默认 40）：一条 read 里不合格碱基占比超过这个百分比，整条丢弃
- `-l`（最小长度，默认 50）：处理后短于这个长度的 read 丢弃
- `-n`（N 碱基上限，默认 10）：read 里未知碱基 N 超过这个数就丢弃

调大 `-q`、调小 `-u` 都会更严格（留得更少但更干净）；调大 `-l` 会丢更多短 read。**默认值经实践验证，一般不用动**，除非你对数据质量有明确更高/更低要求。

### 文件模式 | File patterns

**通俗理解|In plain words:** 控制「怎么认出 R1/R2、怎么抽样本名」。
- `--read1-suffix`/`--read2-suffix`：默认自动检测（支持 `_1.fq.gz`/`_2.fq.gz`、`_1.fastq.gz`/`_2.fastq.gz`）。你的命名不规范时才需要手动指定
- `--single-end`：强制单末端模式。单文件输入或目录无配对命名时会自动识别，一般无需手动加

### 配对修复 | Pair repair

**通俗理解|In plain words:** **默认只跑 fastp，不做配对检查**——fastp 处理双端数据时本来就一条 pair 一条 pair 地过滤，某条 read 被丢时它的搭档也一起丢，所以输出的 R1/R2 天然严格成对，不需要额外修复。v1.0.x 版本默认会多跑一步 seqkit pair 全量复查（把「落单」reads 分离到 `unpaired/`），它能保证配对但会把整个流程耗时翻倍，v1.1.0 起改为默认关闭。确有需求（比如怀疑上游数据本身已错位）时用 `--enable-pair` 显式开启。

### 日志选项 | Logging

**通俗理解|In plain words:** 控制终端和日志的详细程度。`-v` 更详细、`--quiet` 只报错误、`--log-level`/`--log-file` 精确指定级别或日志位置。排查问题时用 `-vv` 看 DEBUG 信息，平时默认即可。

### 执行控制 | Execution control

**通俗理解|In plain words:** `--force` 强制覆盖已存在的输出（否则已完成的样本会跳过）；`--dry-run` 只打印将要执行的命令、不真跑，用来先核对命令是否正确。

## 分析流程 | Pipeline

```text
输入目录/文件
  |
  v
1. 校验 fastp 可执行、创建输出目录
  |
  v
2. 查找样本配对(双末端/单末端)
  |
  v
3. 检测模拟数据(质量值全同 -> 自动把质量阈值降为 0)
  |
  v
4. 逐样本处理:
   双末端: fastp 质控 -> 最终 clean 文件 (--enable-pair 时多一步 seqkit pair 修复)
   单末端: fastp 质控 -> 最终 clean 文件
   (已存在最终输出则跳过,除非 --force)
  |
  v
5. 生成汇总报告 fastp_processing_summary.txt
```

## 输出 | Output

```text
clean_data/
├── sample1_1.clean.fq.gz          # 样本1的 clean R1(双末端)
├── sample1_2.clean.fq.gz          # 样本1的 clean R2(双末端)
├── sample2.clean.fq.gz            # 单末端样本的 clean 数据
├── fastp_reports/
│   ├── sample1.html               # 质控报告(浏览器打开)
│   └── sample1.json               # 质控报告(机器可读)
├── unpaired/                      # 配对修复挑出的"落单"reads(仅 --enable-pair 时)
├── fastp_processing_summary.txt   # 总览报告
└── fastp_processing.log           # 运行日志
```

（双末端样本命名 `{sample}_1.clean.fq.gz` + `{sample}_2.clean.fq.gz`；单末端命名 `{sample}.clean.fq.gz`。）

- **clean 文件**：质控后的干净数据，是下游分析的直接输入
- **fastp_reports/*.html**：fastp 的图形化质控报告，含碱基质量、GC 含量、接头含量、过滤统计等
- **fastp_reports/*.json**：同一份报告的机器可读版本
- **unpaired/**：seqkit pair 分离出的无法配对的 reads（仅 `--enable-pair` 时产生），一般不参与下游分析，可留作记录或删除
- **fastp_processing_summary.txt**：总样本数、成功/失败数、成功率、所用参数的总览

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 跑完看两处——终端日志的「成功/失败样本数」确认没出错，再打开每个样本的 `fastp_reports/*.html` 看过滤前后 reads 数变化和质量曲线。

- **成功 / 失败数**：日志末尾「成功处理 / 失败样本」应为 0 失败；有失败样本看日志里对应样本的报错
- **过滤比例**：HTML 报告里「过滤前后 reads 数」。过滤掉的比例过高（如 >50%）提示数据质量差或阈值过严
- **质量曲线**：HTML 报告里的碱基质量分布，好的数据整体质量值高、末端略有下降属正常
- **好坏判据**：失败样本为 0、R1/R2 的 clean reads 数相等（fastp 双端处理天然保证相等）、过滤比例在合理范围

## 参数选择建议 | Parameter Guidance

- **标准双末端批量质控**：默认参数即可（默认只跑 fastp，输出天然严格配对）
- **单末端长读（HiFi）**：`-i file.fq -o out`，程序自动识别单末端
- **数据质量差、想更严格**：调高 `-q`（如 30→35）、调低 `-u`（如 40→20）
- **想保留更多短 read**：调低 `-l`（如 50→30）
- **怀疑上游数据本身错位、需要强制配对**：加 `--enable-pair` 启用 seqkit pair 修复
- **先看命令不真跑**：加 `--dry-run`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入原始FASTQ数据目录或文件｜Input raw FASTQ data directory or file |
| `--output-dir, -o` | 必填 |  | 输出清洁FASTQ数据目录｜Output clean FASTQ data directory |
| `--fastp-path` | `fastp` |  | fastp可执行文件路径｜fastp executable path |
| `--threads, -t` | `12` | int | 线程数(fastp worker上限16,超出自动钳制)｜Number of threads (fastp worker threads capped at 16; excess is clamped automatically) |
| `--quality-threshold, -q` | `30` | int | 质量阈值｜Quality threshold |
| `--min-length, -l` | `50` | int | 最小长度｜Minimum length |
| `--unqualified-percent, -u` | `40` | int | 不合格碱基百分比阈值｜Unqualified base percentage threshold |
| `--n-base-limit, -n` | `10` | int | N碱基数量限制｜N base count limit |
| `--read1-suffix` | — |  | Read1文件后缀（单末端模式也使用此参数）。默认自动检测，支持_1.fq.gz和_1.fastq.gz｜Read1 file suffix (also used for single-end mode). Auto-detect by default, supports _1.fq.gz and _1.fastq.gz |
| `--read2-suffix` | — |  | Read2文件后缀。默认自动检测，支持_2.fq.gz和_2.fastq.gz｜Read2 file suffix. Auto-detect by default, supports _2.fq.gz and _2.fastq.gz |
| `--single-end` | — |  | 单末端模式（单文件输入时自动检测，无需手动指定）｜Single-end mode (auto-detected for single file input, no need to specify manually) |
| `--enable-pair` | `False` |  | 启用seqkit pair配对修复步骤（默认关闭,fastp双端输出本身已严格配对）｜Enable seqkit pair step (disabled by default; fastp PE output is already strictly paired) |
| `--disable-pair, --disable-repair` | `False` |  | 显式禁用seqkit pair配对修复步骤（与默认行为一致,向后兼容保留）｜Explicitly disable seqkit pair step (same as default, kept for backward compatibility) |
| `--seqkit-path` | `seqkit` |  | seqkit可执行文件路径｜seqkit executable path |
| `--verbose, -v` | — |  | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--force, -f` | — |  | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — |  | 模拟运行(不实际执行)｜Dry run without execution |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入原始FASTQ数据目录｜Input raw FASTQ data directory |
| `-o, --output-dir` | 必填 |  | 输出清洁FASTQ数据目录｜Output clean FASTQ data directory |
| `--fastp-path` | `fastp` |  | fastp可执行文件路径｜fastp executable path |
| `-t, --threads` | `12` | int | 线程数(fastp worker上限16,超出自动钳制)｜Number of threads (fastp worker threads capped at 16; excess is clamped automatically) |
| `-q, --quality-threshold` | `30` | int | 质量阈值｜Quality threshold |
| `-l, --min-length` | `50` | int | 最小长度｜Minimum length |
| `-u, --unqualified-percent` | `40` | int | 不合格碱基百分比阈值｜Unqualified base percentage threshold |
| `-n, --n-base-limit` | `10` | int | N碱基数量限制｜N base count limit |
| `--read1-suffix` | — |  | Read1文件后缀（单末端模式也使用此参数）。默认自动检测，支持_1.fq.gz和_1.fastq.gz｜Read1 file suffix (also used for single-end mode). Auto-detect by default, supports _1.fq.gz and _1.fastq.gz |
| `--read2-suffix` | — |  | Read2文件后缀。默认自动检测，支持_2.fq.gz和_2.fastq.gz｜Read2 file suffix. Auto-detect by default, supports _2.fq.gz and _2.fastq.gz |
| `--single-end` | — | store_true | 单末端模式（单文件输入时自动检测，无需手动指定）｜Single-end mode (auto-detected for single file input, no need to specify manually) |
| `--enable-pair` | `False` | store_true | 启用seqkit pair配对修复步骤（默认关闭,fastp双端输出本身已严格配对）｜Enable seqkit pair step (disabled by default; fastp PE output is already strictly paired) |
| `--disable-pair, --disable-repair` | `False` | store_true | 显式禁用seqkit pair配对修复步骤（与默认行为一致,向后兼容保留）｜Explicitly disable seqkit pair step (same as default, kept for backward compatibility) |
| `--seqkit-path` | `seqkit` |  | seqkit可执行文件路径｜seqkit executable path |
| `-v, --verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run without execution |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- fastp（质控核心）
- seqkit（仅启用 pair 配对修复时，默认需要）

两者通过 PATH 查找，并自动检测所在 conda 环境（`conda run -n <env> --no-capture-output` 包装）；也支持用 `--fastp-path`/`--seqkit-path` 指定路径。

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
本工具按「最终 clean 文件是否存在」做断点续传：已完成的样本会跳过。想用新参数重跑，要么加 `--force` 强制覆盖，要么删掉旧的 clean 文件和报告再跑。

**Q2：质量值全是同一种字符，fastp 把数据全丢了？**
这通常是模拟数据（如 wgsim 输出，质量值全为 `!`）。程序会自动检测并临时把质量阈值降为 0，避免全部被丢，日志里会有「检测到模拟数据」的告警。

**Q3：为什么我的文件被当成单末端处理了？**
当输入是单文件、或目录里找不到 `_1`/`_2` 配对命名时，程序会自动切单末端。若你其实是双末端但命名不规范，用 `--read1-suffix`/`--read2-suffix` 手动指定后缀即可。

**Q4：unpaired 目录里的文件是什么？**
启用 pair 时，质控会让某些 read 被丢、它的搭档还在，seqkit pair 会把这些「落单」read 挑出来放到 `unpaired/`。它们不参与下游成对分析，通常可以忽略或删除。

**Q5：R1/R2 的 read 名字带 /1、/2 或空格格式，pair 能处理吗？**
能。seqkit pair 步骤会自动检测 read 名配对格式（`/1``/2`、`.1``.2`、`_1``_2`、Casava 空格格式、无后缀等）并生成对应的 `--id-regexp`，无需手动指定。

**Q6：给了很大的 `-t`（如 88），为什么作业"假死"——不报错、几乎不占 CPU、永远没结果？**
fastp 的 worker 线程上限是 16，实测超过后（如 `-w 88`）线程互相卡死：作业显示在跑、CPU 几乎为零、输出 0 字节。v1.1.1 起模块自动把 `-t` 钳制到 16 并打印 WARNING，无需手动处理；也建议流程调用时直接给 `-t 16` 以内的值。
