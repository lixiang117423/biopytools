# pair_fastq - FASTQ 双端配对修复 | FASTQ Pair Fixing

一句话理解：**批量把「左右端对不上号」的双端测序文件修复好——把配得上对的 reads 挑进 paired/，配不上的单条 reads 单独放进 single/，让下游流程不再因为配对混乱报错。**

## 功能概述 | Overview { #overview }

- 扫描输入目录，按 R1/R2 后缀自动配对样本（默认 _1.fq.gz / _2.fq.gz）
- 逐样本修复配对，支持两种工具：repair.sh（bbmap，默认）与 seqkit pair
- 输出成对 reads 到 paired/，未配对的 singleton 到 single/
- 报告只有 R1 或只有 R2 的样本（告警，不处理）
- 支持 --dry-run 只打印命令不执行

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools pair-fastq -i raw_data -o fixed_data -t 16
```

最小输入：一个装双端 FASTQ 的目录（文件名用 _1.fq.gz / _2.fq.gz 区分左右端）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 双端测序(paired-end) | 一条 DNA 片段的两端各测一段，得到 R1、R2 两个文件，两条 read 本应一一对应 |
| R1 / R2 | 左端 / 右端读段；文件名里用后缀（如 _1 / _2）区分 |
| 配对(paired) | R1 和 R2 里能一一对应上的 read |
| singleton(单条 read) | 只有一端、找不到搭档的 read，修复时被单独挑出来 |
| 配对混乱 | 测序/质控后某些 read 在 R1 或 R2 缺失，导致两边数量对不上 |
| repair.sh / seqkit pair | 两个不同的修复工具，都干「重新配对」这件事 |

## 输入 | Input { #input }

- 输入目录(-i)：内含成对 FASTQ 文件，R1/R2 通过后缀识别（默认 _1.fq.gz / _2.fq.gz，可用 --suffix1 / --suffix2 改）。

示例目录：

```text
raw_data/
|-- S1_1.fq.gz
|-- S1_2.fq.gz
|-- S2_1.fq.gz
-- S2_2.fq.gz
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** -i 是输入目录，-o 是输出目录（两者都必填）。-t 线程数默认 12，调大能加快大文件处理。

### 文件配对 | Filename pairing

**通俗理解|In plain words:** --suffix1 / --suffix2 决定程序怎么从文件名认出 R1/R2 并配对。默认 _1.fq.gz / _2.fq.gz；只有你的文件名不是这个后缀时才需要改。后缀改了但两边的「样本名」必须一致（如 S1_1.fq.gz 与 S1_2.fq.gz 的 S1 是样本名）。

### 工具选择 | Tool selection

**通俗理解|In plain words:** --tool 选修复工具：repair（默认，走 bbmap 的 repair.sh）或 seqkit（走 seqkit pair）。--seqkit-bin 是 seqkit 路径，--repair-sh 是 repair.sh 脚本名，--repair-conda-env 是 repair.sh 所在 conda 环境（默认 bbmap_v.39.81），--repair-memory 是给 repair.sh 的内存（默认 300g）。一般用默认 repair 即可。

### 调试选项 | Debugging

**通俗理解|In plain words:** --dry-run 只打印将要执行的命令、不真跑（适合先看一遍）；--verbose 更啰嗦；--log-file / --log-level 控制日志输出位置和级别。一般不用动。

## 分析流程 | Pipeline { #pipeline }

```text
扫描输入目录，按后缀找出 R1/R2 并配对样本
    |
    v
对每个样本：
  repair 模式 -> conda run -n bbmap_v.39.81 repair.sh in= in2= out= out2= outs=
  seqkit 模式 -> seqkit pair -1 -2 -O paired/
    |
    v
成对 reads 写入 paired/，singleton 写入 single/
    |
    v
汇总成功/失败样本数
```

## 输出 | Output { #output }

```text
fixed_data/
|-- paired/
|   |-- S1_1.fq.gz     # 修复后的左端
|   -- S1_2.fq.gz     # 修复后的右端
|-- single/
|   -- S1_single.fq.gz # 未配对的单条 reads
-- logs/
    -- pair_fastq_YYYYMMDD_HHMMSS.log
```

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 跑完后看 paired/ 里成对文件是否齐全、single/ 里单条 reads 有多少。single/ 越大说明配对混乱越严重。

- paired/：修复后左右端一一对应的文件，文件名与输入一致，可直接进入下游流程
- single/：配不上对的单条 reads；文件大小为 0 或不存在说明该样本配对完好
- 日志：末尾有「成功 X/总数」「失败 Y/总数」汇总；失败的样本名会在日志中列出

好坏判据：paired/ 里每对 R1/R2 都存在且非空即成功；single/ 占原始数据比例越低越好（通常只占很小比例）。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规双端数据：全部默认（tool=repair、后缀 _1/_2.fq.gz）
- 后缀不同：按实际改 --suffix1 / --suffix2（如 .R1.fastq.gz / .R2.fastq.gz）
- 想先预览：加 --dry-run 看命令不执行
- 大数据量：调大 -t；repair.sh 模式内存默认 300g，一般足够

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入目录（包含FASTQ文件）｜Input directory containing FASTQ files |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--suffix1` | `_1.fq.gz` |  | R1文件后缀｜R1 file suffix |
| `--suffix2` | `_2.fq.gz` |  | R2文件后缀｜R2 file suffix |
| `--tool` | `repair` | seqkit/repair | 工具选择｜Tool selection (seqkit or repair) |
| `--seqkit-bin` | `seqkit` |  | seqkit二进制路径｜seqkit binary path |
| `--repair-sh` | `repair.sh` |  | repair.sh脚本名称｜repair.sh script name |
| `--repair-conda-env` | `bbmap_v.39.81` |  | repair.sh的conda环境名称｜conda environment for repair.sh |
| `--repair-memory` | `300g` |  | repair.sh内存参数｜repair.sh memory (-Xmx) |
| `--dry-run` | — |  | 仅显示命令不执行｜Show commands without executing |
| `--verbose` | — |  | 详细输出｜Verbose output |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入目录（包含FASTQ文件）｜Input directory containing FASTQ files |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--suffix1` | `_1.fq.gz` |  | R1文件后缀｜R1 file suffix (default: "_1.fq.gz") |
| `--suffix2` | `_2.fq.gz` |  | R2文件后缀｜R2 file suffix (default: "_2.fq.gz") |
| `--tool` | `repair` | seqkit/repair | 工具选择｜Tool selection (default: repair) |
| `--seqkit-bin` | `seqkit` |  | seqkit二进制路径｜seqkit binary path (default: "seqkit") |
| `--repair-sh` | `repair.sh` |  | repair.sh脚本名称｜repair.sh script name (default: "repair.sh") |
| `--repair-conda-env` | `bbmap_v.39.81` |  | repair.sh的conda环境名称｜conda environment name for repair.sh (default: "bbmap_v.39.81") |
| `--repair-memory` | `300g` |  | repair.sh内存参数｜repair.sh memory parameter (default: "300g") |
| `--dry-run` | — | store_true | 仅显示命令不执行｜Show commands without executing |
| `--verbose` | — | store_true | 详细输出｜Verbose output |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level (default: INFO) |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- bbmap（repair.sh，默认工具），通过 conda 环境 bbmap_v.39.81 调用
- seqkit（仅 --tool seqkit 时需要）
- conda（repair 模式用 conda run 包装）

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持。每次运行都会全量处理所有样本，已存在的输出会被覆盖。

Q2：--tool seqkit 为什么没生效？
这是当前实现的已知问题：main.py 只把 --seqkit-bin 传给了配置，--tool / --repair-sh / --repair-conda-env / --repair-memory 这几个参数虽然被 argparse 解析，但没有真正传入配置，因此始终用默认值（tool=repair、环境 bbmap_v.39.81、内存 300g）。要用 seqkit 需等代码修复。

Q3：只有 R1 没有 R2 的样本会怎样？
不会处理，只在日志里告警（提示「只有 R1 的样本」）。程序只对 R1/R2 都齐全的样本执行修复。

Q4：后缀怎么改？
用 --suffix1 / --suffix2。注意 R1/R2 文件名除后缀外的前缀（样本名）必须完全一致，否则配不上对。

Q5：repair.sh 内存不够怎么办？
--repair-memory 当前实际不生效（见 Q2），默认 -Xmx300g。若内存确实不足，需要直接改源码或降低输入数据量。