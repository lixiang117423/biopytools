# HiCanu 基因组组装 | HiCanu Genome Assembly

一句话理解：**把 PacBio HiFi 测序数据拼成基因组**。输入一份 HiFi reads（FASTA/FASTQ）和预估的基因组大小，输出 contigs 序列、统计信息和每条 contig 由哪些 reads 组成的映射。

## 功能概述 | Overview

- 基于 Canu（HiCanu），专为 PacBio HiFi reads 优化的组装
- 默认 `assemble` 阶段：HiFi 数据精度高，跳过纠错/修剪，直接组装（更快）
- 输出 contigs / unitigs FASTA、Canu 报告、以及 contig-reads 映射文件
- 内置组装统计（总长、contig 数、N50、L50）
- 断点续传默认启用（见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools hicanu -i reads.fastq -g 120m -p sample1 -o hicanu_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HiFi reads | PacBio 的「高保真」长读长，又长又准，是当下主流的组装原料 |
| contig | 拼装出的连续序列片段，没有空洞（gap） |
| unitig | 组装过程中「唯一路径」的片段，比 contig 更「原始」，HiFi 组装里通常用 contigs |
| 基因组大小 | 预估的基因组总长度，用 `120m`/`1g`/`5000k` 表示，像报「预算」 |
| N50 | 把 contigs 从长到短排，累加到总长一半时那条 contig 的长度；越大说明拼得越长越完整 |
| L50 | 达到 N50 所需的 contig 条数；越小越好 |
| 断点续传 | 中断后重跑，已完成的步骤自动跳过，不从头再来 |

## 输入 | Input

### reads 文件

HiFi reads，FASTA 或 FASTQ 格式（也可 gzip 压缩）：

```text
>read_1
ACGTACGT...
>read_2
...
```

- 必须同时用 `-g` 给出预估基因组大小，格式：数字 + 单位（`g`/`m`/`k`，如 `120m`、`1g`、`5000k`）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是测序数据，`-g` 是「预算」（预估基因组大小），`-p` 是给结果起的前缀名。

### 组装参数 | Assembly parameters

**通俗理解|In plain words:** `--stage` 决定「拼装的工序」。默认 `assemble` 直接拼（HiFi 数据够准，不用先纠错）；若数据质量一般，可改用 `correct`（纠错）或 `trim-assemble`（修剪+拼装）走完整流程。`--min-read-length`、`--min-overlap-length` 是过滤短 reads 和最小重叠的阈值，**默认值对 HiFi 数据已合适，一般不用动**。

### 计算资源 | Computing resources

**通俗理解|In plain words:** `-t` 线程数、`-m` 内存上限。内存按你机器的实际可用量给（组装吃内存），`--use-grid` 是提交到集群调度系统时才开。

### 执行控制 | Execution control

**通俗理解|In plain words:** `--dry-run` 只打印命令不真跑（先看看命令长什么样）；`--keep-intermediate` 保留中间文件（默认不留，省空间）。

## 分析流程 | Pipeline

```text
HiFi reads + 基因组大小
    │
    ▼
步骤1: Canu 组装(-assemble 默认, -pacbio-hifi reads)
    │
    ▼
步骤2: 复制 FASTA 到 02_fasta
    │
    ▼
步骤3: 生成 contig-reads 映射(contig_reads.tsv)
```

## 输出 | Output

```text
hicanu_output/
└── {prefix}/
    ├── 01_raw_output/                 # Canu 原始输出
    │   ├── {prefix}.contigs.fasta     # 最终 contigs(主结果)
    │   ├── {prefix}.unitigs.fasta     # unitigs(HiFi 组装常不存在,正常)
    │   ├── {prefix}.report            # Canu 报告
    │   └── {prefix}.log               # Canu 日志
    ├── 02_fasta/                      # 整理后的结果
    │   ├── {prefix}.contigs.fasta     # 复制的 contigs
    │   ├── {prefix}.unitigs.fasta     # 复制的 unitigs
    │   └── {prefix}_contig_reads.tsv  # contig→reads 映射
    ├── 03_logs/                       # 模块日志(hicanu_时间戳.log)
    └── 04_statistics/                 # 统计目录
```

- `01_raw_output/{prefix}.contigs.fasta`：最终组装结果，最常用
- `02_fasta/{prefix}_contig_reads.tsv`：两列（`contig_id`、`read_name`），记录每条 contig 由哪些 reads 组成，供下游追溯

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看 `contigs.fasta` 的「拼得好不好」——总长接近预估基因组大小、contig 条数少、N50 大，就是好组装。

- **总长 vs 基因组大小**：总长接近 `-g` 给的预估值、又不过分超出（超出可能含污染或冗余），属正常
- **N50 / L50**：N50 越大（越接近染色体长度）、L50 越小，连续性越好；HiFi 植物组装 N50 常见几 Mb 到几十 Mb
- **contig 条数**：越少越好；几百条内可接受，成千上万条说明组装不理想
- 运行日志里会打印解析出的组装大小、contig 数、N50、L50 摘要

## 参数选择建议 | Parameter Guidance

- `--stage`：HiFi 数据默认 `assemble`（最快）；数据质量差或非 HiFi 数据改用 `correct`/`trim-assemble`
- `--genome-size`：报大一点没关系（浪费些内存），报小会严重影响组装，宁大勿小
- `--min-read-length` / `--min-overlap-length`：**一般不动**，只有 reads 异常短时才调
- `--memory`：按机器实际内存给，组装吃内存，建议 >= 64G

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--reads, -i` | 必填 |  | 输入reads文件(FASTA/FASTQ)｜Input reads file path (FASTA/FASTQ format) |
| `--genome-size, -g` | 必填 |  | 基因组大小(如120m, 1g)｜Genome size (e.g., 120m, 1g) |
| `--prefix, -p` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--output-dir, -o` | `./hicanu_output` | Path | 输出目录路径｜Output directory path |
| `--canu-path` | `~/miniforge3/envs/asm/bin/canu` |  | Canu可执行文件路径｜Path to Canu executable |
| `--min-read-length` | `1000` | int | 最小read长度｜Minimum read length |
| `--min-overlap-length` | `500` | int | 最小重叠长度｜Minimum overlap length |
| `--corrected-error-rate` | — | float | 纠错后错误率｜Corrected error rate |
| `--raw-error-rate` | — | float | 原始错误率｜Raw error rate |
| `--max-input-coverage` | — | int | 最大输入覆盖度｜Maximum input coverage |
| `--stage` | `assemble` | haplotype/correct/trim/assemble/trim-assemble | 组装阶段｜Assembly stage |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory, -m` | `100G` |  | 内存限制｜Memory limit |
| `--use-grid` | — |  | 使用网格引擎｜Use grid engine |
| `--grid-options` | — |  | 网格引擎选项｜Grid engine options |
| `--dry-run` | — |  | 测试运行(不执行)｜Dry run (do not execute) |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --reads` | 必填 |  | 输入reads文件路径(FASTA/FASTQ格式)｜Input reads file path (FASTA/FASTQ format) |
| `-g, --genome-size` | 必填 |  | 基因组大小(如120m, 1g)｜Genome size (e.g., 120m, 1g) |
| `-p, --prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-o, --output-dir` | `./hicanu_output` |  | 输出目录路径｜Output directory path |
| `--canu-path` | `~/miniforge3/envs/asm/bin/canu` |  | Canu可执行文件路径｜Path to Canu executable |
| `--min-read-length` | `1000` | int | 最小reads长度｜Minimum read length |
| `--min-overlap-length` | `500` | int | 最小重叠长度｜Minimum overlap length |
| `--corrected-error-rate` | — | float | 纠错后错误率｜Corrected error rate |
| `--raw-error-rate` | — | float | 原始错误率｜Raw error rate |
| `--max-input-coverage` | — | int | 最大输入覆盖度｜Maximum input coverage |
| `--stage` | `assemble` | haplotype/correct/trim/assemble/trim-assemble | 组装阶段｜Assembly stage |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-m, --memory` | `80G` |  | 内存限制｜Memory limit |
| `--use-grid` | — | store_true | 使用网格调度｜Use grid engine |
| `--grid-options` | — |  | 网格调度选项｜Grid engine options |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run (do not execute) |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--no-resume` | — | store_true | 禁用断点续传（强制重新运行所有步骤）｜Disable resume mode (force rerun all steps) |
| `--resume` | — | store_true | 启用断点续传（默认已启用）｜Enable resume mode (enabled by default) |
| `-v, --verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Canu（默认路径 `~/miniforge3/envs/asm/bin/canu`，conda 环境 `asm`，Canu 2.3）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持，且**默认启用**。程序按 `contigs.fasta`、复制的 fasta、`contig_reads.tsv` 是否存在来判断哪些步骤已完成并跳过。注意：click 包装器未暴露 `--no-resume`/`--resume`，如需强制重跑请删除旧产物或直调模块（见参数速查「模块直调参数」）。

**Q2：换参数重跑，结果没变？**
断点续传按输出文件存在性判断。换 `--stage`、`--min-read-length` 等参数重跑前，先删除旧输出目录（或对应产物），否则会复用旧结果。

**Q3：`unitigs.fasta` 不存在正常吗？**
正常。HiFi 组装默认走 `-assemble` 阶段，通常只产出 contigs；程序已把 unitigs 视为可选文件。

**Q4：`contig_reads.tsv` 没生成？**
该文件依赖 Canu 的 `{prefix}.contigs.layout.readToTig` 和 `{prefix}.seqStore/readNames.txt` 两个中间文件；若这些中间文件缺失，映射将跳过（不影响 contigs.fasta 本身）。

