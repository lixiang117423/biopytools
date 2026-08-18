# RNA-Bloom 转录组从头组装 | RNA-Bloom De Novo Transcriptome Assembly

一句话理解：**在没有参考基因组的情况下，把一大堆 RNA 测序读段像拼图一样拼成完整的「转录本序列」**，解决「没参考、又想拿到基因转录本序列」的问题。

## 功能概述 | Overview

- 基于 RNA-Bloom 的无参考（de novo）转录组组装
- 支持多种输入：短读段双端（`--left/--right`）、短读段单端（`--sef/--ser`）、长读段（`--long`，ONT/PacBio）、单细胞混合（`--cell-list`）
- 支持长读 + 短读混合组装、链特异性数据、PacBio 特殊模式
- 可选参考转录本引导组装（`--ref`）
- 输出主转录本、短转录本、去冗余转录本三套 FASTA
- 无断点续传：重复运行会重新组装

## 快速开始 | Quick Start

```bash
biopytools rnabloom --left reads_1.fq --right reads_2.fq -o ./results
```

最小输入：双端短读段 FASTQ 各一份（或任一其他读段类型）+ 输出目录。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 从头组装 de novo | 不靠参考基因组，纯从读段自己拼出序列 |
| 转录本 transcript | 一个基因转录出来的 RNA 序列（去掉了内含子） |
| 短读段 short read | 一次读一两百碱基的测序数据（Illumina） |
| 长读段 long read | 一次读几千到几万碱基（ONT/PacBio） |
| Bloom filter | 一种省内存的「集合去重」结构，RNA-Bloom 用它高效存储海量 k-mer |
| k-mer | 从序列里滑窗取出的固定长度小片段，组装的「积木」 |
| 链特异性 stranded | 能区分读段来自 DNA 哪条链，信息更精确 |
| 去冗余 non-redundant | 去掉高度相似的重复序列，只留代表 |

## 输入 | Input

### 读段文件（至少指定一种）

| 输入方式 | 参数 | 说明 |
|---|---|---|
| 双端短读段 | `--left` + `--right` | 必须成对提供 |
| 单端短读段（正向） | `--sef` | single-end forward |
| 单端短读段（反向） | `--ser` | single-end reverse |
| 长读段 | `--long` | ONT/PacBio |
| 单细胞混合 | `--cell-list` | 列表文件（pooled assembly） |

- FASTQ 支持 `.fq`/`.fastq`，压缩或未压缩均可
- 单细胞模式（`--cell-list`）与 `--left/--right/--long` 互斥
- 参考引导组装（`--ref`）不支持长读段数据

### 参考转录本（可选）

- `--ref`：已有的参考转录本 FASTA，用于引导组装

## 参数说明 | Parameters

### 输入与输出 | Input & output

**通俗理解|In plain words:** 告诉程序「用什么数据拼」和「结果放哪」。`-o` 是唯一必填项（输入至少一种读段）。

- `--left/--right`、`--sef/--ser`、`--long`、`--cell-list`：各类读段输入
- `-o/--output-dir`：输出目录
- `--ref`：参考转录本（引导组装，长读段不可用）

### Bloom filter 内存配置 | Bloom filter

**通俗理解|In plain words:** RNA-Bloom 用 Bloom filter 存海量 k-mer，内存大小决定能装多少、以及误判率。`--mem` 是给它分配的总内存（GB），`--fpr` 是允许的假阳性率，`--nk` 是唯一 k-mer 数量的预估。**这三个一般不用手动设**——程序默认会自动估算合适的 Bloom filter 大小；只有当默认估算不准（如异常大的数据集）或想精细控制内存时，才手动指定。

- `--mem`：Bloom filter 总大小（GB）
- `--fpr`：假阳性率（0-1）
- `--nk/--num-kmers`：唯一 k-mer 数量

### 数据类型 | Data type

**通俗理解|In plain words:** 这几项描述「你的数据有什么特殊之处」。链特异性数据加 `--stranded`；若建库时左/右端 reads 方向反了，用 `--revcomp-left/--revcomp-right` 翻转；PacBio 长读段加 `--pacbio`（默认按 ONT 处理）。**只有数据确实具备这些特性时才加，普通 Illumina 双端数据一个都不用加。**

- `--stranded`：链特异性数据
- `--revcomp-left` / `--revcomp-right`：反向互补左/右端读段
- `--pacbio`：PacBio 数据（默认 ONT）

### 输出选项 | Output options

**通俗理解|In plain words:** `--min-length` 设定最短保留转录本长度，`--uracil` 让输出用 U（RNA 字母）而非 T（DNA 字母），`--no-nr` 不导出去冗余转录本。详见 FAQ：其中 `--min-length` 与 `--no-nr` 目前仅做参数校验、尚未透传到 RNA-Bloom 命令，实际过滤以 RNA-Bloom 内部默认值为准。

- `--min-length`：最小转录本长度（默认 200）
- `--uracil`：输出 U 而非 T
- `--no-nr`：不导出去冗余转录本

### 处理控制 | Processing control

**通俗理解|In plain words:** `-t` 是线程数；`--stage` 可在指定阶段后停止（调试用）；`--rnabloom-path` 指定 RNA-Bloom 路径，装好并能直接调用时不用动。

- `-t/--threads`：线程数（默认 12）
- `--stage`：停止阶段（1-3）
- `--rnabloom-path`：RNA-Bloom 工具路径（默认 `rnabloom`）

## 分析流程 | Pipeline

```text
检查依赖(Java 11/17、minimap2、ntCard 可选)
    |
    v
定位 RNA-Bloom(命令或 JAR)
    |
    v
构建 RNA-Bloom 参数(输入/内存/数据类型/输出)
    |
    v
执行组装(内部含 k-mer 建库 -> 图构建 -> 延伸 -> 转录本输出)
    |
    v
验证输出 + 统计转录本数量
```

## 输出 | Output

```text
output_dir/
+-- rnabloom.transcripts.fa         # 主转录本 FASTA(主要产物)
+-- rnabloom.transcripts.short.fa   # 短转录本
+-- rnabloom.transcripts.nr.fa      # 去冗余转录本
+-- rnabloom_assembly.log           # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 主转录本（`rnabloom.transcripts.fa`）

**通俗理解|In plain words:** 组装出来的「最终成果」——一条条转录本序列，每条以 `>` 开头的名称 + 序列。转录本数量、总长度、N50 是评价组装质量的核心指标（可用 `biopytools assembly-stats` 或 `busco` 进一步评估）。

### 2. 短转录本与去冗余转录本

- `rnabloom.transcripts.short.fa`：较短、可能不够完整的转录本，单独列出便于排查
- `rnabloom.transcripts.nr.fa`：去冗余后的代表转录本，若下游做注释/定量，常用这套以降低冗余

### 3. 组装日志（`rnabloom_assembly.log`）

运行结束后日志会打印各输出文件及其转录本数量（用 `grep -c '^>'` 统计），若某文件缺失或为 0 条，说明组装未成功，需看日志中的报错。

## 参数选择建议 | Parameter Guidance

- **短读段**：Illumina 双端数据直接 `--left/--right`，其余参数全默认
- **长读段**：ONT 数据用 `--long`；PacBio 数据再加 `--pacbio`；可 `--long + --sef` 混合组装提升完整性
- **单细胞**：`--cell-list` 走 pooled 组装，注意它与其他输入互斥
- **链特异性**：确认建库是 stranded 才加 `--stranded`，否则可能组装出方向错误的序列
- **内存**：默认自动估算即可；超大转录组或报内存不足时，再用 `--mem` 手动调大

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--left, -1` | — |  | 左端reads文件｜Left reads file (paired-end) |
| `--right, -2` | — |  | 右端reads文件｜Right reads file (paired-end) |
| `--sef` | — |  | 单端正向reads文件｜Single-end forward reads file |
| `--ser` | — |  | 单端反向reads文件｜Single-end reverse reads file |
| `--long` | — |  | 长reads文件｜Long reads file (ONT/PacBio) |
| `--cell-list` | — |  | 单细胞列表文件｜Single-cell list file (pooled assembly) |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--rnabloom-path` | `rnabloom` |  | RNA-Bloom工具路径｜RNA-Bloom tool path |
| `--mem` | — | float | Bloom filter总大小(GB)｜Total Bloom filter size in GB |
| `--fpr` | — | float | 假阳性率｜False positive rate (0-1) |
| `--nk, --num-kmers` | — | int | 唯一k-mer数量｜Number of unique kmers |
| `--stranded` | — |  | 链特异性数据｜Strand-specific data |
| `--revcomp-left` | — |  | 反向互补左端reads｜Reverse-complement left reads |
| `--revcomp-right` | — |  | 反向互补右端reads｜Reverse-complement right reads |
| `--pacbio` | — |  | PacBio数据｜PacBio data (default: ONT) |
| `--ref, --reference` | — |  | 参考转录本文件｜Reference transcript file |
| `--min-length` | `200` | int | 最小转录本长度｜Minimum transcript length |
| `--uracil` | — |  | 输出尿嘧啶(U)而非胸腺嘧啶(T)｜Write uracil (U) instead of thymine (T) |
| `--no-nr` | — |  | 不导出去冗余转录本｜Do not export non-redundant transcripts |
| `--stage` | — | 1/2/3 | 停止阶段｜Stop at stage (1-3) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--left, -1` | — |  | 左端reads文件｜Left reads file (paired-end) |
| `--right, -2` | — |  | 右端reads文件｜Right reads file (paired-end) |
| `--sef` | — |  | 单端正向reads文件｜Single-end forward reads file |
| `--ser` | — |  | 单端反向reads文件｜Single-end reverse reads file |
| `--long` | — |  | 长reads文件｜Long reads file (ONT/PacBio) |
| `--cell-list` | — |  | 单细胞列表文件｜Single-cell list file (pooled assembly) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--rnabloom-path` | `rnabloom` |  | RNA-Bloom工具路径｜RNA-Bloom tool path (default: rnabloom) |
| `--mem` | — | float | Bloom filter总大小(GB)｜Total Bloom filter size in GB |
| `--fpr` | — | float | 假阳性率｜False positive rate (0-1) |
| `--nk, --num-kmers` | — | int | 唯一k-mer数量｜Number of unique kmers |
| `--stranded` | — | store_true | 链特异性数据｜Strand-specific data |
| `--revcomp-left` | — | store_true | 反向互补左端reads｜Reverse-complement left reads |
| `--revcomp-right` | — | store_true | 反向互补右端reads｜Reverse-complement right reads |
| `--pacbio` | — | store_true | PacBio数据｜PacBio data (default: ONT) |
| `--ref, --reference` | — |  | 参考转录本文件｜Reference transcript file |
| `--min-length` | `200` | int | 最小转录本长度｜Minimum transcript length (default: 200) |
| `--uracil` | — | store_true | 输出尿嘧啶(U)而非胸腺嘧啶(T)｜Write uracil (U) instead of thymine (T) |
| `--no-nr` | — | store_false | 不导出去冗余转录本｜Do not export non-redundant transcripts |
| `--stage` | — | 1/2/3 | 停止阶段｜Stop at stage (1-3) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Java 11 或 17（运行 RNA-Bloom JAR 必需）
- RNA-Bloom（经 conda run 自动检测包装，可用 --rnabloom-path 或环境变量 RNABLOOM_PATH 覆盖；未安装于 conda 环境时回退 PATH 直接调用）
- minimap2（长读段组装时必需；自动解析 align 域环境并经 conda run 调用，可用环境变量 MINIMAP2_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- ntCard（可选，用于 k-mer 计数；无对应功能域环境，可用环境变量 NTCARD_PATH 覆盖）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。每次运行都会删除旧日志、从头执行组装；中断后重跑会重新开始。RNA-Bloom 本身也不提供续传。

**Q2：`--min-length` 和 `--no-nr` 设了为什么没效果？**
这两个参数目前只做了校验，尚未透传到实际的 RNA-Bloom 命令行（`build_command_args` 未包含它们）。实际的最短长度与去冗余由 RNA-Bloom 内部默认决定。如确有需求，请在下游用 `seqkit` 等工具自行过滤。

**Q3：为什么要 Java？**
RNA-Bloom 是 Java 程序。请确保 `java -version` 能输出 11 或 17；若报「Java 未安装」，先装 JDK 再重试。

**Q4：`--ref` 引导组装有什么限制？**
`--ref`（参考转录本引导）不支持长读段数据（长读段与参考转录本互斥）。短读段有近缘物种的参考转录本时，用它可提升组装质量。

**Q5：单细胞模式能和普通读段混用吗？**
不能。`--cell-list` 与 `--left/--right/--long` 互斥，单细胞 pooled 组装需单独运行。
