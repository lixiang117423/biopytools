# TGS-GapCloser Gap 填充 | Gap Filling with TGS-GapCloser

**使用三代长读长（ONT/PacBio/Hifi）数据填充基因组组装中的 N-gap | Fill N-gaps in genome assemblies using TGS long reads**

## 功能概述 | Overview

`gap-fill` 是对 TGS-GapCloser（及其 v2 版本）的封装，使用三代长读长测序数据（ONT、PacBio CLR 或 PacBio HiFi）对 scaffold 级别的基因组组装中的 gap（N 区段）进行填充。流程会调用 minimap2 将长读段比对到 scaffold 上，找到跨 gap 的 reads，并基于这些 reads 的共有序列填充 gap，显著提升组装连续性和完整度。

本模块支持三种 TGS 类型（`ont`/`pb`/`hifi`，对应不同 minimap2 参数和默认同一性阈值），并可选启用 Racon（长读段纠错）或 Pilon（NGS 短读段精修）模式进一步提升填充质量。对于难以一次填满的 gap，工具内置了第二轮 quartet 填充流程（基于 hifiasm unitig 或额外 contig 文件），可显著降低残留 gap 数量。典型场景：染色体级别组装前的 gap 填充、组装质量提升（提升 N50、降低 N 比例）、与 Hi-C 挂载流程衔接。

## 快速开始 | Quick Start

```bash
# 基本用法：ONT reads 填充 scaffold 中的 gap
biopytools gap-fill -s scaffolds.fa -t ont -ir ont_reads.fq.gz -o gapclosed

# HiFi reads，启用 Racon 纠错，64 线程
biopytools gap-fill -s asm.fa -t hifi -ir hifi.fastq.gz -o out \
    -m racon --racon-path /path/to/racon -threads 64
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-s, --scaff-file` | 输入 scaffold 文件（FASTA）|
| `-t, --tgstype` | TGS 类型：`ont` / `pb` / `hifi` |
| `-ir, --reads-file` | 输入 TGS reads 文件（FASTA/FASTQ，可 .gz）|
| `-o, --output-prefix` | 输出前缀 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --mode` | `none` | 纠错模式：`none` / `racon` / `pilon` |
| `--tgsgapcloser-path` | 自动检测 | TGS-GapCloser 可执行文件路径 |
| `-threads` | `12` | 线程数 |
| `-chunk` | `3` | 分块数量（并行切分输入）|
| `-idy, --min-idy` | 自动 | 最小同一性（按 TGS 类型自动设置）|
| `-l, --min-match` | 自动 | 最小匹配长度（按 TGS 类型自动设置）|
| `-min-nread` | `1` | 填充一个 gap 所需的最少 reads 数 |
| `-max-nread` | `-1` | 最多使用 reads 数（-1 为不限）|
| `-max-candidate` | `200` | 每个 gap 的最大候选 reads 数 |
| `-g-check` | 关 | 启用 gap 大小差异检查 |
| `-racon, --racon-path` | 无 | Racon 可执行文件路径（`-m racon` 时使用）|
| `-racon-round` | `3` | Racon 纠错轮数 |
| `-pilon, --pilon-path` | 无 | Pilon 可执行文件路径（`-m pilon` 时使用）|
| `-ngs, --ngs-file` | 无 | NGS reads 文件（`-m pilon` 时使用）|
| `-java, --java-path` | 无 | Java 可执行文件路径（Pilon 依赖）|
| `-samtools, --samtools-path` | 无 | samtools 可执行文件路径 |
| `-pilon-mem` | `300G` | Pilon 内存 |
| `-pilon-round` | `3` | Pilon 纠错轮数 |
| `-minmap-arg` | 无 | 自定义 minimap2 参数（覆盖默认）|
| `-ug, --unitig-file` | 无 | hifiasm unitig/contig 文件（第二轮填充）|
| `-fl, --flanking-len` | `5000` | 第二轮填充的 flanking 序列长度（bp）|
| `-al, --min-align-len` | `1000` | 第二轮最小比对长度（bp）|
| `-ai, --min-identity` | `40` | 第二轮最小比对同一性（%）|
| `-mfl, --max-filling-len` | `1000000` | 最大填充长度（bp）|

（运行 `biopytools gap-fill -h` 查看完整参数列表）

### TGS 类型默认参数 | TGS Type Defaults

不同类型自动应用不同 minimap2 预设和同一性阈值：ONT 使用 `-x map-ont`、PacBio CLR 使用 `-x map-pb`、HiFi 使用 `-x map-hifi`。第二轮 quartet 填充对 HiFi 使用 `-x asm20`，其它使用 `-x asm5`。

## 输出 | Output

```
./
├── {prefix}.gapcloser.fa          # 填充后的最终基因组（FASTA）
├── {prefix}.scaff_seqs            # TGS-GapCloser v2 原始输出
├── {prefix}.round1                # 第 1 轮结果备份（进入第 2 轮时）
├── tgsgapcloser.log               # 运行日志
└── done_step*                     # 各步骤完成标记文件
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --scaff-file` | 必填 |  | 输入scaffold文件｜Input scaffold file |
| `-t, --tgstype` | 必填 | ont/pb/hifi | TGS类型｜TGS type (ont/pb/hifi) |
| `-ir, --reads-file` | 必填 |  | 输入TGS reads文件｜Input TGS reads file |
| `-o, --output-prefix` | 必填 |  | 输出前缀｜Output prefix |
| `-m, --mode` | `none` | none/racon/pilon | 纠错模式｜Error correction mode (default: none) |
| `--tgsgapcloser-path` | — |  | TGS-GapCloser路径｜TGS-GapCloser path (default: auto-detect) |
| `-idy, --min-idy` | — | float | 最小同一性｜Min identity (auto-set) |
| `-l, --min-match` | — | int | 最小匹配长度｜Min match length (auto-set) |
| `-threads, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-chunk` | `3` | int | 分块数量｜Chunk count (default: 3) |
| `-g-check` | — |  | 启用Gap大小差异检查｜Enable gap size difference check |
| `-min-nread` | `1` | int | 最小reads数量｜Min read count (default: 1) |
| `-max-nread` | `-1` | int | 最大reads数量｜Max read count (default: -1) |
| `-max-candidate` | `200` | int | 最大候选数｜Max candidates (default: 200) |
| `-racon, --racon-path` | — |  | Racon路径｜Racon path |
| `-racon-round` | `3` | int | Racon轮数｜Racon rounds (default: 3) |
| `-pilon, --pilon-path` | — |  | Pilon路径｜Pilon path |
| `-ngs, --ngs-file` | — |  | NGS reads文件｜NGS reads file |
| `-java, --java-path` | — |  | Java路径｜Java path |
| `-samtools, --samtools-path` | — |  | Samtools路径｜Samtools path |
| `-pilon-mem` | `300G` |  | Pilon内存｜Pilon memory (default: 300G) |
| `-pilon-round` | `3` | int | Pilon轮数｜Pilon rounds (default: 3) |
| `-minmap-arg` | — |  | 自定义minimap2参数｜Custom minimap2 arguments |
| `-ug, --unitig-file` | — |  | hifiasm unitig/contig文件（第2轮填充）｜hifiasm unitig/contig file (round 2) |
| `-fl, --flanking-len` | `5000` | int | Flanking序列长度（bp）｜Flanking sequence length (bp) (default: 5000) |
| `-al, --min-align-len` | `1000` | int | 最小比对长度（bp）｜Min alignment length (bp) (default: 1000) |
| `-ai, --min-identity` | `40` | int | 最小比对同一性（%）｜Min alignment identity (%) (default: 40) |
| `-mfl, --max-filling-len` | `1000000` | int | 最大填充长度（bp）｜Max filling length (bp) (default: 1000000) |
| `-mgl, --min-gap-length` | `100` | int | 第2轮识别/填充的最小gap长度(bp)｜Min gap length (bp) for round-2 (default: 100) |
| `-f, --force` | — |  | 忽略断点续传强制重跑｜Force rerun, ignore checkpoint |
| `--dry-run` | — |  | 只打印命令不执行｜Dry run, print commands only |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- TGS-GapCloser / TGS-GapCloser2（核心填充工具，自动检测路径或用 `--tgsgapcloser-path`）
- minimap2（比对，TGS-GapCloser 内部调用）
- （可选）Racon（`-m racon` 模式）
- （可选）Pilon + Java + samtools + NGS reads（`-m pilon` 模式）

可设置环境变量 `TGSGAPCLOSER_PATH` 指定 TGS-GapCloser 路径。

## 引用 | Citation

- Xu, M. et al. TGS-GapCloser: a fast and accurate gap closing tool for scaffolding next-generation and third-generation genome assemblies. *Quantitative Biology* 18, 1-12 (2020).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [TGS-GapCloser 仓库](https://github.com/BGI-Qingdao/TGS-GapCloser)
