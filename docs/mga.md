# MGA 共识基因组组装 | MGA Consensus Genome Assembly

一句话理解：**用 HiFi reads 跑 MGA 共识组装，得到一个高质量基因组**。输入一份 HiFi reads，输出组装好的 `assembly.fasta`，内部自动完成依赖检查并记录软件版本。

## 功能概述 | Overview

- 封装 MGA（consensusLJA）共识组装，面向 HiFi reads
- 显式用 `conda run -n mga` 包装（MGA 二进制不在 conda 环境内，但依赖环境内的 minimap2/samtools/python）
- 断点续传：最终产物存在则整体跳过
- 支持 `--dry-run` 只打印命令不执行
- read name 含空格会给出警告（见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools mga -r reads.fastq.gz -o out_dir/
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HiFi reads | PacBio 高保真长读长，又长又准 |
| 共识组装(consensus) | 把多条 reads 的「共识」拼成一条可靠序列，精度高 |
| contig | 拼出的连续序列片段 |
| conda 环境 | 一套隔离的软件运行环境，装着工具及其依赖 |

## 输入 | Input

### reads 文件

HiFi reads，支持 FASTA/FASTQ，可 gzip 压缩（`.gz`）：

```text
>read_1
ACGTACGT...
>read_2
...
```

## 参数说明 | Parameters

### 必需参数与线程 | Required & threads

**通俗理解|In plain words:** `-r` 是 HiFi reads，`-o` 是输出目录，`-t` 线程数（默认 50，按机器调整）。

### 软件路径与执行控制 | Software path & execution

**通俗理解|In plain words:** `--mga-path` 是 MGA 二进制路径（默认 `~/software/MGA/consensusLJA/bin/MGA`），`--conda-env` 是运行环境名（默认 `mga`）。**这两个一般不用动**；`--dry-run` 只打印命令不真跑，适合先看命令。

## 分析流程 | Pipeline

```text
HiFi reads
    │
    ▼
检查最终产物(断点续传) → 扫描 read name(空格警告) → 依赖检查(MGA/minimap2/samtools/python)
    │
    ▼
conda run -n mga <MGA> --reads ... --output ... --threads ...
    │
    ▼
结果检查 + 写 software_versions.yml
```

## 输出 | Output

```text
out_dir/
├── 00_pipeline_info/
│   └── software_versions.yml        # 软件版本与运行参数记录
├── 99_logs/
│   └── mga.log                      # 运行日志
└── 5_polishing/
    └── assembly.fasta               # 最终组装(断点续传判断依据)
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 最终结果就是 `5_polishing/assembly.fasta`，看它的总长和序列条数是否合理。

- `5_polishing/assembly.fasta`：MGA 完成抛光(polishing)后的最终组装
- 总长应接近该物种预估基因组大小；序列条数少、连续性好为佳
- `00_pipeline_info/software_versions.yml`：记录 MGA/minimap2/samtools/python 依赖版本，写论文 Methods 时直接抄

## 参数选择建议 | Parameter Guidance

- `-t`：默认 50，按机器核数调整
- `--mga-path` / `--conda-env`：默认即可，仅在非标准安装路径时修改
- `--dry-run`：首次运行先加 `--dry-run` 看一遍实际命令，确认无误再去掉重跑

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --reads` | 必填 |  | HiFi reads(fasta/fastq,可gz)｜HiFi reads (fasta/fastq, may be gz) |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `-t, --threads` | `50` |  | 线程数｜Threads |
| `--mga-path` | — |  | MGA二进制路径｜MGA binary path |
| `--conda-env` | `mga` |  | conda环境名｜conda env name |
| `--dry-run` | `False` |  | 只打印命令不执行｜Print command without executing |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --reads` | 必填 |  | [FILE] HiFi reads(fasta/fastq,可gz)｜HiFi reads (fasta/fastq, may be gz) |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-t, --threads` | `50` | int | 线程数(默认50)｜Threads (default 50) |
| `--mga-path` | — |  | MGA二进制路径｜MGA binary path |
| `--conda-env` | `mga` |  | conda环境名(默认mga)｜conda env name (default mga) |
| `--dry-run` | — | store_true | 只打印命令不执行｜Print command without executing |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- MGA 二进制（默认路径 `~/software/MGA/consensusLJA/bin/MGA`）
- conda 环境 `mga`：需含 minimap2、samtools、以及 Python 依赖 pysam / biopython / numpy（MGA 运行所需）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。若 `5_polishing/assembly.fasta` 已存在则整体跳过。换参数重跑前需先删除该产物。

**Q2：read name 含空格会怎样？**
程序会扫描 reads 首行，若 read name 含空格会打印警告（MGA/LJA 可能因此出错）。建议先用 trim_header 类工具去掉空格。

**Q3：MGA 二进制不在 conda 环境里，能跑吗？**
能。程序显式用 `conda run -n mga --no-capture-output <MGA二进制> ...` 包装：MGA 本体在 `~/software/MGA/...`，但运行时依赖环境 `mga` 内的 minimap2/samtools/python。

**Q4：`--dry-run` 有什么用？**
只打印将要执行的命令、不真正运行，适合先确认命令和参数无误。

