# panHiTE 群体基因组转座子分析 | panHiTE Pan-genome TE Analysis

一句话理解：用 HiTE 把多个基因组的转座子一起鉴定、聚类成一份"泛转座子库"，再比较它们在群体里的分布(哪些是核心、哪些是私有)。

## 功能概述 | Overview

- 通过 Singularity 容器调用 HiTE 的 panHiTE 流程做泛基因组转座子检测与比较
- 生成跨基因组的泛转座子库(panTE library)与各基因组 TE 注释
- 可选基因注释(GFF/GTF)与 RNA-seq 数据，用于 TE 与基因/表达的关联分析(TIDELs)
- 支持只建库(跳过分析)、断点续跑(--recover)、调试模式(--debug)
- 可指定检测的 TE 类型(LTR/TIR/Helitron/non-LTR/全部)

## 快速开始 | Quick Start

```bash
biopytools panhite -p pan_genomes/ -i genome_list.txt -o results
```

`-p` 是放所有基因组 FASTA 的目录，`-i` 是基因组列表文件。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 转座子(TE) | "跳跃基因"，能在基因组里复制并搬家的 DNA 片段 |
| 泛转座子库(panTE library) | 把多个基因组里鉴定的 TE 聚类去重后的一份代表库 |
| LTR/TIR/Helitron/non-LTR | 转座子的几大"家族"，按复制机制分类 |
| 中性突变率(miu) | 每年每个碱基位点积累多少突变，用来估算 TE 插入的"年龄" |
| TIDELs | TE 插入导致的表达差异位点，结合 RNA-seq 分析 TE 对基因表达的影响 |
| Singularity 容器 | 把整套软件+依赖打包成一个"集装箱"，免去逐个安装的麻烦 |
| 断点续跑(recover) | 从上次中断处继续，而不是从头再来 |

## 输入 | Input

### 必需

- **`-p` 泛基因组目录**：包含所有基因组文件的目录
- **`-i` 基因组列表文件**：列出要分析的基因组

基因组列表支持三种格式(见下)，由是否有基因注释/RNA-seq 决定：

```text
# 仅 TE 检测
genome1.fa
genome2.fa

# 含基因注释(第二列 GFF)
genome1.fa    annotation1.gff
genome2.fa    annotation2.gff

# 含 RNA-seq(用于 TIDELs 分析)
genome1.fa    annotation1.gff    1    rna1_1.fq.gz    rna1_2.fq.gz
genome2.fa    annotation2.gff    1    rna2_1.fq.gz    rna2_2.fq.gz
```

### 可选

- `--genes-dir`：基因注释文件目录(GFF/GTF)
- `--rna-dir`：RNA-seq 数据目录

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-p` 给"材料库"(所有基因组文件)，`-i` 给"清单"(分析哪些基因组)，`-o` 给输出目录(默认 `./panhite_output`)。

### 可选数据 | Optional data

**通俗理解|In plain words:** `--genes-dir`(基因注释)和 `--rna-dir`(RNA-seq)都是可选的。想做"TE 影响基因表达"这类分析才需要传，只做 TE 鉴定就不传。

### 容器与环境 | Container & environment

**通俗理解|In plain words:** `--singularity-cmd` 和 `--sif-file` 指向 Singularity 可执行文件和 HiTE 镜像，一般用默认值即可，只在非标准安装时改。

### 分析与恢复 | Analysis & recovery

**通俗理解|In plain words:** `--te-type` 选只检测哪类 TE(默认 all 全检)；`--miu` 是估算 TE 年龄用的突变率，一般用默认值；`--skip-analyze` 只生成 panTE 库不做后续比较；`--recover` 开启断点续跑(中途失败后加此参数从断点继续)；`--debug` 保留中间文件排错用。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--pan-genomes-dir, -p` | 必填 | Path | 泛基因组目录(包含所有基因组文件)｜Pan-genomes directory containing all genome files |
| `--input, -i` | 必填 | Path | 基因组列表文件｜Genome list file path |
| `--genes-dir` | — | Path | 基因注释文件目录(GFF/GTF)｜Genes annotation files directory (GFF/GTF) |
| `--rna-dir` | — | Path | RNA-seq数据目录｜RNA-seq data directory |
| `--output-dir, -o` | `./panhite_output` | Path | 输出目录｜Output directory path |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--singularity-cmd` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity命令路径｜Singularity executable path |
| `--sif-file` | `~/software/singularity/hite_3.3.3.sif` |  | HiTE SIF镜像路径｜HiTE SIF image path |
| `--te-type` | `all` | ltr/tir/helitron/non-ltr/all | 检测的TE类型｜TE type to detect |
| `--miu` | `1.3e-08` | float | 中性突变率(per bp per year)｜Neutral mutation rate |
| `--skip-analyze` | — |  | 跳过分析，仅生成panTE库｜Skip analysis, only generate panTE library |
| `--recover` | — |  | 是否启用断点续跑｜Whether to enable recovery mode |
| `--debug` | — |  | 是否开启调试模式｜Whether to enable debug mode |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 把输入装进 Singularity 容器跑 panHiTE，再把结果从容器复制出来并汇总。

```text
pan_genomes 目录 + genome_list 文件
    │
    ▼
步骤1: 构建 panHiTE 参数
    │
    ▼
步骤2: Singularity 容器内运行 python /HiTE/panHiTE.py
    │
    ▼
步骤3: 从容器复制结果到主机输出目录
    │
    ▼
步骤4: 扫描并整理输出文件
    │
    ▼
步骤5: 生成结果汇总(panhite_summary.json)并打印
```

## 输出 | Output

```text
panhite_output/
├── panTE_library.fa        # 泛转座子库
├── panTE_annotation.gff    # 泛转座子注释
├── *_stats.txt             # 统计信息
├── *_cluster.fa            # TE 聚类结果
├── presence_absence*.txt   # 各基因组 TE 存在/缺失表
└── panhite_summary.json    # 结果汇总(含输出文件清单与配置)
```

关键文件说明：

- **`panTE_library.fa`**：跨基因组聚类去重后的泛转座子库，是后续比较分析的核心
- **`panTE_annotation.gff`**：TE 在各基因组上的注释坐标
- **`presence_absence*.txt`**：哪些 TE 在哪些基因组存在/缺失，用于判断核心/私有 TE
- **`panhite_summary.json`**：本次运行的输出文件清单与参数快照

## 结果解读 | Interpreting Results

### 1. 泛转座子库（`panTE_library.fa`）

**通俗理解|In plain words:** 这份库是"这个群体里所有转座子的代表名单"。库的大小反映 TE 多样性。

- 库越大说明群体内 TE 种类越丰富；冗余越低说明聚类去重越有效
- 可作为其它工具(如 RepeatMasker)的输入库重复使用

### 2. 存在/缺失表（`presence_absence*.txt`）

**通俗理解|In plain words:** 像一张"签到表"，看每个 TE 出现在哪些基因组里。

- 几乎所有基因组都有的 TE = 核心/古老插入
- 只在个别基因组出现的 TE = 私有/近期插入，常与表型差异相关

### 3. 基因/TE 关联（含 RNA-seq 时）

- 若提供基因注释与 RNA-seq，可进一步分析 TE 插入对邻近基因表达的影响(TIDELs)

## 参数选择建议 | Parameter Guidance

- **`--te-type`**：只关心某一类 TE(如 LTR)可缩小范围，速度更快、结果更聚焦
- **`--skip-analyze`**：只想要 panTE 库、暂时不做比较时开启，省去后续分析时间
- **`--recover`**：长时间运行中途失败后，加此参数重跑即可续跑，不必从头再来
- **`--miu`**：默认 1.3e-8(水稻等常用值)；物种不同可查文献调整，影响 TE 年龄估算
- **`--threads`**：容器内 TE 检测并行度，按机器核数调整

## 依赖 | Dependencies

- Singularity（默认路径 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`）
- HiTE 镜像 `~/software/singularity/hite_3.3.3.sif`（容器内含 HiTE 3.3.3）
- Python 3

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持，但需要显式加 `--recover`。默认不开启，加该参数后 HiTE 会从上次中断处继续。

**Q2：找不到 Singularity 或 SIF 镜像？**
确认默认路径 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` 和 `~/software/singularity/hite_3.3.3.sif` 存在，或分别用 `--singularity-cmd`、`--sif-file` 指定。

**Q3：基因组列表的列数怎么定？**
只做 TE 检测=1 列(基因组名)；加基因注释=2 列(基因组名 + GFF)；加 RNA-seq 做 TIDELs=5 列(基因组名 + GFF + 1 + 双端 reads)。

**Q4：`--genes-dir`/`--rna-dir` 和列表里的注释列冲突吗？**
两者配合使用：`--genes-dir`/`--rna-dir` 指定这些文件所在的目录，列表里写文件名(或路径)。只做 TE 检测时不传这两个目录。

**Q5：TE 年龄是怎么算的？**
由 `--miu`(中性突变率)按 TE 内部的序列分歧度反推插入时间；`--miu` 越准确，年龄估算越可靠。
