# 基因分组核苷酸多样性 | Nucleotide Diversity per Gene Group (pi4gene)

一句话理解：**把一组序列按分组（基因家族/同源基因/等位基因）提出来、做多序列比对、算出每个分组的核苷酸多样性 π**，用来比较不同基因家族之间的变异程度。

## 功能概述 | Overview

- 按分组提取序列 → MAFFT 多序列比对 → 计算 π（Nei & Li, 1979）
- 每个分组输出一个 π 值和序列条数
- 断点续传：已提取/已比对的分组自动跳过
- 分组内序列数 <2 时自动跳过（无法计算 π）

## 快速开始 | Quick Start

```bash
biopytools pi4gene -i genes.fasta -d groups.txt -o pi4gene_output
```

最小输入：一个序列 FASTA + 一个分组 ID 文件（两列：分组名、序列 ID）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| π（核苷酸多样性） | 组内序列两两之间平均有多少个碱基不一样；越大=这组序列变异越大 |
| 多序列比对(MSA) | 把多条序列「对齐」，让同源位置排到同一列，才能公平地比差异 |
| gap（缺口） | 比对时补进去的空位；计算 π 时全 gap 位点跳过、部分 gap 只算有数据的部分 |
| 分组 ID 文件 | 一张表告诉程序「哪条序列属于哪一组」 |

## 输入 | Input

### 序列 FASTA

标准 FASTA，序列 ID 须与分组 ID 文件第二列一致：

```text
>geneA_sample1
ATGCCGTAA
>geneA_sample2
ATGTCGTAA
```

### 分组 ID 文件

每行两列「分组名 序列ID」，分隔符可自动识别（TAB / 逗号 / 空格均可）：

```text
geneA    geneA_sample1
geneA    geneA_sample2
geneB    geneB_sample1
```

- 序列 ID 在 FASTA 里不存在的会告警并跳过；某分组所有 ID 都找不到时该分组被跳过
- 每个分组至少 2 条序列才能算 π

## 分析流程 | Pipeline

```text
输入 FASTA + 分组 ID 文件
    │
    ▼
步骤1: 按分组提取序列 → 01_mafft/{group}.fasta
    │
    ▼
步骤2: MAFFT 多序列比对 → 01_mafft/{group}.aligned.fasta
    │
    ▼
步骤3: 计算 π(Nei & Li 1979) → pi_results.tsv
```

## 输出 | Output

```text
pi4gene_output/
├── pi_results.tsv                    # π 结果汇总(核心)
├── 00_pipeline_info/
│   └── software_versions.yml         # 软件版本与参数
├── 01_mafft/
│   ├── {group}.fasta                 # 各分组提取的序列
│   └── {group}.aligned.fasta         # 各分组 MAFFT 比对结果
└── 99_logs/
    └── pi4gene_analysis.log          # 运行日志
```

## 结果解读 | Interpreting Results

- **`pi_results.tsv`（核心表）**：三列 `Group / Pi / N_seq`，每个分组一行
- **π 越大 = 该分组序列之间变异越大**；比较不同基因家族时，π 高的家族进化更快/多样性更高
- **`N_seq`（序列条数）**：太少（如仅 2 条）时 π 估计代表性有限，解读需谨慎
- **`.aligned.fasta`**：可直接用别的软件（如建树）继续分析，是可靠的中间产物

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 这个工具参数极少，基本只有线程数和 MAFFT 路径可能需要动。

- **`-t/--threads`**：默认 12；序列特别多或每条很长时加大可加速 MAFFT
- **`--mafft-path`**：MAFFT 不在默认 conda 环境（phylo）时，用这个指定路径
- 其余无需调整

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 序列FASTA文件路径｜Input sequence FASTA file path |
| `-d, --id-file` | 必填 |  | 分组ID文件路径（第一列分组，第二列序列ID）｜Group ID file path (col1: group, col2: seq_id) |
| `-o, --output-dir` | `./pi4gene_output` | Path | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--mafft-path` | — |  | MAFFT路径｜MAFFT path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 序列FASTA文件路径｜Input sequence FASTA file path |
| `-d, --id-file` | 必填 |  | 分组ID文件路径（第一列分组，第二列序列ID）｜Group ID file path (col1: group, col2: seq_id) |
| `-o, --output-dir` | `./pi4gene_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--mafft-path` | — |  | MAFFT路径｜MAFFT path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- MAFFT（默认 conda 环境 phylo 的 `~/miniforge3/envs/phylo/bin/mafft`，可用 `MAFFT_PATH` 环境变量或 `--mafft-path` 覆盖）
- Python 库：biopython（读 FASTA、算 π）

## 常见问题 | FAQ

**Q1：某些分组没出现在结果里？**
两类情况会被跳过：分组内序列数 <2（无法算 π），或该分组所有序列 ID 在 FASTA 里都找不到（日志会告警）。请核对分组 ID 文件与 FASTA 的序列 ID 是否完全一致。

**Q2：换序列/分组重跑，结果没变？**
断点续传按 `01_mafft/{group}.fasta`（提取）和 `{group}.aligned.fasta`（比对）是否存在判断。改了输入后需删除 `01_mafft/` 和 `pi_results.tsv`，否则会复用旧结果。

**Q3：π 的范围是多少？**
对单条比对，π 取值在 0（完全相同）到约 0.75（每对序列都完全不同，DNA 有 4 种碱基）之间；实际基因通常远小于 0.1。
