# 装配质量 QV 计算 | Assembly Quality QV Calculation

一句话理解：**用量化打分告诉你「这次组装错得有多离谱」**——通过 k-mer 光谱分析算出组装 QV 值，QV 越高、错误越少，是衡量基因组组装准确度的核心指标。

## 功能概述 | Overview { #overview }

- 封装 Merqury（meryl + qv.sh），用 k-mer 光谱对比 reads 与组装，算出 QV 值
- 支持 Illumina 双端 / 单端、HiFi / PacBio / ONT 等数据类型（`--data-type auto` 可自动识别）
- 自动识别并配对 R1/R2 文件；k-mer 大小不指定时自动选 21
- 输出 `reads.meryl`（k-mer 数据库）和 `qv_result.qv`（QV 结果，已加表头）
- 步骤简单：建 meryl 数据库 → 跑 qv.sh → 解析 QV，数据库已存在时自动跳过

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools assembly-qv -i fastq_dir/ -g genome.fa
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| QV 值 | 组装「准确度」的打分；QV40 约等于每 1 万个碱基错 1 个，数值越大越准 |
| 错误率（error rate） | QV 的反面：多少比例碱基可能是错的 |
| k-mer | 把序列切成固定长度 k 的小片段，像把一句话拆成若干 k 个字的小纸条 |
| k-mer 光谱 | 统计每种小纸条出现的次数；reads 与组装的「纸条谱」对不上，就说明组装有错 |
| meryl | 高效统计 k-mer 的工具，Merqury 用它建数据库 |

## 输入 | Input { #input }

- `-i` 输入：FASTQ 文件或目录，支持 `.fastq` / `.fq`（可带 `.gz` 及大写扩展名）
- `-g` 输入：基因组 FASTA，支持 `.fa` / `.fasta`（可带 `.gz`）

配对识别规则：文件名含 `_R1` / `.R1` / `_read1` 视为 R1，对应 R2 为样本名配对的 R2 文件；双端优先取第一对，单端 / HiFi 最多取前 5 个文件参与建库。

## 参数说明 | Parameters { #parameters }

### 数据类型 | Data type

**通俗理解|In plain words:** 告诉程序你的 reads 是哪一类，以便选合适的建库方式。默认 `auto` 会按文件名特征猜（含 hifi/ccs/subreads 判为 HiFi，含 R1/R2 判为 Illumina）。**绝大多数情况保持 `auto` 即可，识别错了再手动指定。**

相关参数：`--data-type`（auto / illumina / hifi）。

### k-mer 与线程 | K-mer and threads

**通俗理解|In plain words:** `-k` 是 k-mer 长度，越短越灵敏、越长越特异，程序默认自动选 21，一般不用动（手动指定须在 15–31 之间）；`-t` 是线程数，只影响建 meryl 数据库的速度，调大更快但吃内存。**默认值足够日常使用。**

相关参数：`-k, --kmer-size`、`-t, --threads`。

### 软件路径 | Software path

**通俗理解|In plain words:** `--conda-env` 指定 Merqury 所在的 conda 环境路径，默认 `~/miniforge3/envs/merqury_v.1.3/bin/`。**只要按默认路径装好 Merqury 就不用动它**，装到别处才需要改。

相关参数：`--conda-env`。

## 分析流程 | Pipeline { #pipeline }

```text
输入 FASTQ + 基因组 FASTA
    │
    ▼
步骤1: 构建 meryl k-mer 数据库 (reads.meryl, 已存在则跳过)
    │
    ▼
步骤2: 运行 merqury qv.sh 计算 QV
    │
    ▼
步骤3: 解析并输出 qv_result.qv (自动加表头)
```

## 输出 | Output { #output }

```text
输出目录/
├── reads.meryl/      # meryl k-mer 数据库
├── qv_result.qv      # QV 结果（带表头）
└── merqury_qv.log    # 运行日志
```

`qv_result.qv` 列含义：序列名、N50/contig 数、总碱基数、QV 值、错误率。

## 结果解读 | Interpreting Results { #interpreting-results }

- `QV值` 列是核心结论：≥40 优秀，≥30 良好，20–30 可用，<20 组装质量较差（对应错误率约 >1%）
- `错误率` 列与 QV 互补：QV 每降低约 3.32，错误率翻 10 倍（QV40 ≈ 0.01% 错误）
- `总碱基数` 反映参与评估的组装规模，过小（如远小于基因组）说明评估覆盖不全

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规评估：`-i reads_dir -g genome.fa`，其余全默认
- HiFi 数据文件名不含 hifi/ccs 字样时：加 `--data-type hifi`
- 想复现论文参数：`-k 21`（HiFi 常用）；大基因组提速：`-t 24`
- 换参数重跑前先删 `reads.meryl`，否则会复用旧数据库

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | FASTQ文件或目录｜FASTQ file or directory |
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./assembly_qv_output` |  | 输出目录｜Output directory |
| `-k, --kmer-size` | — | int | K-mer大小｜K-mer size (default: auto-select) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--conda-env` | `~/miniforge3/envs/merqury_v.1.3/bin/` |  | Conda环境路径｜Conda environment path |
| `--data-type` | `auto` | auto/illumina/hifi | 数据类型｜Data type |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | FASTQ文件或目录｜FASTQ file or directory |
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./merqury_qv_output` |  | 输出目录｜Output directory |
| `-k, --kmer-size` | — | int | K-mer大小｜K-mer size |
| `-t, --threads` | `24` | int | 线程数｜Number of threads |
| `--conda-env` | `~/miniforge3/envs/merqury_v.1.3/bin/` |  | Conda环境路径｜Conda environment path |
| `--data-type` | `auto` | auto/illumina/hifi | 数据类型｜Data type |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Merqury v1.3（含 meryl、qv.sh），conda 环境名 / 路径 `merqury_v.1.3`（默认 `~/miniforge3/envs/merqury_v.1.3`）
- 需能定位 Merqury 安装：默认从 conda 环境的 `share/merqury` 找 `merqury.sh`，或设置 `MERQURY` 环境变量
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：报「无法找到 Merqury 安装」？**
程序按 `--conda-env` 的 `../share/merqury` 找 `merqury.sh`。若装到别处，设 `MERQURY` 环境变量指向 merqury 目录，或核对 `--conda-env` 路径。

**Q2：不指定 `-o` 时结果到底写到哪？**
命令行层显示的默认是 `./assembly_qv_output`，但模块内部实际默认是 `./merqury_qv_output`。为避免混淆，建议显式用 `-o` 指定输出目录。

**Q3：为什么只用了部分 FASTQ 文件？**
双端数据只取第一对 R1/R2，单端 / HiFi 最多取前 5 个文件，是为避免命令过长。多样本建议逐个样本目录分别运行。

**Q4：会断点续传吗？**
部分支持：`reads.meryl` 已存在时跳过建库（换 k 值重跑前需删掉它）；QV 计算步骤不跳过，会重新运行。
