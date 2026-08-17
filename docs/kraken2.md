# Kraken2 - 宏基因组物种分类 | Kraken2 Metagenomic Classification

一句话理解：**给一批宏基因组测序 reads「报户口」——判断每条 read 来自哪个物种，再用 Bracken 把条数换算成各物种的相对丰度占比**。

## 功能概述 | Overview { #overview }

- 自动扫描目录识别双端 FASTQ 样本（默认后缀 _1.clean.fq.gz / _2.clean.fq.gz）
- 逐样本运行 Kraken2（配对、gzip 输入、--use-names 输出学名）
- 可选 Bracken 丰度重估，把「分类条数」纠正为「物种相对丰度」
- 汇总所有样本为一张分类总表与一张物种丰度矩阵
- 按样本断点续传：某样本输出已存在即跳过

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools kraken2 -i ./fastq/ -d ~/database/kraken2_db -o ./kraken2_output
```

最小输入：一个装双端 FASTQ 的目录 + 一个建好的 Kraken2 数据库目录。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 宏基因组(metagenomics) | 不挑物种、直接把环境里所有微生物的 DNA 混在一起测序 |
| k-mer | 把 DNA 切成固定长度的小片段，像给每条序列拍「指纹」 |
| Kraken2 | 用小片段指纹去「物种词典」里查，判断每条 read 最像谁 |
| 分类(classification) | 把每条 read 归到某个物种（或「未分类」） |
| 置信度(confidence) | 有多少指纹指向同一物种，比例越高越确定 |
| Bracken | 把 Kraken2 的「条数」按物种基因组大小/命中概率纠正成「真实占比」 |
| 分类级别 | D域/P门/C纲/O目/F科/G属/S种/S1株，越往下越细 |
| read 长度 | 单端读长（bp），Bracken 估算时需要知道 |

## 输入 | Input { #input }

- **FASTQ 目录**（`-i`）：内含成对的双端文件，R1/R2 通过后缀配对（默认 `_1.clean.fq.gz` / `_2.clean.fq.gz`，可用 `--r1-suffix` / `--r2-suffix` 改）。
- **Kraken2 数据库目录**（`-d`）：已用 kraken2-build 建好（含 hash.k2d、opts.k2d、taxo.k2d 等）。

示例目录：

```text
fastq/
├── S1_1.clean.fq.gz
├── S1_2.clean.fq.gz
├── S2_1.clean.fq.gz
└── S2_2.clean.fq.gz
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** `-i` 输入 FASTQ 目录、`-d` 数据库目录、`-o` 输出目录（默认 ./kraken2_output）。`-t` 线程数默认 12。

### Kraken2 分类 | Kraken2 classification

**通俗理解|In plain words:** `--confidence` 是「指纹一致性」门槛，默认 0.0（不做过滤，全部归类）；想只保留高置信结果可设如 0.2。**默认 0 最省心，追求严格再调高。**

### Bracken 丰度重估 | Bracken abundance

**通俗理解|In plain words:** `--read-len` 必须和你的实际读长一致（默认 150），Bracken 才估算得准；`--bracken-level` 决定输出到哪个分类级别（默认 S=种，S1=株更细）；`--bracken-threshold` 是过滤低丰度的最小读数（默认 10）。**只要种级结果用默认即可**；`--no-bracken` 可跳过 Bracken 只看 Kraken2 原始分类。

### 样品命名 | Sample naming

**通俗理解|In plain words:** `--r1-suffix` / `--r2-suffix` 决定程序怎么从文件名认出 R1/R2 并配对。**只有文件名不是默认后缀时才需要改。**

## 分析流程 | Pipeline { #pipeline }

```text
扫描输入目录 → 按后缀配对样本
    │
    ▼
逐样本: Kraken2 分类(--paired --gzip-compressed --report --use-names)
    │
    ▼
逐样本: Bracken 丰度重估(可选,--no-bracken 跳过)
    │
    ▼
汇总: kraken2_summary.tsv + bracken_species.tsv
```

## 输出 | Output { #output }

```text
kraken2_output/
├── 00_pipeline_info/
│   └── software_versions.yml          # 版本与参数存档
├── 01_kraken2/
│   ├── S1.kraken                      # Kraken2 逐条分类结果
│   └── S1.kraken_report.txt           # Kraken2 分类报告(分级汇总)
├── 02_bracken/
│   ├── S1.bracken.txt                 # Bracken 丰度估算表
│   └── S1.bracken_report.txt          # Bracken 新报告
├── 03_summary/
│   ├── kraken2_summary.tsv            # 各样本分类率/优势种总表
│   └── bracken_species.tsv            # 物种丰度矩阵(样本×物种)
└── 99_logs/
    └── kraken2_pipeline.log
```

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. S1.kraken_report.txt（Kraken2 分类报告）

**通俗理解|In plain words:** 每行一个分类单元，从根到物种逐级汇总，回答「多少 reads 落到了谁头上」。

每列含义：`percentage`（占根的比例%）、`clade_reads`（该分支累积 reads 数）、`taxid`（NCBI 分类号）、`rank`（级别代码 U=未分类/R=根/D..S）、`name`（学名）。taxid=0 那行是未分类 reads。

### 2. S1.bracken.txt（Bracken 丰度表）

**通俗理解|In plain words:** 这是最常用的结果——每个物种的「真实占比」估算。

每列：`name`（物种名）、`taxid`、`kraken_assigned`（Kraken2 原始条数）、`bracken_assigned`（Bracken 重估条数）、`new_est_reads`、`fraction_total_reads`（占全部 reads 的比例，核心列）。

### 3. 03_summary/ 汇总表

**通俗理解|In plain words:** `kraken2_summary.tsv` 一屏看完所有样本的「总 reads、已分类、未分类、优势种」；`bracken_species.tsv` 是标准样本×物种丰度矩阵，直接拿去做下游（堆叠柱状图、热图、多样性分析）。

### 4. 好坏判据

- **分类率**：已分类 reads 比例（classified %）。宏基因组通常 30–80%，太低多半是数据库与样品物种差距大。
- **未分类偏高**：常见原因——样品含真核/病毒/未知物种，或数据库不含这些类群。
- **fraction_total_reads**：Bracken 后所有物种比例之和应接近 1；若某物种占比异常高，结合 kraken_report 确认是否污染。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **常规流程**：全部默认（confidence=0、bracken-level=S、threshold=10、read-len=150）。
- **要严格结果**：`--confidence 0.1`~0.2，减少低置信误分类。
- **读长非 150bp**：务必改 `--read-len`（否则 Bracken 估算偏差）。
- **要属级或株级**：`--bracken-level G`（属）或 `S1`（株）。
- **只要物种分类、不关心丰度**：加 `--no-bracken` 跳过，更快。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-dir, -i` | 必填 |  | 输入FASTQ目录｜Input FASTQ directory |
| `--db-path, -d` | 必填 |  | Kraken2数据库路径｜Kraken2 database path |
| `--output-dir, -o` | `./kraken2_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--read-len` | `150` | int | 读长(用于Bracken)｜Read length (for Bracken) |
| `--confidence` | `0.0` | float | Kraken2置信度阈值｜Kraken2 confidence score threshold |
| `--bracken-level` | `S` | D/P/C/O/F/G/S/S1 | Bracken分类级别｜Bracken taxonomic level |
| `--bracken-threshold` | `10` | int | Bracken最小读数阈值｜Bracken minimum read threshold |
| `--no-bracken` | — |  | 跳过Bracken分析｜Skip Bracken analysis |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 输入FASTQ目录｜Input FASTQ directory |
| `-d, --db-path` | 必填 |  | Kraken2数据库路径｜Kraken2 database path |
| `-o, --output-dir` | `./kraken2_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--read-len` | `150` | int | 读长(用于Bracken)｜Read length (for Bracken) |
| `--confidence` | `0.0` | float | Kraken2置信度阈值｜Kraken2 confidence score threshold |
| `--bracken-level` | `S` | D/P/C/O/F/G/S/S1 | Bracken分类级别｜Bracken taxonomic level |
| `--bracken-threshold` | `10` | int | Bracken最小读数阈值｜Bracken minimum read threshold |
| `--no-bracken` | — | store_true | 跳过Bracken分析｜Skip Bracken analysis |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Kraken2（命令行 kraken2，通过 conda 环境自动检测并 `conda run` 调用）
- Bracken（命令行 bracken，同样走 conda 环境）
- Kraken2 数据库（需预先用 kraken2-build 构建或下载，由 `-d` 指定）

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持，按样本粒度。某样本的 `01_kraken2/<样本>.kraken` 与报告已存在即跳过 Kraken2，`02_bracken/<样本>.bracken.txt` 与报告已存在即跳过 Bracken。

**Q2：只支持双端吗？**
是。当前实现只按 R1/R2 后缀配对双端，单端数据需先改后缀或另行处理；未找到任何配对样本会直接报错退出。

**Q3：为什么 Bracken 结果和 Kraken2 条数对不上？**
正常。Bracken 会根据物种基因组大小和「reads 命中分布」重新分配条数，`new_est_reads` 与 `fraction_total_reads` 才是更接近真实占比的估算。

**Q4：read-len 设错会怎样？**
Bracken 的丰度重估会偏差（因为它按读长建模）。务必与实际读长一致。

**Q5：confidence 设多少合适？**
默认 0.0（全部归类）。宏基因组常见 0.0~0.1；要求高置信度物种名单时用 0.2 左右。设太高会漏掉真实低丰度物种。