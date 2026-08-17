# kmc - KMC k-mer 统计与丰度矩阵分析 | KMC K-mer Counting and Abundance Matrix

一句话理解：**把一堆测序数据（FASTQ/FASTA）切成固定长度的小片段（k-mer）数一数各出现多少次，再汇总成一张「k-mer × 样本」的丰度大表，支持查某个 k-mer 在哪些样本里有多少**。

## 功能概述 | Overview

- 用 KMC 软件为每个样本建立 k-mer 数据库（计数），支持目录自动识别双末端/单末端
- 把多个样本的计数结果合并成跨样本丰度矩阵（HDF5 存储，自动导出 TSV）
- 支持从 FASTA 批量查询 k-mer 在各样本的丰度（百万级）
- 支持增量式添加新样本（add）与矩阵导出（export）
- 稀疏存储默认开启，内存占用可控；可选内存/RocksDB 索引加速查询

## 快速开始 | Quick Start

```bash
biopytools kmc count -d raw_data/ -o kmc_output
```

然后构建丰度矩阵：

```bash
biopytools kmc matrix -o kmc_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| k-mer | 从序列上滑窗切出的、长度为 k 的小片段，像把长句拆成固定字数的词 |
| k-mer 计数 | 数每个「词」在数据里出现多少次，出现次数即「丰度」 |
| 丰度矩阵 | 一张大表：行是 k-mer，列是样本，格子填该 k-mer 在该样本的出现次数 |
| 双末端测序 | 一个样本有 R1/R2 两个文件，互为配对，需一起统计 |
| 最小计数阈值 | 出现次数低于该值的 k-mer 视为测序错误丢弃，像「少数服从多数」 |
| 稀疏矩阵 | 只存非零的格子（大多数 k-mer 只出现在少数样本），省内存 |
| HDF5 | 一种能高效存大矩阵的文件格式，程序内部用它，导出时才变 TSV |

## 输入 | Input

count 步骤支持两种输入方式：

```text
目录模式（自动识别配对）：-d raw_data/   # 按 read1/read2 后缀自动配对
文件列表模式：-i s1.fq -i s2.fq       # 逐个指定文件
```

默认后缀：Read1 为 _1.clean.fq.gz，Read2 为 _2.clean.fq.gz（可用 --read1-suffix/--read2-suffix 改）；单末端数据加 --single-end。

query 步骤输入一个 FASTA 文件（每行一条查询序列），输出各序列在各样本的丰度。

## 参数说明 | Parameters

### 计数核心参数 | Counting core

**通俗理解|In plain words:** -k 是 k-mer 长度（切多长的片段），--min-count 是最小计数阈值（低于此数的 k-mer 丢掉），--max-count 是最大计数上限（可选）。k 太小片段重复多、信息少，太大则片段难重复；21 是常用值。--min-count 调大能去更多测序错误、但会丢稀有 k-mer，一般 2 即可。

### 输入与配对参数 | Input and pairing

**通俗理解|In plain words:** -d 是输入目录（自动按后缀找样本并配对），-i 是逐个指定文件，-n 是手动指定样本名。read1/read2 后缀决定怎么配对，换用别的命名规范时才需要改。--single-end 告诉程序这是单末端数据、不要找配对。

### 路径与性能参数 | Paths and performance

**通俗理解|In plain words:** -o 输出目录、--tmp-dir 临时目录、--kmc-path KMC 软件位置，都有默认值。--threads 线程数默认 12，机器核多可调大。--memory-limit 限制 KMC 内存（如 12G），机器内存紧张时才需要设。

### 矩阵构建参数 | Matrix building

**通俗理解|In plain words:** matrix 步骤把各样本的计数合成一张大表。--max-memory 是矩阵构建允许用的内存上限（GB，默认 500），--dense-matrix 换成密集矩阵（全表存满，仅样本少时建议），--matrix-format 选存储格式（默认 hdf5）。--no-export 关闭自动导出 TSV；--no-keep-dump 不保留中间 dump 文件省空间。

### 查询与导出参数 | Query and export

**通俗理解|In plain words:** query 用 -f 给一个 FASTA 查询文件，--db-dir 指到建好的矩阵目录，结果可写文件或打印到屏幕。export 用 --min-abundance 过滤掉丰度太低的 k-mer，--format full/sparse 选输出格式。

## 分析流程 | Pipeline

```text
count：FASTQ/FASTA → KMC 计数 → 每样本数据库（kmc_databases/）
    │
    ▼
matrix：合并各样本数据库 → 全局 k-mer 字典 + 丰度矩阵（HDF5）
    │
    ▼
自动导出 abundance_matrix.tsv + presence_matrix.tsv
    │
    ▼
query：FASTA 序列 → 查各样本丰度（可选，随时执行）
```

## 输出 | Output

```text
kmc_output/
├── kmc_databases/                  # 每样本 KMC 数据库（<样本>.kmc_pre/.kmc_suf）
├── abundance_matrix.h5             # 丰度矩阵（稀疏或密集存储）
├── kmer_dictionary.h5              # k-mer 序列字典（k-mer ↔ 编号映射）
├── kmer_metadata.json              # 元数据（k-mer 长度、索引类型等）
├── global_kmers.txt                # 合并后的全部 k-mer 列表
├── global_kmers_with_id.txt        # 带编号的 k-mer 列表
├── global_kmers.kmc_pre/.kmc_suf   # 合并后的全局 KMC 数据库
├── abundance_matrix.tsv            # 自动导出的丰度矩阵（TSV）
├── presence_matrix.tsv             # 自动导出的存在/缺失矩阵（0/1）
├── dump_files/                     # 每样本 dump 中间文件（可删，--no-keep-dump 不保留）
└── kmc_analysis.log                # 运行日志
```

- abundance_matrix.h5：核心结果，k-mer × 样本 的丰度矩阵，供后续分析用。
- kmer_dictionary.h5：k-mer 序列与内部编号的映射，查询时要用。
- abundance_matrix.tsv：人类可读的丰度表，第一列 k-mer、后续列样本。
- presence_matrix.tsv：0/1 存在矩阵，只关心「有没有」而非「多少」时看它。
- kmc_databases/：每样本的 KMC 计数数据库，是 matrix 步骤的输入。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 核心看 abundance_matrix.tsv——某一行是某个 k-mer，某一列是某个样本，格子的数字就是「这个 k-mer 在这个样本里出现了几次」。

- 丰度矩阵行数（k-mer 总数）反映数据的 k-mer 多样性；样本列数等于参与分析的样本数。
- presence_matrix.tsv 里 1 表示该样本含这个 k-mer，0 表示没有；常用于「某 k-mer 在哪些样本出现」的筛选。
- query 结果中「未找到」的 k-mer 说明该序列不在任何样本的数据库中（可能因 --min-count 被过滤或数据本身没有）。
- 日志里的矩阵统计会打印 k-mer 数、样本数、k-mer 长度、存储类型和稀疏度，可据此核对参数是否正确。

## 参数选择建议 | Parameter Guidance

- k-mer 长度 -k：常规 21 或 31，视基因组大小/分析目的；做 k-mer 谱评估常用 21，去重/聚类可更大。
- --min-count：默认 2 即可去除大部分测序错误；高覆盖数据可适当调大（如 5）进一步降噪。
- 样本极多、k-mer 极多时：保持默认稀疏存储；仅当样本很少（几个）且要全表运算时用 --dense-matrix。
- 想省磁盘：matrix 步骤加 --no-keep-dump（不保留 dump_files/）。
- 只需 0/1 存在信息：直接用 presence_matrix.tsv，不必看丰度矩阵。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-dir, -d` | — |  | 输入目录(自动识别双末端测序)｜Input directory (auto-detect paired-end) |
| `--input, -i` | — |  | 输入文件(FASTQ/FASTA)，可多次使用｜Input files (FASTQ/FASTA), can be used multiple times |
| `--sample-names, -n` | — |  | 样本名称，可多次使用(默认使用文件名)｜Sample names, can be used multiple times (default: use filename) |
| `--read1-suffix` | `_1.clean.fq.gz` |  | Read1文件后缀｜Read1 file suffix |
| `--read2-suffix` | `_2.clean.fq.gz` |  | Read2文件后缀｜Read2 file suffix |
| `--single-end` | — |  | 单末端测序模式｜Single-end sequencing mode |
| `--input-dir, -i` | — |  | 包含kmc_databases的目录(即count步骤的-o参数)｜Directory containing kmc_databases (i.e., -o from count step) |
| `--max-memory` | `500` | int | 最大内存使用量(GB)｜Maximum memory usage (GB) |
| `--matrix-format` | `hdf5` | hdf5/tsv/sqlite | 矩阵存储格式｜Matrix storage format |
| `--dense-matrix` | — |  | 使用密集矩阵(默认稀疏)｜Use dense matrix (default: sparse) |
| `--no-export` | — |  | 不导出TSV文件(默认自动导出)｜Do not export TSV files (auto export by default) |
| `--keep-dump` | — |  | 保留dump文件到dump_files目录(默认保留)｜Keep dump files in dump_files directory (default: True) |
| `--no-keep-dump` | — |  | 不保留dump文件(节省空间)｜Do not keep dump files (save space) |
| `--input-fasta, -f` | 必填 |  | 查询的FASTA文件(支持百万级批量查询)｜FASTA file for query (supports millions) |
| `--output-file, -o` | — |  | 输出文件(可选，默认输出到屏幕)｜Output file (optional, default to screen) |
| `--db-dir, --output-dir` | `./kmc_output` |  | KMC数据库目录(包含abundance_matrix.h5等文件)｜KMC database directory (contains abundance_matrix.h5, etc.) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--format` | `sparse` | full/sparse | 输出格式｜Output format: full: 完整矩阵(所有k-mer x 所有样本)｜full matrix (all k-mers x all samples) sparse: 稀疏格式(只包含非零值)｜sparse format (non-zero values only) |
| `--min-abundance` | `1` | int | 最小丰度阈值｜Minimum abundance threshold |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-m, --mode` | `count` | count/matrix/query/add/export | 操作模式｜Operation mode: count: 统计k-mer｜count k-mers matrix: 构建丰度矩阵｜build abundance matrix query: 查询k-mer｜query k-mer add: 添加新样本｜add new samples export: 导出矩阵为TSV｜export matrix to TSV |
| `-d, --input-dir` | — |  | 输入目录(自动识别双末端测序)｜Input directory (auto-detect paired-end) |
| `-i, --input` | — |  | 输入文件列表(FASTQ/FASTA)｜Input files list (FASTQ/FASTA) |
| `-n, --sample-names` | — |  | 样本名称(默认使用文件名)｜Sample names (default: use filename) |
| `--read1-suffix` | `_1.clean.fq.gz` |  | Read1文件后缀｜Read1 file suffix |
| `--read2-suffix` | `_2.clean.fq.gz` |  | Read2文件后缀｜Read2 file suffix |
| `--single-end` | — | store_true | 单末端测序模式｜Single-end sequencing mode |
| `-k, --kmer-size` | `21` | int | k-mer大小｜k-mer size |
| `--min-count` | `2` | int | 最小计数阈值｜Minimum count threshold |
| `--max-count` | — | int | 最大计数阈值｜Maximum count threshold |
| `--kmc-path` | `~/miniforge3/envs/kmc_v.3.2.4/bin` |  | KMC软件路径｜KMC software path |
| `-o, --output-dir` | `./kmc_output` |  | 输出目录｜Output directory |
| `--tmp-dir` | `./kmc_tmp` |  | 临时文件目录｜Temporary directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--memory-limit` | — |  | 内存限制(如: 12G)｜Memory limit (e.g., 12G) |
| `--max-memory` | `500` | int | 最大内存使用量(GB)｜Maximum memory usage (GB) |
| `--matrix-format` | `hdf5` | hdf5/tsv/sqlite | 矩阵存储格式｜Matrix storage format |
| `--dense-matrix` | — | store_true | 使用密集矩阵(默认稀疏)｜Use dense matrix (default: sparse) |
| `--no-export` | — | store_true | 不自动导出TSV文件(默认自动导出丰度和存在文件)｜Do not auto export TSV files (default: auto export abundance and presence files) |
| `--no-keep-dump` | — | store_true | 不保留dump文件(默认保留到dump_files目录)｜Do not keep dump files (default: keep in dump_files directory) |
| `--index-mode` | `auto` | auto/memory/db | k-mer索引模式｜K-mer index mode:\nauto: 自动选择(默认)｜auto: Auto select based on file size (default)\nmemory: 使用内存索引(快但占用大内存)｜memory: Use in-memory index (fast but high memory)\ndb: 使用数据库索引(省内存但稍慢)｜db: Use database index (low memory but slightly slower) |
| `--index-threshold` | `1.0` | float | 自动选择索引的文件大小阈值(GB)｜File size threshold for auto index selection (GB) |
| `-f, --input-fasta` | — |  | 查询的FASTA文件(批量查询)｜FASTA file for batch query |
| `--output-file` | — |  | 查询结果输出文件/导出TSV文件｜Query result output file / Export TSV file |
| `--format` | `sparse` | full/sparse | 导出格式｜Export format: full (完整矩阵｜full matrix) or sparse (稀疏格式｜sparse format) |
| `--min-abundance` | `1` | int | 最小丰度阈值｜Minimum abundance threshold |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- KMC 软件（kmc / kmc_tools / kmc_dump，默认路径 ~/miniforge3/envs/kmc_v.3.2.4/bin）
- Python 库：h5py、numpy（矩阵存储）；python-rocksdb 或 plyvel（数据库索引，可选）；biopython（query 读 FASTA）
- 无需额外 conda 环境（KMC 通过 conda run 自动包装）

## 常见问题 | FAQ

**Q1：有没有断点续传？**
count 步骤无跳过逻辑（重跑会重数）；matrix 步骤会复用已存在的 dump 文件（dump_files/<样本>_dump.txt 存在就跳过重新 dump），add 步骤是增量式。中断后重跑建议保留原输出目录。

**Q2：换 k-mer 大小重跑要注意什么？**
k-mer 大小会影响所有下游结果。换 -k 后旧矩阵和数据库不兼容，建议清空输出目录或换新目录重跑，避免新旧结果混用。

**Q3：为什么 query 很多 k-mer 未找到？**
可能原因：该 k-mer 在计数时低于 --min-count 被过滤；或查询序列与建库时的 k-mer 长度不一致（查询前需确认 -k 一致）。

**Q4：export 子命令为什么报「无效的操作模式」？**
帮助里列出的 export 子命令当前未被配置校验覆盖（有效模式只含 count/matrix/query/add），直接用会报错。需要导出 TSV 时，matrix/add 步骤默认已自动生成 abundance_matrix.tsv 和 presence_matrix.tsv，一般够用。

**Q5：内存不够怎么办？**
调小 --max-memory（如 200）；样本多、k-mer 多时用稀疏存储（默认）配合数据库索引（RocksDB）可显著降低查询内存占用。
