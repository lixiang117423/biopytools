# kmertools - K-mer 工具集（建库、查询、分析） | K-mer Tools Suite (Build, Query and Analysis)

一句话理解：**围绕 k-mer 的一整套工具箱——把 FASTQ 建成可查询的 k-mer 库、按库查样本、从 FASTA 抽 k-mer、把 k-mer 丰度转成 VCF、比较两个矩阵等，十个子命令各管一件事**。

## 功能概述 | Overview

- build：用 kmtricks（默认）或 kmindex 把 FASTQ 建成可查询的 k-mer 库（RocksDB）
- query：用查询 FASTA 对 k-mer 库做查询，输出样本命中情况（matrix/json）
- extract：从 FASTA 序列中提取 k-mer（unikmer/pyfastx），输出 k-mer FASTA、位置与 BED
- count：基于 jellyfish 做 k-mer 丰度分析（滑窗统计）
- intersect：从 k-mer 丰度矩阵中提取目标 k-mer 的丰度
- kmer2vcf：把 k-mer 丰度矩阵转成 VCF 变异文件
- compare：比较两个 k-mer 矩阵文件
- gen-fof：生成 FOF 文件（样本名到文件路径的映射表）
- import-db：把 k-mer 矩阵导入 RocksDB
- split-fasta：把一个多序列 FASTA 按序列拆成单个文件

## 快速开始 | Quick Start

```bash
biopytools kmertools build -i ./fastq_dir -o ./kmer_db
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| k-mer | 从序列上滑窗切出的固定长度小片段，像把长句拆成定长词 |
| k-mer 库 | 把样本里所有 k-mer 预先建好的「字典+索引」，查起来比重新扫描快得多 |
| kmtricks / kmindex | 两个第三方 k-mer 建库/查询软件，本工具默认用 kmtricks、可选 kmindex |
| RocksDB | 一种嵌入式键值数据库，k-mer 库就存在里面 |
| FOF 文件 | File Of Files，记录「样本名 → R1/R2 文件路径」的清单，建库前先要它 |
| minimizer | 从每个窗口挑出的代表性 k-mer，用来压缩存储、加速 |
| 布隆过滤器 | 一种省空间的「可能包含」判定结构，用来快速过滤不存在的 k-mer |
| VCF | 记录变异（SNP/INDEL 等）的标准格式文件 |

## 输入 | Input

不同子命令输入不同，按需准备：

```text
build   : -i 输入目录（含 FASTQ/FASTA）
query   : -d RocksDB 库目录 + -q 查询 FASTA（或 -i kmindex 索引 + --use-kmindex）
extract : -i 输入 FASTA
count   : -i 输入目录 + -p 文件模式 + -k k-mer 库 FASTA
intersect : -m k-mer 矩阵 + -k 目标 k-mer FASTA
kmer2vcf  : -i k-mer 丰度矩阵（TSV）
compare   : -f1 矩阵1 + -f2 矩阵2
gen-fof   : --dir 输入目录
import-db : -i k-mer 矩阵（可 .gz 压缩）
split-fasta : 输入 FASTA + 输出目录（位置参数）
```

FASTQ 默认按 _1.clean.fq.gz / _2.clean.fq.gz 配对识别，其他命名可用后缀参数调整。

## 参数说明 | Parameters

### build 建库参数 | build

**通俗理解|In plain words:** -k 是 k-mer 长度（默认 51，必须 8-127 之间）；--hard-min 是最小丰度、--recurrence-min 是最小重现次数，低于阈值的 k-mer 丢弃（调大更干净但丢稀有 k-mer）；--mode/--minimizer-size 是 kmtricks 的内部参数，一般不用动。--use-kmindex 切换成 kmindex 建库。

### query 查询参数 | query

**通俗理解|In plain words:** -q 是查询序列，-o 是输出目录，-d 指到建好的 RocksDB 库（或用 -i + --use-kmindex 指 kmindex 索引）。-k 必须与建库时一致，否则查不出。--format 选 matrix（矩阵，默认）或 json。

### extract 提取参数 | extract

**通俗理解|In plain words:** 从 FASTA 里把每条序列的 k-mer 都抽出来。-k 是 k-mer 长度；--method 选提取引擎（unikmer 默认 / pyfastx）；--no-bed 不输出 BED 文件。输出文件名默认自动生成，也可用 --kmer-output/--kmer-pos-output 指定。

### count 丰度分析参数 | count

**通俗理解|In plain words:** 用 jellyfish 对一批样本统计指定 k-mer 库中每个 k-mer 的丰度。-p 是文件匹配模式，-k 是 k-mer 库 FASTA；-w/--step-size 控制滑窗统计粒度（窗口越大越粗、越快）；-C 正反链都统计。

### kmer2vcf 转换参数 | kmer2vcf

**通俗理解|In plain words:** 把 k-mer 丰度矩阵当「变异」写成 VCF。--fast-mode（默认）单遍快速处理，适合大规模；--standard-mode 三遍处理+排序更规范。--chr-length/--chr-number 用于快速模式给虚拟染色体定长度，--min-freq 过滤低频 k-mer。

### intersect / compare / import-db 参数 | intersect, compare, import-db

**通俗理解|In plain words:** intersect 从大矩阵里挑出目标 k-mer 的丰度（-w 每 N 个 k-mer 统计一次）；compare 对比两个矩阵差异；import-db 把矩阵导入 RocksDB 供 query 用（--force-overwrite 覆盖已存在库）。这些默认值基本够用。

### gen-fof / split-fasta 参数 | gen-fof, split-fasta

**通俗理解|In plain words:** gen-fof 扫描目录生成 FOF 清单（build 的输入）；split-fasta 把多序列 FASTA 拆成每序列一个 .fa 文件。都按默认后缀/规则处理，一般不用改。

## 分析流程 | Pipeline

```text
gen-fof（可选）：扫描目录 → FOF 文件
    │
    ▼
build：FASTQ + FOF → kmtricks/kmindex 建库 → RocksDB
    │
    ▼
query：查询 FASTA → 对库查询 → 样本命中矩阵
    │
    ▼
（独立工具）extract / count / kmer2vcf / intersect / compare / import-db / split-fasta
```

## 输出 | Output

各子命令都遵循 §12 目录规范，输出目录下生成 00_pipeline_info/software_versions.yml 与 99_logs/<子命令>.log：

```text
build 输出：
<输出目录>/rocksdb/              # k-mer 库（RocksDB）
<输出目录>/00_pipeline_info/software_versions.yml
<输出目录>/99_logs/build.log

extract 输出：
<输出目录>/<基名>_kmer_<k>.fa        # k-mer 序列 FASTA
<输出目录>/<基名>_kmer_<k>_pos.txt   # k-mer 位置文件
<输出目录>/<基名>.bed                 # k-mer 区间 BED（默认）
```

其他子命令（query/count/kmer2vcf/intersect/compare 等）把结果写到各自的 -o 输出文件或目录，同时写软件版本与日志。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** build 产出一个可查询的库；query 产出的矩阵回答「查询序列在哪些样本里有多少共享 k-mer」；kmer2vcf 产出的 VCF 可交给下游变异工具看。

- build 成功标志：输出目录出现 rocksdb/ 目录，日志显示「建库成功」。
- query 的 matrix 结果：行是查询序列、列是样本、值为共享 k-mer 数（或按 --threshold 过滤后的命中），值越高说明该样本与查询序列越接近。
- kmer2vcf 的 VCF：把 k-mer 当作变异位点写出，可配合其他 VCF 工具过滤/可视化。
- compare：输出两矩阵的差异统计，用于看两批数据的 k-mer 组成变化。

## 参数选择建议 | Parameter Guidance

- 建库首选默认 kmtricks；只有已用 kmindex 建好索引时才用 --use-kmindex。
- k-mer 长度：query 必须与 build 一致，否则查不出；改长度需重建库。
- 大数据建库：适当调大 -t（默认 64）并确保磁盘够放 RocksDB 与临时目录。
- kmer2vcf 默认 --fast-mode 即可；需要规范排序、可复现时用 --standard-mode。
- 只想挑一部分 k-mer：用 intersect 而非 query，直接对已有矩阵裁剪更快。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入目录(包含FASTQ/FASTA文件)｜Input directory (containing FASTQ/FASTA files) |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `--use-kmindex` | — | store_true | 使用kmindex模式 (默认: 使用kmtricks)｜Use kmindex mode (default: use kmtricks) |
| `-k, --kmer-size` | `51` | int | k-mer大小 (默认: 51)｜K-mer size (default: 51) |
| `--hard-min` | `2` | int | 最小丰度 (默认: 2)｜Minimum abundance (default: 2) |
| `--recurrence-min` | `1` | int | 最小重现次数 (默认: 1)｜Minimum recurrence (default: 1) |
| `--mode` | `kmer:pa:bin` |  | kmtricks模式 (默认: kmer:pa:bin)｜kmtricks mode (default: kmer:pa:bin) |
| `--minimizer-size` | `10` | int | minimizer大小 (默认: 10)｜Minimizer size (default: 10) |
| `--tmp-dir` | `` |  | kmtricks临时目录 (默认: 输出目录/tmp)｜kmtricks temporary directory (default: output_dir/tmp) |
| `--nb-partitions` | `0` | int | 分区数 (默认: 0=自动计算, -1=kmtricks默认)｜Partitions (default: 0=auto, -1=kmtricks default) |
| `--fof-file` | `` |  | 预存的FOF文件路径｜Pre-existing FOF file path |
| `--header-file` | `` |  | 预存的header文件路径｜Pre-existing header file path |
| `--index-name` | `kmer_index` |  | kmindex索引名称 (默认: kmer_index)｜kmindex index name (default: kmer_index) |
| `--bloom-size` | `1000000000000` | int | kmindex布隆过滤器大小 (默认: 1000000000000)｜kmindex bloom filter size (default: 1000000000000) |
| `-t, --threads` | `64` | int | 线程数 (默认: 64)｜Thread count (default: 64) |
| `--kmtricks-path` | — |  | kmtricks路径 (默认按 KMTRICKS_PATH环境变量>配置文件>内置默认 解析)｜kmtricks path (resolved via KMTRICKS_PATH env>config>built-in) |
| `--kmindex-path` | — |  | kmindex路径 (默认按 KMINDEX_PATH环境变量>配置文件>内置默认 解析)｜kmindex path (resolved via KMINDEX_PATH env>config>built-in) |
| `--bgzip-path` | — |  | bgzip路径 (默认按 BGZIP_PATH环境变量>配置文件>内置默认 解析)｜bgzip path (resolved via BGZIP_PATH env>config>built-in) |
| `--fof-suffix-1` | `_1.clean.fq.gz` |  | R1文件后缀 (默认: _1.clean.fq.gz)｜R1 file suffix (default: _1.clean.fq.gz) |
| `--fof-suffix-2` | `_2.clean.fq.gz` |  | R2文件后缀 (默认: _2.clean.fq.gz)｜R2 file suffix (default: _2.clean.fq.gz) |
| `-d, --database` | — |  | RocksDB数据库目录 (kmtricks模式必需，默认模式)｜RocksDB database directory (required for kmtricks mode, default) |
| `-i, --index` | — |  | kmindex索引目录 (kmindex模式必需，需加--use-kmindex)｜kmindex index directory (required for kmindex mode, requires --use-kmindex) |
| `-q, --query` | 必填 |  | 查询FASTA文件｜Query FASTA file |
| `--header-db-key` | `kmer_header` |  | 数据库中的header key (默认: kmer_header)｜Header key in database (default: kmer_header) |
| `--bed` | — |  | BED文件路径（用于生成位置丰度文件）｜BED file path (for generating position-abundance file) |
| `--zvalue` | `0` | int | findere算法z值 (默认: 0)｜findere z-value (default: 0) |
| `--threshold` | `0.0` | float | 共享k-mer阈值 (默认: 0.0)｜Shared k-mer threshold (default: 0.0) |
| `--format` | `matrix` | json/matrix | 输出格式 (默认: matrix)｜Output format (default: matrix) |
| `--dir` | 必填 |  | 输入目录｜Input directory |
| `--suffix-1` | `_1.clean.fq.gz` |  | R1文件后缀 (默认: _1.clean.fq.gz)｜R1 file suffix (default: _1.clean.fq.gz) |
| `--suffix-2` | `_2.clean.fq.gz` |  | R2文件后缀 (默认: _2.clean.fq.gz)｜R2 file suffix (default: _2.clean.fq.gz) |
| `--input-delimiter` | `	` |  | 输入分隔符 (默认: tab)｜Input delimiter (default: tab) |
| `--batch-size` | `20000` | int | 批量写入大小 (默认: 20000)｜Batch write size (default: 20000) |
| `--bloom-bits` | `15` | int | Bloom filter位数 (默认: 15)｜Bloom filter bits per key (default: 15) |
| `--force-overwrite` | — | store_true | 强制覆盖已存在的数据库｜Force overwrite existing database |
| `--method` | `unikmer` | unikmer/pyfastx | 提取方法 (默认: unikmer)｜Extraction method (default: unikmer) |
| `--unikmer-path` | — |  | unikmer路径 (默认按 UNIKMER_PATH环境变量>配置文件>内置默认 解析)｜unikmer path (resolved via UNIKMER_PATH env>config>built-in) |
| `--kmer-output` | — |  | kmer FASTA文件 (默认: output_dir/basename_kmer_k.fa)｜Kmer FASTA file |
| `--kmer-pos-output` | — |  | kmer位置文件 (默认: output_dir/basename_kmer_k_pos.txt)｜Kmer position file |
| `--no-bed` | — | store_true | 不输出BED文件｜Do not output BED file |
| `-p, --pattern` | 必填 |  | 文件模式，支持FASTQ和FASTA格式｜File pattern (FASTQ/FASTA) |
| `-k, --kmer-lib` | 必填 |  | K-mer库文件(FASTA格式)｜K-mer library file (FASTA format) |
| `-b, --bed-file` | — |  | BED文件路径｜BED file path |
| `-m, --kmer-size` | `51` | int | K-mer长度 (默认: 51)｜K-mer size (default: 51) |
| `-s, --hash-size` | `1000M` |  | 哈希表大小 (默认: 1000M)｜Hash table size (default: 1000M) |
| `-w, --window-size` | `500000` | int | 滑动窗口大小bp (默认: 500000)｜Window size in bp (default: 500000) |
| `--step-size` | — | int | 滑动窗口步长bp (默认: window-size/5)｜Step size in bp (default: window-size/5) |
| `-C, --canonical` | — | store_true | 统计正向和反向互补链｜Count both strands |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep temporary files |
| `--keep-binary` | — | store_true | 保留0/1存在缺失矩阵｜Keep 0/1 matrix |
| `--jellyfish-path` | — |  | Jellyfish路径 (默认按 JELLYFISH_PATH环境变量>配置文件>内置默认 解析)｜Jellyfish path (resolved via JELLYFISH_PATH env>config>built-in) |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |
| `-i, --input-matrix` | 必填 |  | 输入kmer丰度矩阵文件｜Input kmer abundance matrix file (TSV format) |
| `-o, --output-vcf` | 必填 |  | 输出VCF文件路径｜Output VCF file path (.vcf or .vcf.gz) |
| `--fast-mode` | `True` | store_true | 快速模式（单次遍历，默认）｜Fast mode (single pass, default) |
| `--standard-mode` | — | store_true | 标准模式（3遍处理+排序）｜Standard mode (3-pass processing + sorting) |
| `--chr-length` | `100000000` | int | 每条染色体长度（快速模式，默认100M）｜Chromosome length for fast mode (default: 100M) |
| `--chr-number` | `0` | int | 染色体数量（如果设置则优先使用）｜Number of chromosomes (if set, takes priority over --chr-length) |
| `--min-freq` | `0` | int | 最小出现频次过滤（快速模式）｜Minimum frequency filter for fast mode (default: 0=no filter) |
| `--kmer-length` | `51` | int | Kmer长度，用于VCF INFO字段｜Kmer length for VCF INFO field (default: 51) |
| `--no-header` | — | store_true | 输入文件没有header行（第一行就是样本名）｜Input file has no header line (first line is sample names) |
| `-m, --min-agg-count` | `3` | int | 最小聚合频次阈值（标准模式）｜Minimum aggregated count threshold for standard mode (default: 3) |
| `-T, --temp-dir` | — |  | 临时文件目录（标准模式）｜Temporary directory for standard mode (default: ./temp) |
| `-m, --kmer-matrix` | 必填 |  | kmer矩阵文件｜Kmer matrix file (TSV format) |
| `-k, --kmer-fasta` | 必填 |  | 目标kmer的fasta文件｜Target kmer fasta file |
| `-o, --output-file` | 必填 |  | 输出文件路径｜Output file path |
| `-r, --use-reverse-complement` | `True` | store_true | 使用反向互补查询 (默认: 启用)｜Use reverse complement query (default: enabled) |
| `--no-reverse-complement` | — | store_false | 不使用反向互补查询｜Do not use reverse complement query |
| `--keep-not-found` | `True` | store_true | 保留未找到的kmer (默认: 保留)｜Keep kmers not found (default: enabled) |
| `--no-keep-not-found` | — | store_false | 不保留未找到的kmer｜Do not keep kmers not found |
| `-f, --output-format` | `tsv` | tsv/csv | 输出格式 (默认: tsv)｜Output format (default: tsv) |
| `-f1, --file1` | 必填 |  | 第一个kmer矩阵文件｜First kmer matrix file |
| `-f2, --file2` | 必填 |  | 第二个kmer矩阵文件｜Second kmer matrix file |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- kmtricks（默认 ~/miniforge3/envs/pan/bin/kmtricks）、kmindex（~/miniforge3/envs/pan/bin/kmindex）、bgzip
- unikmer（extract 默认方法）、jellyfish（count）、seqkit 等按需
- Python 库：python-rocksdb（query/import-db）、biopython（split-fasta）、pandas（count/kmer2vcf）、yaml
- 工具路径按「环境变量 > 配置文件 ~/.config/biopytools/config.yml > 内置默认」解析

## 常见问题 | FAQ

**Q1：有没有断点续传？**
各子命令基本是一次性运行，无跨运行断点续传。中断后重跑会从头执行；建库/导入可用 --force-overwrite 明确覆盖旧结果。

**Q2：query 结果全是 0 或找不到？**
最常见原因是 -k 与建库时不一致，或 -d 指错了 RocksDB 目录。请确认查询 k-mer 长度与建库一致、路径正确。

**Q3：build 很慢或磁盘爆了？**
kmtricks 建库本身吃 CPU/内存/磁盘。确认 -t 合理、临时目录（默认 输出目录/tmp）所在分区有足够空间，必要时用 --tmp-dir 指到大空间分区。

**Q4：kmtricks 找不到？**
检查 kmtricks 是否安装，或用 --kmtricks-path 指定路径；也可设环境变量 KMTRICKS_PATH 或写入 ~/.config/biopytools/config.yml。

**Q5：gen-fof 生成的样本名为什么被「清洗」了？**
kmtricks 要求样本名只能是字母数字，gen-fof 会自动去掉下划线等特殊字符，属正常行为，不是数据丢失。
