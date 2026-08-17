# 启动子提取 | Promoter Extractor

一句话理解：**根据基因注释(GFF3)和基因组序列，把每个基因上游(反义链则是下游)固定长度的「启动子」DNA 片段批量切出来，输出 FASTA 序列和 BED 坐标。**

## 功能概述 | Overview { #overview }

- 从 GFF3 的 gene 特征 + 基因组 FASTA 提取启动子序列
- 自动处理基因方向：正链取上游、反链取下游并做反向互补
- 边界自动截断，支持最小长度过滤，避免输出太短的片段
- 可选只提取指定基因 ID 列表
- 纯 Python 实现，无外部生信软件依赖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools promoter-extractor -g genes.gff3 --genome genome.fa -o promoters
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 启动子(promoter) | 基因前面的「开关区」，控制基因何时、在哪儿表达 |
| GFF3 | 记录基因位置与结构的标准注释格式 |
| 正链 / 反义链 | DNA 是双链，基因可能落在任一链上 |
| 反向互补 | 反义链基因的启动子要按相反方向读出，才是正确的 5 端到 3 端 |
| BED | 记录「基因组上某区间」的坐标格式 |

## 输入 | Input { #input }

两个必需文件 + 一个可选列表：

- `--gff`：GFF3 注释文件，只解析 `gene` 特征，且 gene 必须带 `ID`(或 `gene_id`)属性
- `--genome`：基因组 FASTA，序列名需与 GFF3 的 seqid 一致
- `--gene-list`(可选)：每行一个基因 ID，只提取列表里的基因

```text
# 基因列表文件示例(可选)
gene1
gene2
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required { #parameters-required }

**通俗理解|In plain words:** 两个必需输入(GFF3 + 基因组 FASTA)。`-o` 注意是「输出前缀」而不是目录，产物会写成 `<前缀>.fa`、`<前缀>.bed` 等，放在当前目录(或前缀里带路径时放对应目录)。

### 启动子参数 | Promoter { #parameters-promoter }

**通俗理解|In plain words:** `--promoter-length` 决定切多长，默认 2000bp(常见值)。`--min-length` 是最小接受长度，默认 0 表示「再短也要」；把它调大(如 500)可以过滤掉那些因靠近染色体末端被截得只剩一丁点的启动子。

### 基因选择 | Gene selection { #parameters-gene }

**通俗理解|In plain words:** `--gene-list` 只提取指定基因。不传就是「全部基因都提」。

### 输出开关 | Output toggles { #parameters-output }

**通俗理解|In plain words:** 默认 FASTA + BED + 统计三样都出。`--no-bed` / `--no-stats` 分别关掉 BED 和统计文件。

### 运行控制 | Run control { #parameters-run }

**通俗理解|In plain words:** `--force` 强制覆盖已有文件；`--dry-run` 只演练不真正执行；`-v` / `--quiet` 调节日志详略。`--threads` 目前默认 1，实际提取是纯 Python 顺序处理，未真正并行，一般不用管。

## 分析流程 | Pipeline { #pipeline }

```text
输入 GFF3 + 基因组 FASTA
    │
    ▼
解析 GFF3 的 gene 特征(ID / 位置 / 链方向)
    │
    ▼
计算启动子区间：正链 [start-N, start-1]，反链 [end+1, end+N]
    │
    ▼
边界截断到染色体范围 + 最小长度过滤
    │
    ▼
从基因组切序列(反义链做反向互补)
    │
    ▼
输出 FASTA + BED + 统计
```

## 输出 | Output { #output }

```text
promoters.fa                 # 启动子序列 FASTA(每行 80 bp 折行)
promoters.bed                # 启动子坐标 BED(0-based)
promoters_stats.txt          # 提取统计
promoter_extractor.log       # 运行日志
```

- `promoters.fa`：每条记录名为 `<基因ID>_promoter`，序列已按链方向校正(反义链为反向互补)
- `promoters.bed`：BED6 格式(seqid、start、end、基因ID、长度、链)
- `promoters_stats.txt`：总基因数、成功提取数、边界截断数、长度不足数、成功率

## 结果解读 | Interpreting Results { #interpreting-results }

- 看 `promoters_stats.txt` 的「成功率」：接近 100% 正常；「边界截断」多说明很多基因靠近染色体末端，启动子切不满 2000bp
- 「长度不足」数量由 `--min-length` 决定，默认 0 时基本为 0
- 反义链基因的启动子序列已做反向互补，拿到手直接就是 5 端到 3 端方向，可直接做下游 motif 分析
- 若某基因没出现在输出里：检查 GFF3 里它是不是 `gene` 特征、有没有 `ID` 属性、seqid 是否与基因组 FASTA 序列名一致

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规启动子分析用默认 2000bp 即可
- 只想研究近端启动子(如转录因子结合核心区)可把 `--promoter-length` 调到 500~1000
- 想剔除被边界截断太短的：把 `--min-length` 设成接近 `--promoter-length` 的值(如 1500)
- 只提取一批目标基因：用 `--gene-list` 传 ID 列表

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --gff` | 必填 | Path | 输入GFF3文件路径｜Input GFF3 file path |
| `--genome` | 必填 | Path | 输入基因组FASTA文件路径｜Input genome FASTA file path |
| `-o, --output` | `promoters` |  | 输出前缀｜Output prefix (default: promoters) |
| `-p, --promoter-length` | `2000` | int | 启动子长度（bp）｜Promoter length in bp (default: 2000) |
| `--min-length` | `0` | int | 最小接受长度（bp）｜Minimum acceptable length in bp (default: 0) |
| `--gene-list` | — | Path | 基因ID列表文件｜Gene ID list file (one gene ID per line) |
| `--no-bed` | — |  | 不输出BED格式文件｜Do not output BED format file |
| `--no-stats` | — |  | 不输出统计文件｜Do not output statistics file |
| `-t, --threads` | `1` | int | 线程数｜Number of threads (default: 1) |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--force, -f` | — |  | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — |  | 模拟运行(不实际执行)｜Dry run without execution |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --gff` | 必填 |  | 输入GFF3文件路径｜Input GFF3 file path |
| `-g, --genome` | 必填 |  | 输入基因组FASTA文件路径｜Input genome FASTA file path |
| `-o, --output` | `promoters` |  | 输出前缀｜Output prefix (default: promoters) |
| `-p, --promoter-length` | `2000` | int | 启动子长度（bp）｜Promoter length in bp (default: 2000) |
| `--min-length` | `0` | int | 最小接受长度（bp）｜Minimum acceptable length in bp (default: 0) |
| `--gene-list` | — |  | 基因ID列表文件｜Gene ID list file (one gene ID per line) |
| `--no-bed` | — | store_true | 不输出BED格式文件｜Do not output BED format file |
| `--no-stats` | — | store_true | 不输出统计文件｜Do not output statistics file |
| `-t, --threads` | `1` | int | 线程数｜Number of threads (default: 1) |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level (default: INFO) |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run without execution |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- 纯 Python 3 实现，无外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
不支持。这是单步流程，重新运行会覆盖同名输出文件。

**Q2：`-o` 是目录还是前缀？**
是「输出前缀」，不是目录。产物写到当前目录；想放子目录就写 `-o out/promoters`，会写到 `out/` 下。

**Q3：为什么有的基因没被提取？**
常见原因：GFF3 里它不是 `gene` 特征、缺 `ID` 属性、seqid 与基因组序列名对不上、或长度不足被 `--min-length` 过滤。

**Q4：反义链基因的启动子方向对吗？**
对。反义链基因的启动子会做反向互补，输出即正确的 5 端到 3 端方向。

**Q5：`--threads` 能加速吗？**
目前提取是纯 Python 顺序处理，`--threads` 未真正用于并行，默认 1 即可。