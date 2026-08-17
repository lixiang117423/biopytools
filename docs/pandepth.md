# PanDepth 覆盖度计算 | PanDepth Coverage Calculation

一句话理解：极快地算出 BAM 比对文件里"每个染色体/每个基因/每个区域被读了多少遍、覆盖了多少比例"，用来判断测序够不够深、哪里没测到。

## 功能概述 | Overview

- 极速覆盖度计算，支持单个 BAM 或整个目录的 BAM 批量处理
- 四种覆盖目标：染色体、基因(需 GFF/GTF)、自定义区域(需 BED)、滑动窗口
- 批量模式自动把各样本结果合并成一张大表，方便横向比较
- 支持比对质量/FLAG/深度过滤，可启用 GC 含量计算
- 单样本或多样本目录自动识别

## 快速开始 | Quick Start

```bash
biopytools pandepth -i sample.bam -o coverage_results
```

`-i` 可以是单个 BAM 文件，也可以是一个放多个 BAM 的目录（目录即批量模式）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BAM | 把测序 reads "贴"到参考基因组上的结果文件，记录了每条 read 落在哪里 |
| 覆盖深度(depth) | 某个位置被多少条 read 盖住，像"这里被照了几遍" |
| 覆盖比例(coverage) | 一段区域里"至少被盖到一次"的长度占比，判断有没有漏测 |
| MAPQ | read 比对的"可信度评分"，越低越可能是乱贴的 |
| FLAG | 每条 read 身上的一串"状态标签"(是否重复、是否次要比对等)，用数字位表示 |
| CDS | 基因里真正编码蛋白质的那段，做"基因覆盖度"时的默认统计对象 |
| 滑动窗口(window) | 把染色体切成一段段固定长度，逐段报深度，用于找大段的缺失/重复 |

## 输入 | Input

- **BAM 文件**：单个 `-i sample.bam`，或目录 `-i bam_files/`(批量)
- 可选 **GFF/GTF**(`-g`)：算每个基因的覆盖度
- 可选 **BED**(`-b`)：算自定义区域的覆盖度
- 可选 **窗口大小**(`-w`)：算滑动窗口覆盖度
- 三者互斥，一次只能指定一种(不指定则算染色体整体)
- 可选 `-r` 参考基因组：CRAM 解码或启用 GC 计算时必需

```text
# BED 示例(制表符分隔: 染色体 起点 终点)
chr1    1000    2000
chr1    5000    6000
```

## 参数说明 | Parameters

### 目标区域 | Target region

**通俗理解|In plain words:** 决定"覆盖度算到哪个粒度"。不给任何目标=算染色体整体；`-g`=算到每个基因(默认统计 CDS)；`-b`=算到你指定的几个区域；`-w`=切窗口逐段算。**一次只能选一种。**

### 过滤 | Filtering

**通俗理解|In plain words:** 决定"哪些 reads 不算数"。`-q` 提高=只信比对质量高的 reads；`-d` 是"统计时深度达到多少才计入覆盖"；`-x` 用 FLAG 位排除重复/质检失败/次要比对/未比对的 reads(默认 1796 即排除这四类)。**默认值一般不用动。**

### 其他 | Other

**通俗理解|In plain words:** `-r` 给参考基因组(算 GC 或解 CRAM 才需要)；`-c` 开启 GC 含量计算；`-a` 输出每个位点的深度(数据量巨大，仅在需要时开)。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入BAM文件或BAM文件目录｜Input BAM file or BAM file directory |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--gff, -g` | — | str | GFF/GTF文件用于基因覆盖度｜GFF/GTF file for gene coverage |
| `--feature, -f` | `CDS` | CDS/exon | GFF/GTF特征类型｜GFF/GTF feature type |
| `--bed, -b` | — | str | BED文件用于特定区域覆盖度｜BED file for specific region coverage |
| `--window, -w` | — | int | 滑动窗口大小(bp)｜Sliding window size in bp |
| `--min-mapq, -q` | `0` | int | 最小比对质量｜Minimum mapping quality |
| `--min-depth, -d` | `1` | int | 最小深度用于统计｜Minimum depth for statistics |
| `--exclude-flag, -x` | `1796` | int | 排除reads的FLAG标志｜FLAG bits to exclude reads |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--reference, -r` | — | str | 参考基因组文件｜Reference genome file |
| `--enable-gc, -c` | — |  | 启用GC含量计算｜Enable GC content calculation |
| `--all-sites, -a` | — |  | 输出所有位点深度｜Output all site depths |
| `--pandepth-path` | `~/software/PanDepth-2.26-Linux-x86_64/pandepth` | str | PanDepth程序路径｜PanDepth program path |
| `--verbose, -v` | — |  | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入BAM文件或BAM文件目录｜Input BAM file or BAM file directory |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-g, --gff` | — |  | GFF/GTF文件用于基因覆盖度｜GFF/GTF file for gene coverage |
| `-f, --feature` | `CDS` | CDS/exon | GFF/GTF特征类型｜GFF/GTF feature type (default: CDS) |
| `-b, --bed` | — |  | BED文件用于特定区域覆盖度｜BED file for specific region coverage |
| `-w, --window` | — | int | 滑动窗口大小(bp)｜Sliding window size in bp |
| `-q, --min-mapq` | `0` | int | 最小比对质量｜Minimum mapping quality (default: 0) |
| `-d, --min-depth` | `1` | int | 最小深度用于统计｜Minimum depth for statistics (default: 1) |
| `-x, --exclude-flag` | `1796` | int | 排除reads的FLAG标志｜FLAG bits to exclude reads (default: 1796) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-r, --reference` | — |  | 参考基因组文件(用于CRAM解码或GC计算)｜Reference genome file for CRAM decode or GC calculation |
| `-c, --enable-gc` | — | store_true | 启用GC含量计算｜Enable GC content calculation |
| `-a, --all-sites` | — | store_true | 输出所有位点深度｜Output all site depths |
| `--pandepth-path` | `~/software/PanDepth-2.26-Linux-x86_64/pandepth` |  | PanDepth程序路径｜PanDepth program path |
| `-v, --verbose` | `0` | count | 增加输出详细程度｜Increase output verbosity |
| `--quiet` | — | store_true | 静默模式，仅输出错误信息｜Quiet mode, only output errors |
| `--log-file` | — | str | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 逐个 BAM 调用 PanDepth 算覆盖度，再把同类型结果合并成一张总表。

```text
输入 BAM(单个或目录)
    │
    ▼
(批量) 扫描目录里的 .bam 文件
    │
    ▼
逐个样本调用 PanDepth 计算
    ├─ 无目标: 染色体统计 → {sample}.chr.stat.gz
    ├─ -g:      基因统计   → {sample}.gene.stat.gz
    ├─ -b:      区域统计   → {sample}.bed.stat.gz
    └─ -w:      窗口统计   → {sample}.win.stat.gz
    │
    ▼
合并同类结果 → merged_chr_statistics.txt 或 merged_gene_statistics.txt
```

## 输出 | Output

```text
coverage_results/
├── sample1.chr.stat.gz            # 单样本染色体统计(PanDepth 原生输出)
├── sample2.chr.stat.gz
├── merged_chr_statistics.txt      # 合并的染色体统计(多样本一张表)
├── sample1.gene.stat.gz           # 基因统计(使用 -g 时)
├── merged_gene_statistics.txt     # 合并的基因统计(使用 -g 时)
├── sample1.bed.stat.gz            # 区域统计(使用 -b 时)
├── sample1.win.stat.gz            # 窗口统计(使用 -w 时)
└── ...                            # 其余样本文件
```

关键文件说明：

- **`{sample}.chr.stat.gz`**：PanDepth 原生染色体统计(gzip 压缩)，含各染色体的覆盖比例与平均深度
- **`merged_chr_statistics.txt`**：把所有样本的染色体统计合并，每行一条染色体、每个样本两列(覆盖比例、平均深度)，适合直接做横向比较
- **`merged_gene_statistics.txt`**：使用 `-g` 时的基因级合并表，每行一个基因
- 合并文件名由是否有 `-g` 决定：有 `-g` 合并基因表，否则合并染色体表

## 结果解读 | Interpreting Results

### 1. 染色体统计（`merged_chr_statistics.txt`）

**通俗理解|In plain words:** 一行一条染色体，看每个样本在每条染色体上的"覆盖比例"和"平均深度"。

- **覆盖比例低(<95%)**：这条染色体有大量区域没测到，可能是缺失、异染色质区或比对困难
- **平均深度**：与测序量直接相关，深度过低的样本下游分析(如变异检测)不可靠
- 各样本间同一染色体的深度异常差异，提示样本间测序量不均

### 2. 基因统计（`merged_gene_statistics.txt`）

**通俗理解|In plain words:** 逐基因看覆盖，找出"哪些基因根本没被测到/测得很浅"。

- 覆盖比例极低的基因：可能被删除、拷贝数低，或该区域难比对
- 用于辅助判断基因存在/缺失(presence/absence)

### 3. 窗口统计（`*.win.stat.gz`）

**通俗理解|In plain words:** 把染色体切成窗口逐段报深度，深度陡降的窗口就是"缺失或大片段删除"的候选位置。

## 参数选择建议 | Parameter Guidance

- **`-q/--min-mapq`**：做覆盖度评估一般用默认 0；若担心错配干扰，可提到 20~30
- **`-d/--min-depth`**：默认 1(≥1 即算覆盖)；若只要"足够深"的区域，可提到 5~10
- **`-x/--exclude-flag`**：默认 1796 排除重复/QC失败/次要比对/未比对；一般不要改
- **`-w/--window`**：找大片段缺失/重复时常用 100000~1000000(bp)，窗口越小越细但越慢
- **`-c/--enable-gc`**：需要分析 GC 偏差时开，注意必须先给 `-r` 参考基因组

## 依赖 | Dependencies

- PanDepth 独立二进制（非 conda，默认路径 `~/software/PanDepth-2.26-Linux-x86_64/pandepth`）
- Python 3（pandas）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。重跑会覆盖已有结果；批量模式下每个 BAM 顺序执行，中途失败不会自动跳过已完成的样本。

**Q2：为什么报"只能指定一个目标区域选项"？**
`-g`、`-b`、`-w` 互斥，一次只能用一种覆盖目标。

**Q3：启用 `-c`(GC) 为什么报错？**
GC 计算需要参考基因组，必须同时给 `-r`。

**Q4：输入目录模式怎么识别样本名？**
每个 BAM 文件名(去扩展名)作为样本名，出现在合并表的列名里。

**Q5：PanDepth 程序找不到？**
确认默认路径 `~/software/PanDepth-2.26-Linux-x86_64/pandepth` 存在，或用 `--pandepth-path` 指定。
