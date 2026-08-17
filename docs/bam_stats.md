# bam_stats - BAM 批量统计分析 | BAM File Batch Statistics

一句话理解：**给一批 BAM 各做一次"体检"，输出一张覆盖比对率、测序深度、GC、插入片段、重复率、变异/污染预警等指标的总表**，让你一眼看清每个样本测得好不好、有没有污染。

## 功能概述 | Overview

- 批量处理多个 BAM（单文件或目录），用多进程并行加速
- 六类统计：**比对(alignment)、覆盖度(coverage)、序列特征(sequence)、插入片段(insert)、重复(duplicate)、变异/特异性(variation)**，每类可用 `--skip-*` 单独跳过
- 输出三份：全局长表（样本为行、指标为列）、染色体级别表、基因组级 JSON 汇总
- 内置污染预警（低比对率、非主要染色体异常高覆盖）
- 无断点续传（每次重算）

## 快速开始 | Quick Start

~~~bash
biopytools bam-stats -i ./bam_files -o result.summary.tsv
~~~

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 比对率(Map Rate) | 读段里有多少能"贴回"参考基因组，低说明样品/参考不匹配或污染 |
| Proper Pair | 一对读段方向、距离都正常的"合法配对"，比例低提示文库质量或插入片段异常 |
| 覆盖度(深度) | 每个碱基平均被测多少次，深=测得多 |
| 插入片段 | 一对读段之间夹着的那段 DNA 长度，反映文库构建的片段大小 |
| 重复率 | 同一位置出现多条一模一样的读段，PCR 过度或测序过深都会抬高它 |
| 软剪切(soft clip) | 读段两端"比对不上被剪掉"的部分，比例高提示变异或接头污染 |
| GC 偏差 | 不同 GC 含量的区域被测得"不公平"，有参考基因组才能算 |

## 输入 | Input

- 单个 BAM（或 SAM）文件，或包含 BAM 的目录（自动收集 `*.bam`/`*.BAM`）
- 可选 `-g 参考基因组 FASTA`（用于计算 GC 偏差，不给则变异统计里的 GC 项缺失）
- 可选 `--bed-file` 限定只在目标区域统计

~~~text
biopytools bam-stats -i sample.bam -o result.summary.tsv -g reference.fa
biopytools bam-stats -i ./bam_dir -o result.summary.tsv --bed-file targets.bed
~~~

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 是 BAM 文件或目录，`-o` 是汇总表输出路径（默认 `bam_stats.summary.tsv`）。`-g` 参考基因组、`--bed-file` 目标区域是可选增强，给了能算 GC 偏差、能把统计限定在感兴趣区域。

### 过滤与阈值 | Filters

**通俗理解|In plain words:** `--min-mapq`（默认 20）是最低比对质量，`--max-insert`（默认 1000）是插入片段上限。调大 `--min-mapq` 只统计更可靠读段；`--max-insert` 过小会漏掉大片段文库的正常读段。**常规全基因组质控用默认即可。**

### 统计模块开关 | Module switches

**通俗理解|In plain words:** 六个 `--skip-*` 开关分别跳过某类统计（alignment/coverage/sequence/insert/duplicate/variation）。**只想快速看比对率或深度时，跳过不需要的类别能显著提速**；全量分析则一个都不用加。

### 性能 | Performance

**通俗理解|In plain words:** `-t`（默认 12）是给 `samtools` 的线程数，`-p`（默认 16）是同时并行处理的样本数。样本多就调大 `-p`，注意与 `-t` 相乘别超过机器核数。

## 分析流程 | Pipeline

~~~text
输入 BAM(文件/目录) [+ 参考FASTA] [+ BED]
    │
    ▼
1. 检查依赖(samtools/bedtools + pandas/openpyxl/tqdm)
2. 收集 BAM 列表, 缺索引则 samtools index 补建
3. 多进程并行, 逐样本跑六类统计:
   比对 | 覆盖度 | 序列特征 | 插入片段 | 重复 | 变异/特异性
4. 生成报告:
   {prefix}.summary.tsv        (全局长表)
   {prefix}.per_chromosome.tsv (染色体级)
   {prefix}.genome_stats.json  (基因组级JSON)
~~~

## 输出 | Output

~~~text
输出目录/
├── result.summary.tsv               # 全局长表(样本为行, 指标为列) —— 核心结果
├── result.per_chromosome.tsv        # 染色体级统计(每样本每染色体一行)
├── result.genome_stats.json         # 基因组级汇总(JSON)
└── 99_logs/
    └── pipeline.log                 # 运行日志
~~~

全局长表 `result.summary.tsv` 主要列：

| 分组 | 列 |
|------|----|
| 比对 | Total_Reads / Mapped_Reads / Unmapped_Reads / Map_Rate(%) / Proper_Pair(%) / Unique_Mapped / Multi_Mapped / Singletons |
| 覆盖度 | Avg_Depth / Median_Depth / Depth_Std / Depth_CV / Cov_Bases_Rate(%) |
| 重复 | Duplicates_Rate(%) / Duplicates_Marked |
| 插入片段 | Mean_Insert_Size / Median_Insert_Size / Insert_Size_Std |
| 序列 | Avg_Read_Length / Avg_GC(%) |
| 变异 | Mismatch_Rate(%) / SoftClip_Rate(%) / SoftClip_Reads_Rate(%) / GC_Bias(%) / GC_Deviation |
| 污染预警 | Contam_Warning（Normal 或提示信息） |

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 直接看每列"好不好"即可：比对率、proper pair 比例高，重复率、软剪切率低，深度适中且均匀，就是好样本。

- `Map_Rate(%)` < 70% 会触发污染预警；说明样品与参考不匹配或混入了外源 DNA
- `Duplicates_Rate(%)` 过高（如 >30%）提示 PCR 过度或测序过深
- `SoftClip_Rate(%)` 高提示结构变异或接头/污染
- `Avg_Depth` 与 `Median_Depth` 差距大说明深度不均
- `Contam_Warning` 列：`Normal` 表示无预警；出现 `LowMapRate(...)` 或某染色体名则提示该样本/该染色体可疑
- `result.per_chromosome.tsv`：看某条染色体 `Coverage_Depth` 是否异常偏低/偏高，可用于定位污染或倍性异常

## 参数选择建议 | Parameter Guidance

- **标准全基因组质控**：默认参数 + `-g reference.fa`（补全 GC 偏差），全量六类统计
- **快速只看出库质量**：加 `--skip-sequence --skip-insert --skip-duplicate --skip-variation`，只看比对与覆盖度
- **只关心目标区域（如外显子组）**：加 `--bed-file`
- **样本很多**：调大 `-p`（如 32），保持 `-t 12`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | BAM文件或目录｜BAM file or directory |
| `--output, -o` | `bam_stats.summary.tsv` |  | 输出文件｜Output file |
| `--threads, -t` | `12` | int | samtools线程数｜Samtools threads |
| `--processes, -p` | `16` | int | 并行样本数｜Max parallel samples |
| `--reference, -g` | — |  | 参考基因组｜Reference genome FASTA |
| `--bed-file` | — |  | 目标区域BED｜Target regions BED file |
| `--min-mapq` | `20` | int | 最小MAPQ｜Minimum MAPQ |
| `--max-insert` | `1000` | int | 最大插入片段｜Maximum insert size |
| `--skip-alignment` | — |  | 跳过比对统计｜Skip alignment stats |
| `--skip-coverage` | — |  | 跳过覆盖度统计｜Skip coverage stats |
| `--skip-sequence` | — |  | 跳过序列特征｜Skip sequence features |
| `--skip-insert` | — |  | 跳过插入片段｜Skip insert size stats |
| `--skip-duplicate` | — |  | 跳过重复统计｜Skip duplicate stats |
| `--skip-variation` | — |  | 跳过变异统计｜Skip variation stats |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | BAM文件或目录｜BAM file or directory containing BAM files |
| `-o, --output` | `bam_stats.summary.tsv` |  | 输出文件｜Output file (default: bam_stats.summary.tsv) |
| `-t, --threads` | `12` | int | samtools线程数｜Samtools threads (default: 12) |
| `-p, --processes` | `16` | int | 并行处理的样本数｜Max parallel samples (default: 16) |
| `-g, --reference` | — |  | 参考基因组FASTA｜Reference genome FASTA (for GC bias) |
| `--bed-file` | — |  | 目标区域BED文件｜Target regions BED file |
| `--min-mapq` | `20` | int | 最小MAPQ阈值｜Minimum MAPQ threshold (default: 20) |
| `--max-insert` | `1000` | int | 最大插入片段｜Maximum insert size (default: 1000) |
| `--skip-alignment` | — | store_true | 跳过比对统计｜Skip alignment statistics |
| `--skip-coverage` | — | store_true | 跳过覆盖度统计｜Skip coverage statistics |
| `--skip-sequence` | — | store_true | 跳过序列特征统计｜Skip sequence feature statistics |
| `--skip-insert` | — | store_true | 跳过插入片段统计｜Skip insert size statistics |
| `--skip-duplicate` | — | store_true | 跳过重复统计｜Skip duplicate statistics |
| `--skip-variation` | — | store_true | 跳过变异统计｜Skip variation statistics |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `samtools`、`bedtools`（通过路径管理自动定位，可用环境变量 `SAMTOOLS_PATH` / `BEDTOOLS_PATH` 覆盖）
- Python 库：`pandas`、`openpyxl`、`tqdm`（缺任一会在启动时提示）

## 常见问题 | FAQ

**Q1：会不会断点续传？**
不会。每次运行都会重新分析所有样本。BAM 索引会复用（存在则跳过创建），但统计本身不缓存。

**Q2：报"缺少Python库 pandas/openpyxl/tqdm"？**
装齐即可，如 `pip install pandas openpyxl tqdm` 或使用已含这些库的 conda 环境。

**Q3：GC_Bias 列为空？**
GC 偏差需要参考基因组，请加 `-g reference.fa`。不给参考则该类指标缺失（不影响其他统计）。

**Q4：--min-mapq 为什么默认是 20 而不是 0？**
为了默认就只统计"比对可靠"的读段，得到的比对率/深度更有意义。想看全部读段可显式 `--min-mapq 0`。

**Q5：单 BAM 和目录都支持吗？**
都支持。单文件直接给路径，目录自动收集所有 BAM；也接受 SAM 文件。
