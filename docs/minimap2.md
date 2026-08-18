# minimap2 - 序列比对与未比对区间提取 | Minimap2 Alignment and Unmapped Region Extraction

一句话理解：**把两条（组）序列用 minimap2 做全基因组比对，找出「查询序列上有、但目标基因组上比对不上的区间」，把这些区间的坐标（BED）和序列（FASTA）提取出来**。常用于找两基因组间的差异/特有片段。

## 功能概述 | Overview

- 用 minimap2 对目标基因组与查询基因组做全基因组比对，输出 PAF 比对结果
- 从 PAF 中按 tp 类型（primary/secondary）和匹配长度筛选可靠比对
- 用「减法」策略找出查询序列上未被比对覆盖的区间（未比对区间）
- 输出未比对区间的 BED 坐标 + 用 seqkit 提取对应序列（FASTA，60bp 换行）
- 生成总结报告，统计未比对区间数量与总长度

## 快速开始 | Quick Start

```bash
biopytools minimap2 -t target_genome.fasta -q query_genome.fasta -o results/
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 目标基因组 | 作为「参考底板」的序列，比对到它上面 |
| 查询基因组 | 要被拿去比对的序列，找它上面目标基因组里没有的部分 |
| PAF | minimap2 的比对输出格式，逐行记录一段比对的对齐关系 |
| 未比对区间 | 查询序列上「目标基因组对不上」的连续片段，即两基因组差异处 |
| tp 类型 | 比对的主次标签：primary（主要比对）/ secondary（次要比对） |
| BED | 记录区间坐标（染色体、起、止）的标准格式，供提取/可视化 |
| seqkit | 一个序列处理工具，按 BED 坐标从 FASTA 里把序列切出来 |

## 输入 | Input

- 目标基因组 FASTA（-t/--target）：比对底板。
- 查询基因组 FASTA（-q/--query）：要分析的序列。

两者都支持普通 FASTA 文件，路径会自动转为绝对路径。输出文件名以查询基因组的文件名为基名自动生成。

## 参数说明 | Parameters

### 必需与输出参数 | Required and output

**通俗理解|In plain words:** -t 目标、-q 查询是两个必填输入；-o 输出目录（默认 ./minimap2_output）。输出文件名按查询基因组基名自动生成，一般不用管。

### 比对预设参数 | Preset

**通俗理解|In plain words:** -x 是 minimap2 的预设参数，决定「比对严不严、适合多近缘的序列」。asm5 适合高度相似的序列（如不同组装间比对，默认），asm10/asm20 依次放宽（序列差异更大用），map-ont/map-pb 是给纳米孔/PacBio 长读用的。默认 asm5 即可。

### 筛选参数 | Filtering

**通俗理解|In plain words:** -m 是最小匹配长度（比这个短的比对不算数，默认 1000bp），-u 是最小未比对区间长度（比这个短的区间忽略，默认 1000bp），--tp-type 决定保留哪种比对（P=主要比对默认，S=次要，SP=都保留）。调大 -m 更严格、只留长且可靠的比对；调小 -u 会保留更多细小的差异片段。

### 工具路径参数 | Tool paths

**通俗理解|In plain words:** -M 是 minimap2 可执行文件，-S 是 seqkit 可执行文件；默认自动解析功能域环境（align/misc），域环境缺失时回退 PATH 直接调用，装了别名或用其他路径时才需要显式指定。

## 分析流程 | Pipeline

```text
目标基因组 + 查询基因组
    │
    ▼
步骤1：minimap2 比对 → PAF 文件（已存在则跳过）
    │
    ▼
步骤2：读 PAF → 按 tp 类型 + 最小匹配长度筛选
    │
    ▼
步骤3：减法策略找出未比对区间（整条序列减去所有比对片段）
    │
    ▼
步骤4：写 BED 文件 → seqkit 按 BED 提取序列 → 60bp 换行 FASTA
    │
    ▼
步骤5：生成总结报告
```

## 输出 | Output

```text
results/
├── <查询基名>_alignment.paf    # minimap2 原始比对结果
├── <查询基名>_unmapped.bed     # 未比对区间坐标（BED 格式）
├── <查询基名>_unmapped.fa      # 提取的未比对区间序列（FASTA）
├── minimap2_summary.txt        # 总结报告
└── minimap2_analysis.log       # 运行日志
```

- <查询基名>_alignment.paf：完整比对结果，保留所有比对记录供复查。
- <查询基名>_unmapped.bed：三列（序列名、起始、终止，0 基半开区间），即查询序列上对不上的区域。
- <查询基名>_unmapped.fa：每个未比对区间一条序列，序列 ID 形如 序列名:起-止。
- minimap2_summary.txt：未比对区间数量、总长度、各序列统计。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** BED 文件告诉你「差异在哪」，FASTA 文件把差异片段的序列也给你，summary 汇总差异总量。

- 未比对区间数量与总长度：反映查询基因组相对目标基因组的「多出来/对不上」的总量，越大说明两基因组差异越大。
- 未比对区间 FASTA：可直接拿去注释或比对，确认这些片段是什么（如特有基因、重复序列、污染）。
- 若 summary 显示「未找到符合条件的未比对区间」，说明查询序列基本都被目标基因组覆盖（两基因组高度一致），或 -m 阈值太严把比对都过滤了。
- PAF 文件为空：提示 minimap2 没比对上任何东西，检查序列是否同源、预设是否合适。

## 参数选择建议 | Parameter Guidance

- 近缘基因组（同物种不同组装）：默认 asm5。
- 差异更大的物种：改用 asm10 或 asm20，放宽比对条件。
- 只想看大片差异：调大 -m 和 -u（如 5000），过滤琐碎短片段。
- 想保留更细的差异：调小 -u（如 500）。
- 长读数据（ONT/PacBio）比对：-x 用 map-ont/map-pb。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--target, -t` | 必填 | Path | 目标基因组文件路径｜Target genome file path |
| `--query, -q` | 必填 | Path | 查询基因组文件路径｜Query genome file path |
| `--output-dir, -o` | `./minimap2_output` | Path | 输出目录｜Output directory |
| `--preset, -x` | `asm5` | asm5/asm10/asm20/map-ont/map-pb | Minimap2预设参数｜Minimap2 preset parameters |
| `--threads, -p` | `12` | int | 线程数｜Number of threads |
| `--min-match, -m` | `1000` | int | 最小匹配长度阈值｜Minimum match length threshold |
| `--min-unmapped, -u` | `1000` | int | 最小未比对区间长度阈值｜Minimum unmapped region length threshold |
| `--tp-type` | `P` | S/P/SP | 保留的比对类型(tp类型)｜Alignment type to keep: S(secondary), P(primary), SP(both) |
| `--minimap2-path, -M` | `minimap2` | str | minimap2可执行文件路径｜minimap2 executable path |
| `--seqkit-path, -S` | `seqkit` | str | seqkit可执行文件路径｜seqkit executable path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-t, --target` | 必填 |  | 目标基因组文件路径｜Target genome file path |
| `-q, --query` | 必填 |  | 查询基因组文件路径｜Query genome file path |
| `-o, --output-dir` | `./minimap2_output` |  | 输出目录｜Output directory |
| `-x, --preset` | `asm5` | asm5/asm10/asm20/map-ont/map-pb | Minimap2预设参数｜Minimap2 preset parameters |
| `-p, --threads` | `8` | int | 线程数｜Number of threads |
| `-m, --min-match` | `1000` | int | 最小匹配长度阈值｜Minimum match length threshold |
| `-u, --min-unmapped` | `1000` | int | 最小未比对区间长度阈值｜Minimum unmapped region length threshold |
| `--tp-type` | `P` | S/P/SP | 保留的tp类型｜tp type to keep: S(secondary), P(primary), SP(both) - 默认P｜default P |
| `-M, --minimap2-path` | — |  | minimap2可执行文件路径(默认域环境自动解析)｜minimap2 executable path (default: auto domain env) |
| `-S, --seqkit-path` | — |  | seqkit可执行文件路径(默认域环境自动解析)｜seqkit executable path (default: auto domain env) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- minimap2（自动解析 align 域环境并经 conda run 调用；可用 -M 或环境变量 MINIMAP2_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- seqkit（自动解析 misc 域环境并经 conda run 调用；可用 -S 或环境变量 SEQKIT_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- Python 库：pandas（读 PAF）、biopython（seqkit 替代路径，实际提取由 seqkit 完成）

## 常见问题 | FAQ

**Q1：有没有断点续传？**
仅有「PAF 已存在则跳过比对」这一步的复用：若输出目录已有 <查询基名>_alignment.paf，重跑会直接跳过 minimap2 比对、从 PAF 处理开始。其余步骤无断点续传。

**Q2：PAF 文件为空是怎么回事？**
minimap2 没比对上任何片段。可能两序列差异过大、预设过严（如 asm5 对远缘序列），或输入文件格式有问题。换更宽松的 -x 预设再试。

**Q3：为什么筛完没有有效比对？**
-m 阈值太严会把短比对全过滤掉。可调小 -m 或确认序列确实存在同源比对。

**Q4：未比对区间的坐标是哪种？**
BED 文件是 0 基半开区间（[start, end)），可直接喂给 bedtools/IGV；FASTA 里的序列 ID 用的是 1 基闭区间（序列名:起-止），两者同一区间、只是坐标表示不同。

**Q5：minimap2/seqkit 找不到？**
默认调用系统 PATH 里的命令。若装了但不在 PATH，用 -M（minimap2）和 -S（seqkit）指定完整路径。
