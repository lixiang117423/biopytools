# PlotSR 多基因组共线性可视化 | PlotSR Multi-genome Synteny Visualization

一句话理解：**自动跑 minimap2 + SyRI + PlotSR，把多个基因组的染色体排成一列，用彩色连线标出它们之间的同源区段、倒位、易位、重复等结构变异**。

## 功能概述 | Overview

- 全自动流程：minimap2 相邻基因组两两比对 → SyRI 结构注释 → PlotSR 绘图
- 支持 2 个及以上基因组，可指定名称与顺序（map 文件或 `-n`）
- 支持按染色体过滤（`-c`，支持数字或名称）
- 断点续传：染色体长度 / BAM / SyRI / 最终图已存在则跳过（`--force-run` 重跑）
- 输出 pdf / png / svg

## 快速开始 | Quick Start

```bash
biopytools plotsr -i genome1.fa -i genome2.fa -o output/
```

最小输入：2 个基因组 FASTA（可 `-i` 多次）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 结构变异(SV) | 大片段的重排，如倒位、易位、重复 |
| 同源(syntenic) | 两基因组中彼此对应的一段 |
| SyRI | 结构重排识别软件 |
| PlotSR | 把 SyRI 结果画成染色体重排图 |
| 倒位(inversion) | 一段序列方向反了 |
| 易位(translocation) | 一段序列跑到了别的染色体 |

## 输入 | Input

- `-i/--input`：基因组 FASTA，可多次使用；也可是"包含基因组的文件夹"或"map 文件"
- map 文件格式（两列，制表符分隔，可指定名称与顺序）：

```text
Col-0	/path/to/col0.fa
Ler	/path/to/ler.fa
```

## 参数说明 | Parameters

### 输入输出 | Input & Output

**通俗理解|In plain words:** `-i` 至少给 2 个基因组；`-o` 是输出目录。基因组的顺序很重要——**`-i` 的顺序（或 map 文件的顺序）就是图上染色体的排列顺序**，相邻基因组两两比对。`-n` 用逗号分隔给基因组起可读名（如 `Col-0,Ler,Cvi`），不给则自动从文件名提取。

### 比对参数 | Alignment

**通俗理解|In plain words:** `--minimap2-preset`（默认 asm5）按亲缘关系选：近缘 asm5、较远 asm20。`-t` 线程数。`-s/--min-sr-size`（默认 10000）是最小结构变异大小，小于它的 SV 不画——调大过滤杂音，调小保留更细微的重排。**一般用默认**。

### 可视化参数 | Visualization

**通俗理解|In plain words:** `--output-format`（pdf/png/svg）、`-f` 字体大小、`-d` DPI、`--space-ratio`（0.1-0.75，同源染色体间距）。`-v` 垂直排列染色体（默认水平）；`--itx` 只画染色体间互作。**默认 pdf 即可，图太挤时调 font-size 或 space-ratio**。

### 过滤参数 | Filtering

**通俗理解|In plain words:** `--nosyn`/`--noinv`/`--notr`/`--nodup` 分别关闭"同源/倒位/易位/重复"四类标注。**一般都不加**；图太杂乱时可用它们只保留想看的类型。

### 染色体过滤 | Chromosome filter

**通俗理解|In plain words:** `-c` 指定只显示哪些染色体，可多次使用，支持数字（如 `1` 表示第 1 条）或名称（如 `Chr1`）。基因组很大、只想看少数染色体时用。

### 流程控制 | Control

**通俗理解|In plain words:** `--skip-existing`（默认开）断点续传，`--force-run` 强制重跑。**换参数重跑时用 --force-run，日常重跑直接跑**。

## 分析流程 | Pipeline

```text
输入基因组(2+)
    │
    ▼
提取染色体长度(samtools faidx → {name}.chrlen)
    │
    ▼
minimap2 相邻基因组两两比对 → {A}_vs_{B}.bam(.bai)
    │
    ▼
SyRI 结构注释 → {A}_vs_{B}syri.out → 过滤 → syri.filtered.out
    │
    ▼
生成 plotsr/genomes.txt 配置
    │
    ▼
PlotSR 绘图 → plot.{format}
```

## 输出 | Output

```text
output/
├── alignment/{A}_vs_{B}.bam(.bai)         # minimap2 比对结果（相邻基因组）
├── syri/{A}_vs_{B}syri.out                # SyRI 原始结构注释
├── syri/{A}_vs_{B}syri.filtered.out       # 过滤注释行后的结果（PlotSR 输入）
├── plotsr/genomes.txt                     # PlotSR 基因组配置
├── plotsr/{name}.chrlen                   # 各基因组染色体长度
├── plot.pdf                               # 最终多基因组共线性图（格式随 --output-format）
└── plotsr.log                             # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 共线性图（`plot.pdf`）

**通俗理解|In plain words:** 核心结果。每列是一个基因组，染色体竖排，彩色连线连接相邻基因组间的对应区段。

- 连线大体平行、颜色成束：两基因组共线性好
- 连线交叉：该区段有倒位
- 连线跳到别的染色体：易位
- 某区段连线消失：可能有缺失或分化过大

### 2. SyRI 输出（`syri.filtered.out`）

**通俗理解|In plain words:** 结构注释明细，记录了每段是 SYN（同源）/INV（倒位）/TRANS（易位）/DUP（重复）等。供 PlotSR 使用，也便于程序化统计各类 SV 数量。

### 3. 好坏判据

- 图清晰、连线规整：比对与注释正常
- SyRI 失败（流程中断）：多为染色体链方向不一致，需先统一方向（如用 RectChr 校正）

## 参数选择建议 | Parameter Guidance

- 顺序：`-i` 或 map 文件的顺序就是展示顺序，按亲缘或进化关系排
- 命名：用 `-n` 给可读名，避免长路径做标签
- 亲缘：近缘 asm5、较远 asm20（`--minimap2-preset`）
- 简化：杂音多用 `-s` 调大最小 SV 大小，或 `--noinv` 等只留想看的类型
- 只画部分染色体：`-c 1 -c 2` 或 `-c Chr1`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入基因组FASTA文件（可多次使用）、包含基因组的文件夹或map文件｜Input genome FASTA files (multiple), folder containing genomes, or map file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-n, --names` | — |  | 基因组名称（逗号分隔）｜Genome names (comma-separated) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minimap2-preset` | `asm5` | asm5/asm10/asm20 | minimap2预设参数｜minimap2 preset |
| `-s, --min-sr-size` | `10000` | int | 最小结构变异大小｜Minimum structural variant size |
| `--output-format` | `pdf` | pdf/png/svg | 输出格式｜Output format |
| `-f, --font-size` | `6` | int | 字体大小｜Font size |
| `-d, --dpi` | `300` | int | 图片DPI｜Image DPI |
| `--space-ratio` | `0.7` | float | 同源染色体间距(0.1-0.75)｜Space for homologous chromosomes |
| `-v, --vertical` | `False` |  | 垂直排列染色体｜Plot vertical chromosomes |
| `--itx` | `False` |  | 染色体间交互模式｜Inter-chromosomal plotting mode |
| `--nosyn` | `False` |  | 不绘制同源区域｜Do not plot syntenic regions |
| `--noinv` | `False` |  | 不绘制倒位｜Do not plot inversions |
| `--notr` | `False` |  | 不绘制易位｜Do not plot translocations |
| `--nodup` | `False` |  | 不绘制重复｜Do not plot duplications |
| `-c, --chr` | — |  | 指定要显示的染色体（可多次使用，支持数字如1或名称如Chr1）｜Specify chromosomes to display (can be used multiple times, supports number like 1 or name like Chr1) |
| `--skip-existing` | `True` |  | 跳过已完成的步骤（默认启用）｜Skip completed steps (default: enabled) |
| `--force-run` | `False` |  | 强制重新运行所有步骤｜Force re-run all steps |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 | append | 输入基因组FASTA文件（可多次使用）或包含基因组的文件夹｜Input genome FASTA files (can be used multiple times) or folder containing genomes |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-n, --names` | — |  | 基因组名称（逗号分隔）｜Genome names (comma-separated, e.g., Col-0,Ler,Cvi) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads [default: 12] |
| `--minimap2-preset` | `asm5` | asm5/asm10/asm20 | minimap2预设参数｜minimap2 preset [default: asm5] |
| `-s, --min-sr-size` | `10000` | int | 最小结构变异大小｜Minimum structural variant size [default: 10000] |
| `--output-format` | `pdf` | pdf/png/svg | 输出格式｜Output format [default: pdf] |
| `-f, --font-size` | `6` | int | 字体大小｜Font size [default: 6] |
| `-d, --dpi` | `300` | int | 图片DPI｜Image DPI [default: 300] |
| `--space-ratio` | `0.7` | float | 同源染色体间距(0.1-0.75)｜Space for homologous chromosomes [default: 0.7] |
| `-v, --vertical` | — | store_true | 垂直排列染色体｜Plot vertical chromosomes |
| `--itx` | — | store_true | 染色体间交互模式｜Inter-chromosomal plotting mode |
| `--nosyn` | — | store_true | 不绘制同源区域｜Do not plot syntenic regions |
| `--noinv` | — | store_true | 不绘制倒位｜Do not plot inversions |
| `--notr` | — | store_true | 不绘制易位｜Do not plot translocations |
| `--nodup` | — | store_true | 不绘制重复｜Do not plot duplications |
| `-c, --chr` | — | append | 指定要显示的染色体（可多次使用，支持数字如1或名称如Chr1）｜Specify chromosomes to display (can be used multiple times, supports number like 1 or name like Chr1) |
| `--skip-existing` | `True` | store_true | 跳过已完成的步骤（默认启用）｜Skip completed steps (default: enabled) |
| `--force-run` | — | store_false | 强制重新运行所有步骤｜Force re-run all steps |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- minimap2（比对，align 域环境）
- samtools（建索引、提取染色体长度，align 域环境）
- syri（结构注释，可用环境变量 SYRI_PATH 覆盖）
- plotsr（绘图，可用环境变量 PLOTSR_PATH 覆盖）

> 工具路径自动解析功能域环境（`~/miniforge3/envs/<域>/bin/<工具>`）并经 conda run 调用；域环境缺失时回退 PATH（程序用 `which` 检查）。`minimap2 | samtools sort` 管道按规范以域环境二进制直调（管道中禁用 conda run）。

## 常见问题 | FAQ

**Q1：报"至少需要 2 个基因组"？**
`-i` 至少给 2 个 FASTA；用 map 文件或文件夹时确认里面真的解析出了 ≥2 个基因组。

**Q2：SyRI 失败？**
常见原因是两基因组染色体链方向不一致。需先用工具校正方向，或在日志里看 SyRI 的具体报错。

**Q3：换参数重跑结果没变？**
断点续传按各步产物判断。改了 `-s`/`--minimap2-preset` 等参数后要加 `--force-run`，否则复用旧产物。

**Q4：map 文件格式？**
两列制表符分隔：`名称\t路径`，空行和 `#` 开头行忽略。

**Q5：-c 里数字是什么意思？**
数字按该基因组 .fai 里的染色体顺序（1=第一条），名称则直接匹配染色体名。
