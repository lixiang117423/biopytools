# Pan-Blocks 泛基因组共线性块构建 | Pan-Genome Syntenic Block Construction

一句话理解：把多个基因组的染色体两两比对，找出"大家都能对齐的同源大片段"，再按优先级去重，拼成一张每个基因组都能对上的"共线性块图"。

## 功能概述 | Overview

- 用 MUMmer(nucmer) 两两比对多个基因组，找染色体间的共线性区域
- 按优先级顺序迭代去重，把共线性区域划分成互不重叠的 Pan-Blocks
- 每块记录"这段是谁贡献的"，可追溯来源基因组
- 输出 BED 格式块表 + 可对齐的可视化图(共线性连线)
- 分步可独立运行(align/build/plot)或一键全跑(all)，支持断点续传
- `prepare` 子命令可从 FASTA 目录自动生成基因组列表文件

## 快速开始 | Quick Start

```bash
biopytools pan-blocks all -i genome_list.txt -o output_dir/
```

`-i` 是基因组列表文件（`名称<TAB>路径` 两列）；没有的话先用 `biopytools pan-blocks prepare -i fasta_dir/ -o genome_list.txt` 生成。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 共线性(synteny) | 两个基因组里"同一段祖传 DNA"的对应关系，像两版译文里对齐的段落 |
| Pan-Block | 一个"大家都对齐得上"的染色体片段区间，按来源基因组去重后互不重叠 |
| 两两比对(pairwise alignment) | 每次拿两个基因组互相比，找出它们之间所有对齐的段落 |
| delta-filter | 比对完做"只保留可靠、足够长的对齐"的过滤器 |
| contributor | 某个片段"最初来自哪个基因组"的标记，去重时决定归谁 |
| 优先级顺序(genome order) | 去重时谁先挑：排在前面的基因组优先"认领"重叠片段 |

## 输入 | Input

核心输入是**基因组列表文件**（制表符分隔两列：名称、FASTA 路径）：

```text
genomeA    /path/to/genomeA.fa
genomeB    /path/to/genomeB.fa
genomeC    /path/to/genomeC.fa
```

- 至少需要 2 个基因组；FASTA 支持 `.fa`、`.fasta`、`.fna` 及 `.gz` 压缩
- 每行 `#` 开头为注释，空行跳过
- 可选 `--genome-order` 文件(每行一个基因组名)指定去重优先级；不指定则按列表顺序
- 可选 `--chromosome` 只处理某一条染色体

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 给基因组列表，`-o` 给输出目录。分步运行时用 `--step align/build/plot` 只跑某一步，方便中途续跑。

### 比对与性能 | Alignment & performance

**通俗理解|In plain words:** `-t` 是总线程数，`--parallel-alignments` 是同时跑几对基因组比对(每对内部再分线程)；`--min-alignment-length` 是"多短的对齐直接扔掉"的下限，调大=只要大块、结果更干净也更快，调小=保留更多细碎对齐。**默认值一般够用。**

### 构图与可视化 | Layout & plotting

**通俗理解|In plain words:** `--genome-order` 决定去重优先级(排前面的先认领)；`--chromosome` 只画/只建指定染色体；`--plot-format` 选 svg 或 png。**一般不用动。**

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-list` | 必填 |  | 基因组列表文件｜Genome list file (name<TAB>path) |
| `-o, --output-dir` | `./pan_blocks_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--parallel-alignments` | `4` | int | 并行比对数｜Parallel alignments |
| `--min-alignment-length` | `10000` | int | 最小比对长度｜Min alignment length |
| `--genome-order` | — |  | 基因组优先级顺序文件｜Genome priority order file |
| `--chromosome` | — |  | 指定染色体｜Specific chromosome |
| `--plot-format` | `svg` | svg/png | 绘图格式｜Plot format |
| `-i, --input-dir` | 必填 |  | FASTA文件目录｜FASTA files directory |
| `-o, --output` | `./genome_list.txt` |  | 输出文件路径｜Output file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-list` | — |  | 基因组列表文件｜Genome list file (name<TAB>path) |
| `-o, --output-dir` | `./pan_blocks_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--parallel-alignments` | `4` | int | 并行比对数｜Number of parallel alignments (default: 4) |
| `--min-alignment-length` | `10000` | int | 最小比对长度｜Minimum alignment length for delta-filter (default: 10000) |
| `--genome-order` | — |  | 基因组优先级顺序文件｜Genome priority order file (one name per line) |
| `--chromosome` | — |  | 指定染色体｜Specific chromosome to process |
| `--step` | — | align/build/plot | 执行特定步骤｜Run specific step |
| `--nucmer` | — |  | nucmer可执行文件路径｜nucmer executable path |
| `--delta-filter` | — |  | delta-filter可执行文件路径｜delta-filter executable path |
| `--show-coords` | — |  | show-coords可执行文件路径｜show-coords executable path |
| `--bedtools` | — |  | bedtools可执行文件路径｜bedtools executable path |
| `--minimap2` | — |  | minimap2可执行文件路径｜minimap2 executable path |
| `--plot-format` | `svg` | svg/png | 绘图格式｜Plot format (default: svg) |
| `--plot-width` | `20` | int | 绘图宽度｜Plot width (default: 20) |
| `--plot-height` | `10` | int | 绘图高度｜Plot height (default: 10) |
| `--prepare` | — |  | 从目录自动生成genome_list.txt｜Auto-generate genome_list.txt from FASTA directory |
| `--prepare-output` | `./genome_list.txt` |  | prepare输出文件路径｜Output file path for prepare (default: ./genome_list.txt) |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先"两两对齐找共同段落"，再"按优先级去重切成不重叠的块"，最后"画成对齐图"。

```text
基因组列表文件
    │
    ▼
步骤1 align: nucmer 两两比对 → delta-filter 过滤 → show-coords 输出坐标
    │  (输出 01_coords/*.filtered.coords,已完成的比对自动跳过)
    ▼
步骤2 build: 收集共线性区域 → bedtools merge 合并 → 按顺序迭代 subtract 去重
    │  (输出 02_blocks/{chrom}.pan_blocks.bed,已完成的染色体自动跳过)
    ▼
步骤3 plot: 按 Pan-Blocks 画每个基因组的彩色块 + 相邻基因组共线性连线
    │  (输出 03_plots/{chrom}.svg 或 .png)
    ▼
00_pipeline_info/software_versions.yml + 99_logs/ 日志
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml               # 工具版本与运行参数
├── 01_coords/
│   ├── genomeA.vs.genomeB.filtered.coords  # 两两比对坐标(每对一个文件)
│   └── ...
├── 02_blocks/
│   ├── chr1.pan_blocks.bed                 # 每条染色体的 Pan-Blocks(1-based 闭区间)
│   └── ...
├── 03_plots/
│   ├── chr1.svg                            # 共线性块可视化图
│   └── ...
├── tmp/                                    # 临时文件(运行结束清理)
└── 99_logs/
    └── pan_blocks.log                      # 运行日志
```

关键文件说明：

- **`*.filtered.coords`**：show-coords 输出的比对坐标，是构建 Pan-Blocks 的原料，也可直接用于其它共线性分析
- **`{chrom}.pan_blocks.bed`**：四列 `Chr\tStart\tEnd\tGenome`，`Genome` 列即该片段被判定归属的"贡献者"基因组；坐标 1-based、两端闭区间
- **`{chrom}.svg`**：每行一个基因组的染色体条，彩色块=Pan-Block(颜色=贡献者)，相邻行之间的灰色/绿色连线=共线性(绿色表示方向倒置)

## 结果解读 | Interpreting Results

### 1. Pan-Blocks 表（`*.pan_blocks.bed`）

**通俗理解|In plain words:** 这是一张"分区表"——每条染色体被切成了哪些块，每块归谁。块的总数反映基因组间的共线性破碎程度。

- `Genome` 列是贡献者：排在 `--genome-order` 前面的基因组会优先认领重叠片段
- 块越大越连续，说明这些基因组结构越保守；块碎成很多小块，说明重排/变异多
- 想看某段区域的保守性，直接按坐标查这张表即可

### 2. 可视化图（`*.svg`）

**通俗理解|In plain words:** 直观看出"谁的哪一段和谁对齐"。颜色相同的块来自同一贡献者；相邻行之间的绿色连线表示方向倒置的共线性。

- 灰色细线=方向一致的共线性，绿色线=方向倒置(倒位等结构变异)
- 连线密集且连续=高度共线；连线稀疏断裂=这段经历了较多重排

### 3. 比对坐标（`*.filtered.coords`）

- 每行是一个比对片段，含 ref/query 的染色体与起止坐标
- 被 `--min-alignment-length`(默认 10000bp) 过滤过，短于该值的不在文件里

## 参数选择建议 | Parameter Guidance

- **`--min-alignment-length`**：物种近、基因组大(植物/哺乳)可调大(如 50000)减少碎片；物种远或想保留微共线性则调小
- **`--parallel-alignments`**：内存充足可调大加快比对；注意每对比对会再分线程，总并发 = parallel_alignments × (threads/parallel_alignments)
- **`--genome-order`**：把"参考/高质量"基因组放最前面，让它的区块优先被保留，结果更符合预期
- **`--chromosome`**：只想快速看某条关键染色体时指定，可大幅减少计算量
- **`--plot-format png`**：svg 适合矢量编辑，png 适合直接查看/插入报告

## 依赖 | Dependencies

- MUMmer 套件：`nucmer`、`delta-filter`、`show-coords`（conda 环境 `pan`，默认路径 `~/miniforge3/envs/pan/bin/`）
- `bedtools`（默认路径 `~/.local/bin/bedtools`）
- Python 3（matplotlib、pyyaml）
- `minimap2` 路径在配置中定义(默认 conda 环境 `align`)并记录到版本信息，但当前比对流程使用 MUMmer nucmer

## 常见问题 | FAQ

**Q1：中途断了，重跑会从头来吗？**
不会。align 步骤按 `*.filtered.coords` 是否存在跳过已完成比对；build 步骤按 `{chrom}.pan_blocks.bed` 是否存在跳过已完成染色体。直接用 `biopytools pan-blocks all` 重跑即可续跑。

**Q2：想只重跑某一步怎么做？**
用 `--step`：`biopytools pan-blocks plot -i genome_list.txt -o output_dir/` 只重画图，不重算比对和建块。

**Q3：换了 `--min-alignment-length` 为什么比对没变？**
比对结果按 coords 文件存在性续跑。换过滤阈值后要先删 `01_coords/` 里的旧 coords 文件，再重跑 align。

**Q4：报了"需要至少 2 个基因组"？**
检查 genome_list.txt 是否至少两行有效记录，且每行是"名称<TAB>路径"两列（用制表符分隔，不是空格）。

**Q5：画图时某条染色体没有输出？**
plot 依赖 `02_blocks/{chrom}.pan_blocks.bed`；若 build 步骤未生成该染色体(或比对没覆盖到)，该染色体会被跳过。先确认比对坐标里是否有这条染色体。
