# GenomeSyn 基因组共线性分析 | Genome Synteny Analysis

一句话理解：**把两个或多个基因组的染色体「排排站」，自动比对出它们之间哪些段落是「同源、顺序一致」的（共线性），并画成一张直观的对比图**，让你一眼看出物种/个体间大片段的重排、倒位和缺失。

## 功能概述 | Overview

- 一条命令完成「比对 + 绘图」两步：用 minimap2（或 MUMmer）做基因组两两比对，再用 **NGenomeSyn** 画出共线性图
- 三种输入方式任选：直接给参考+查询 FASTA、给一张样本映射表、或给一份可编辑的配置文件（.xlsx/.yaml）
- 支持只分析指定染色体（如 1,2,3 或 1-5），画布尺寸、输出格式（SVG/PNG）可调
- 比对环节有断点续传：已存在且非空的比对文件会自动跳过（见 FAQ）
- 额外生成完整流程脚本与环境检查脚本，方便复现和排查

## 快速开始 | Quick Start

```bash
biopytools genomesyn -i ref.fa -I query.fa -o output
```

最小输入：一个参考基因组 FASTA + 一个或多个查询基因组 FASTA（-I 可重复），自动按文件名推断基因组名并完成比对绘图。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 共线性 | 两个基因组里「同源基因排列顺序一致」的段落；顺序乱了就说明发生过重排/倒位 |
| 比对 | 把两条序列逐段「对齐找相同」，输出它们在哪里匹配、匹配多长 |
| minimap2 | 一种极快的序列比对软件，专为长序列（整条染色体）设计，这里默认用它 |
| MUMmer | 另一种经典比对工具，对高度相似、无重复的大基因组更精细，速度较慢 |
| PAF | minimap2 的比对结果格式，一行一条比对（谁对谁、区间、方向、长度） |
| LINK | NGenomeSyn 绘图用的比对输入格式，由 PAF 转换而来 |
| NGenomeSyn | 画共线性图的软件：给它基因组长度文件（.len）+ 比对文件，它出 SVG/PNG 图 |
| sample map | 一张「文件路径 ↔ 基因组名」的两列对照表，让程序知道每个文件叫什么 |

## 输入 | Input

### 方式一：直接给 FASTA（-i + -I）

-i/--ref 给参考基因组 FASTA，-I/--query 给查询基因组（可重复多个）。程序按文件名（去扩展名）推断基因组名，生成一张临时 sample map 后继续。

### 方式二：样本映射表（--sample-map）

制表符分隔，两列：**基因组文件路径 + 基因组名**。

```text
/path/to/ref.fa<TAB>ref
/path/to/query1.fa<TAB>query1
/path/to/query2.fa<TAB>query2
```

- 以 # 开头的行和空行会被跳过；格式不对的行被警告并跳过
- 至少需要 2 个基因组

### 方式三：配置文件（--config）

.xlsx 或 .yaml。可先用 --generate-config + --sample-map 生成一份默认配置文件，编辑后再用 --config 跑。配置文件不能与 -i/-I 或 --sample-map 同时使用。

## 参数说明 | Parameters

### 输入选择 | Input selection

**通俗理解|In plain words:** 三种输入方式三选一，别混用。只想快速看两个基因组就用 -i/-I；要精细控制（指定画布、区域标注等）就用 sample map 或 config。

- -i, --ref：参考基因组 FASTA
- -I, --query：查询基因组 FASTA（可重复）
- --sample-map, -s：样本映射表（文件路径<TAB>基因组名）
- --config, -c：配置文件（.xlsx 或 .yaml），不能与其他输入参数同时用
- --generate-config：仅生成配置文件（配合 --sample-map），不跑分析

### 比对参数 | Alignment

**通俗理解|In plain words:** --aligner 选比对的「引擎」；--min-length 是「短于多长不算数」的过滤线，调大 = 只保留更可靠的共线性块、图更干净但可能漏掉小的；--no-minimap2-preset 关掉 minimap2 的预设参数，**一般不用动**。

- --aligner, -a：minimap2 / mcscanx / syri / mummer，默认 minimap2（实际实现 minimap2 与 mummer 两种）
- --alignment-mode：chain / star / all_vs_all，默认 chain（当前按样本顺序做相邻两两比对）
- --threads, -t：线程数，默认 12
- --min-length：最小比对长度，默认 5000 bp
- --no-minimap2-preset：不使用 minimap2 的 -x 预设参数

### 染色体过滤 | Chromosome filter

**通俗理解|In plain words:** 只想看某几条染色体时用 --chromosome（按 FASTA 里出现的顺序编号），格式 1,2,3 或 1-5。**一般不用，默认分析全部**。

- --chromosome：指定染色体，如 1,2,3 或 1-5

### 可视化参数 | Visualization

**通俗理解|In plain words:** 控制图多大、出什么格式、要不要在图上标出特殊区域。--canvas-width/height 不指定时按基因组数量自动算，**一般不用动**；--regions 用于把基因/区间高亮标在第一个基因组上。

- --canvas-width：画布宽度（不指定则自动计算）
- --canvas-height：画布高度（不指定则自动计算）
- --output-formats：svg / png，默认 svg png
- --regions：区域标注文件（特殊区域标注到图上）
- --output-dir, -o：输出目录，默认 ./genome_syn_output

## 分析流程 | Pipeline

```text
输入（-i/-I 或 sample map 或 config）
    |
    v
读样本 -> 生成配置文件（genome_syn_config.yaml/xlsx）
    |
    v
步骤1: 两两比对（minimap2 -> PAF -> GetTwoGenomeSyn.pl -> LINK）
    |      （已存在的非空比对文件自动跳过）
    v
步骤2: 生成 .len 文件 + ngenomesyn.conf，运行 NGenomeSyn 绘图
    |
    v
步骤3: 生成完整流程脚本 complete_genome_syn_pipeline.sh
```

## 输出 | Output

```text
output/
├── genome_syn_config.yaml            # 生成的配置（YAML）
├── genome_syn_config.xlsx            # 生成的配置（Excel，可编辑）
├── ref.len                           # 各基因组的染色体长度文件（文件名=基因组名）
├── query1_vs_ref_ref_vs_query1.paf   # minimap2 比对原始结果
├── query1_vs_ref_ref_vs_query1.link  # NGenomeSyn 用的比对文件
├── ngenomesyn.conf                   # NGenomeSyn 绘图配置
├── query1_vs_ref_genome_synteny.svg  # 最终共线性图（矢量）
├── query1_vs_ref_genome_synteny.png  # 最终共线性图（位图）
├── complete_genome_syn_pipeline.sh   # 完整流程脚本（可 bash 复跑）
├── check_environment.sh              # 环境检查脚本
└── *_commands.sh                     # 各步骤执行过的命令脚本
```

- 比对文件名格式：{查询}_vs_{参考}_{参考}_vs_{查询}.paf/.link（程序内部把 query 放前面）；若指定了染色体，文件名会加 _chr1_2 之类后缀
- 最终图名前缀随第一个比对对的 {查询}_vs_{参考} 变化，例如 query1_vs_ref_genome_synteny.svg

## 结果解读 | Interpreting Results

- **.svg / .png 图**：横排每个基因组的染色体（按 .len 顺序），连线表示共线性块。**连线大致平行、颜色分区规整 = 基因组间保守性高**；大量交叉、断裂的连线 = 发生了重排/倒位/缺失，正是要找的结构差异。
- **.link / .paf 文件**：比对原始数据，.link 每行记录一对同源区段（两端坐标 + 方向）。块数量多、覆盖长 = 两基因组相似度高。
- **.len 文件**：每行记录染色体名、起始、长度与填充色，是 NGenomeSyn 画染色体条的来源，也顺带给了每条染色体的长度。
- **complete_genome_syn_pipeline.sh**：完整 shell 脚本，可手动 bash 重跑整个流程（部分步骤需替换真实路径）。
- **好坏判据**：SVG 文件生成且非空即成功；若日志报「缺失比对文件」说明比对步骤没产出 .link/.paf，回头查比对日志。

## 参数选择建议 | Parameter Guidance

- **比对引擎**：一般基因组用默认 minimap2（快）；对高度相似、需要更精细锚点的近缘物种用 mummer（慢但稳）。mcscanx/syri 虽在选项里但当前实现未接入，建议仍用前两者。
- **--min-length**：默认 5000 对多数场景合适；图上噪点（碎片式短连线）太多可调到 10000 以上；想保留小规模共线性则降到 2000–3000。
- **--chromosome**：只关心个别染色体时用它，能显著减少比对时间和图复杂度。
- **--canvas-width/height**：先不指定让程序自动算；图太挤或太空再手动调。
- **--generate-config**：想精细控制比对关系/画布/标注时，先生成 config、编辑后 --config 重跑。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --ref` | — | Path | 参考基因组｜Reference genome file |
| `-I, --query` | — | Path | 比对基因组（可重复）｜Query genome file(s), can be repeated |
| `--sample-map, -s` | — |  | 样本映射文件｜Sample mapping file |
| `--config, -c` | — |  | 配置文件(.xlsx或.yaml)｜Configuration file (.xlsx or .yaml) |
| `--output-dir, -o` | `./genome_syn_output` | Path | 输出目录｜Output directory |
| `--generate-config` | — |  | 仅生成配置文件｜Generate configuration file only |
| `--aligner, -a` | `minimap2` | minimap2/mcscanx/syri/mummer | 比对工具｜Aligner type |
| `--alignment-mode` | `chain` | chain/star/all_vs_all | 比对模式｜Alignment mode |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-length` | `5000` | int | 最小比对长度｜Minimum alignment length |
| `--no-minimap2-preset` | `False` |  | 不使用minimap2的-x预设参数｜Disable minimap2 -x preset parameter |
| `--chromosome` | — | str | 指定染色体｜Specify chromosomes to analyze |
| `--canvas-width` | — | int | 画布宽度｜Canvas width |
| `--canvas-height` | — | int | 画布高度｜Canvas height |
| `--output-formats` | `['svg', 'png']` | svg/png | 输出格式｜Output formats |
| `--regions` | — | Path | 区域标注文件｜Special region annotation file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --sample-map` | — |  | 样本映射文件｜Sample mapping file (tab-separated: genome_file\tgenome_name) |
| `-c, --config` | — |  | 配置文件｜Configuration file (.xlsx or .yaml) |
| `-o, --output-dir` | `./genome_syn_output` |  | 输出目录｜Output directory (default: ./genome_syn_output) |
| `--generate-config` | — | store_true | 仅生成配置文件｜Generate configuration file only |
| `-a, --aligner` | `minimap2` | minimap2/mcscanx/syri/mummer | 比对器类型｜Aligner type (default: minimap2) |
| `--alignment-mode` | `chain` | chain/star/all_vs_all | 比对模式｜Alignment mode (default: chain) |
| `-t, --threads` | `16` | int | 线程数｜Number of threads (default: 16) |
| `--min-length` | `5000` | int | 最小比对长度｜Minimum alignment length (default: 5000) |
| `--no-minimap2-preset` | — | store_true | 不使用minimap2的-x预设参数｜Disable minimap2 -x preset parameter (default: False) |
| `--chromosome` | — | str | 指定要分析的染色体｜Specify chromosomes to analyze (e.g., "1,2,3" or "1-5" or "1") |
| `--canvas-width` | — | int | 画布宽度｜Canvas width (auto-calculated if not specified) |
| `--canvas-height` | — | int | 画布高度｜Canvas height (auto-calculated if not specified) |
| `--output-formats` | `['svg', 'png']` | svg/png | 输出格式｜Output formats (default: svg png) |
| `--regions` | — | str | 区域标注文件｜Special region annotation file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **NGenomeSyn**（必需，需在 PATH 中；绘图用）
- **minimap2**（默认比对器）或 **nucmer + show-coords**（MUMmer 比对器）
- **GetTwoGenomeSyn.pl**（PAF 转 LINK，随 NGenomeSyn 提供；缺失则退回使用 PAF）
- **ImageMagick convert**（仅在输出 PNG 时需要）
- **syri**（仅结构变异分析需要，当前流程默认不调用）
- Python 3 + Biopython、pandas、PyYAML、openpyxl（封装脚本自身）

## 常见问题 | FAQ

**Q1：换参数重跑，比对没重算？**
比对环节的断点续传按「.link 或 .paf 文件存在且非空」判断，命中即跳过。若改了比对参数（如 --min-length、--aligner）想重算，先删除输出目录里对应的 .paf/.link 文件。

**Q2：日志提示「只找到 PAF 文件，NGenomeSyn 可能需要 LINK 格式」？**
说明 GetTwoGenomeSyn.pl 不在 PATH 中，PAF 没转成 LINK。可把 NGenomeSyn 自带脚本目录加入 PATH，或确认比对结果本身正常（PAF 非空）。

**Q3：图上的连线方向/前缀为什么和我想的不一样？**
程序统一把「查询基因组」作为图的第一个、比对文件名前缀用 {查询}_vs_{参考}。图的方向由样本在 sample map 中的顺序决定，想调方向就调 sample map 顺序或改用 config。

**Q4：报「NGenomeSyn 未安装」？**
NGenomeSyn 必须能在 which NGenomeSyn 里找到。把它所在 bin/ 加入 PATH，或确认安装路径后再跑。

**Q5：PNG 没生成？**
PNG 由 ImageMagick 的 convert 从 SVG 转出。缺 convert 时只会有 SVG（日志会有警告）。可只保留 SVG（矢量、可无限放大），或装 ImageMagick 后再转。
