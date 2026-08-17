# NGenomeSyn 基因组共线性可视化 | NGenomeSyn Genome Synteny Visualization

一句话理解：**给一份「基因组文件 ↔ 名字」的清单，自动完成两两比对并调用 NGenomeSyn 画出共线性图**，还可用 SyRI 额外标出结构变异——是把比对到出图一条龙串起来的轻量流程。

## 功能概述 | Overview

- 一条命令完成「比对 + NGenomeSyn 绘图」：默认用 minimap2 比对，也可选 MUMmer
- 输入只需一张样本映射表（genome_file + genome_name），按顺序做相邻两两比对
- 支持只分析指定染色体（1,2,3 或 1-5），输出 SVG/PNG
- 可选 SyRI 结构变异分析（--use-syri），在共线性之外额外标出插入/缺失/倒位等 SV
- 比对环节有断点续传（已存在非空的比对文件跳过）；生成环境检查脚本与命令脚本便于复现

## 快速开始 | Quick Start

```bash
biopytools ngenomesyn -s samples.txt -o output_dir
```

最小输入：一张样本映射表 samples.txt（两列：基因组文件路径 + 基因组名）+ 一个输出目录。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 共线性 | 两个基因组里「同源基因排列顺序一致」的段落；乱了说明发生重排/倒位 |
| minimap2 | 极快的长序列比对软件，默认用它把两条染色体比到一起 |
| MUMmer | 另一种更精细的比对工具，对高度相似大基因组更稳、更慢 |
| PAF | minimap2 的比对结果（一行一条比对），可转成 LINK 供绘图 |
| NGenomeSyn | 画共线性图的软件：给它长度文件 + 比对文件，出 SVG/PNG |
| SyRI | 结构变异检测软件：专门找出两基因组间的插入/缺失/倒位等大差异 |
| sample map | 「文件路径 ↔ 基因组名」两列对照表，程序据此知道每个文件叫什么 |
| 染色体过滤 | 只挑某几条染色体参与分析，按 FASTA 里出现顺序编号 |

## 输入 | Input

### 样本映射表（-s/--sample-map）

制表符分隔，两列：**基因组文件路径 + 基因组名**。

```text
/path/to/genome1.fa<TAB>genomeA
/path/to/genome2.fa<TAB>genomeB
/path/to/genome3.fa<TAB>genomeC
```

- `#` 开头与空行跳过；格式不对的行警告并跳过
- 至少需要 2 个基因组；按顺序做相邻两两比对（genomeA vs genomeB、genomeB vs genomeC ...）

### 配置文件（-c/--config）

参数上支持传入，但当前实现里**配置文件模式尚未接通**，传入会报「暂不支持从配置文件读取模式」（见 FAQ）。实际请用 --sample-map。

## 参数说明 | Parameters

### 输入与输出 | Input & output

**通俗理解|In plain words:** -s 给样本清单，-o 给输出目录。二者一个定位输入、一个定位输出。

- -s, --sample-map：样本映射表（genome_file<TAB>genome_name）
- -c, --config：配置文件（当前暂未支持）
- -o, --output：输出目录（必填）

### 比对参数 | Alignment

**通俗理解|In plain words:** --aligner 选比对引擎（minimap2 快、mummer 稳）。--min-length 是「短于多长不算」的过滤线，调大 = 图更干净但可能漏小片段。--minimap-preset 与 --mummer-* 是各自引擎的调优参数，**一般不用动**。

- -a, --aligner：minimap2 / mummer，默认 minimap2
- -t, --threads：线程数，默认 12
- --min-length：最小比对长度，默认 5000 bp
- --minimap-preset：minimap2 预设（如 asm5），默认 asm5
- --mummer-match-type：MUMmer 匹配类型 mum / mumreference / maxmatch，默认 mumreference
- --mummer-min-match：MUMmer 最小匹配长度，默认 20

### 染色体过滤 | Chromosome filter

**通俗理解|In plain words:** 只想看某几条染色体时用 --chromosome（按 FASTA 里出现顺序编号），格式 1,2,3 或 1-5。**默认分析全部，一般不用**。

- --chromosome：指定染色体，如 1,2,3 或 1-5

### 可视化参数 | Visualization

**通俗理解|In plain words:** 控制出什么格式的图。SVG 是矢量（可无限放大），PNG 是位图（方便直接插入文档）。要 PNG 需要 ImageMagick 的 convert。

- --output-formats：svg / png，默认 svg png
- --ngenomesyn-bin：NGenomeSyn 二进制路径（默认自动在 PATH 与常见位置查找）

### 结构变异分析（SyRI）| Structural variation via SyRI

**通俗理解|In plain words:** 默认只画共线性。想额外标出两基因组间的结构变异（插入/缺失/倒位）就加 --use-syri；SyRI 不在 PATH 时用 --syri-bin 指路径。

- --use-syri：启用 SyRI 结构变异分析（默认关闭）
- --syri-bin：SyRI 二进制路径

## 分析流程 | Pipeline

```text
样本映射表 samples.txt
    |
    v
检查环境（NGenomeSyn + 比对器 [+ SyRI]）
    |
    v
步骤1: 相邻两两比对（minimap2 -> PAF -> GetTwoGenomeSyn.pl -> LINK）
    |      （已存在的非空比对文件自动跳过）
    v
步骤2: 生成 .len 文件 + ngenomesyn.conf，运行 NGenomeSyn 绘图
    |
    v
输出 genome_synteny.svg / .png + 各命令脚本
```

## 输出 | Output

```text
output_dir/
├── genomeA.len                       # 各基因组染色体长度文件（文件名=基因组名）
├── genomeA_vs_genomeB.paf            # minimap2 比对原始结果
├── genomeA_vs_genomeB.link           # NGenomeSyn 用的比对文件
├── ngenomesyn.conf                   # NGenomeSyn 绘图配置
├── genome_synteny.svg                # 最终共线性图（矢量）
├── genome_synteny.png                # 最终共线性图（位图，需 ImageMagick）
├── ngenomesyn_pipeline.log           # 运行日志
├── check_environment.sh              # 环境检查脚本
└── *_commands.sh                     # 各步骤执行过的命令脚本
```

- 若启用 SyRI，比对产物会变为 .sam（minimap2+SyRI）或 .coords（MUMmer+SyRI），并额外生成 SyRI 的结构变异结果
- 指定染色体时，文件名会加 _chr1_2 之类后缀

## 结果解读 | Interpreting Results

- **genome_synteny.svg / .png**：核心结果。横排各基因组的染色体，连线表示共线性块。**连线规整平行 = 保守；交叉/断裂 = 重排或倒位**。
- **.link / .paf 文件**：比对原始数据，块多、覆盖长 = 两基因组相似度高。
- **.len 文件**：每行一条染色体（名字、起始、长度、颜色），也是染色体长度一览。
- **SyRI 结果**（若启用）：在共线性图基础上标出结构变异类型与位置，变异越多说明两基因组结构差异越大。
- **好坏判据**：日志末尾「分析完成」、SVG 文件生成且非空即成功；日志报「缺失比对文件」说明比对步骤没产出，回头查比对日志。

## 参数选择建议 | Parameter Guidance

- **比对引擎**：一般用默认 minimap2（快）；近缘物种要更精细锚点用 mummer（慢）。
- **--min-length**：图里短碎连线太多调到 10000+；想保留小共线性降到 2000–3000。
- **--use-syri**：需要结构变异结论时才加；它会额外耗时，且 SyRI 必须可用（可用 --syri-bin 指路径）。
- **--chromosome**：只关心个别染色体时用它，显著减少时间与图复杂度。
- **--output-formats**：无 ImageMagick 时只保留 svg（矢量、可无限放大）即可。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--sample-map, -s` | — |  | 样本映射文件｜Sample mapping file (genome_file\tgenome_name) |
| `--config, -c` | — |  | 配置文件｜Configuration file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--aligner, -a` | `minimap2` | minimap2/mummer | 比对器类型｜Aligner type |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-length` | `5000` | int | 最小比对长度｜Minimum alignment length |
| `--minimap-preset` | `asm5` |  | Minimap2预设模式｜Minimap2 preset |
| `--mummer-match-type` | `mumreference` | mum/mumreference/maxmatch | MUMmer匹配类型｜MUMmer match type |
| `--mummer-min-match` | `20` | int | MUMmer最小匹配长度｜MUMmer min match |
| `--chromosome` | — | str | 指定染色体｜Specify chromosomes (e.g., "1,2,3" or "1-5") |
| `--output-formats` | `['svg', 'png']` | svg/png | 输出格式｜Output formats |
| `--ngenomesyn-bin` | — |  | NGenomeSyn二进制文件路径｜NGenomeSyn binary path |
| `--use-syri` | `False` |  | 使用SyRI进行结构变异分析｜Use SyRI for structural variation analysis |
| `--syri-bin` | — |  | SyRI二进制文件路径｜SyRI binary path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --sample-map` | — |  | [FILE] 样本映射文件｜Sample mapping file (genome_file\tgenome_name) |
| `-c, --config` | — |  | [FILE] 配置文件｜Configuration file |
| `-o, --output` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-a, --aligner` | `minimap2` | minimap2/mummer | [STR] 比对器类型 (默认: minimap2)｜Aligner type (default: minimap2) |
| `-t, --threads` | `16` | int | [INT] 线程数 (默认: 16)｜Number of threads (default: 16) |
| `--min-length` | `5000` | int | [INT] 最小比对长度 (默认: 5000)｜Minimum alignment length (default: 5000) |
| `--minimap-preset` | `asm5` |  | [STR] Minimap2预设模式 (默认: asm5)｜Minimap2 preset (default: asm5) |
| `--mummer-match-type` | `mumreference` | mum/mumreference/maxmatch | [STR] MUMmer匹配类型 (默认: mumreference)｜MUMmer match type (default: mumreference) |
| `--mummer-min-match` | `20` | int | [INT] MUMmer最小匹配长度 (默认: 20)｜MUMmer min match (default: 20) |
| `--chromosome` | — | str | [STR] 指定染色体 (如: "1,2,3" 或 "1-5")｜Specify chromosomes |
| `--output-formats` | `['svg', 'png']` | svg/png | [STR] 输出格式 (默认: svg png)｜Output formats (default: svg png) |
| `--ngenomesyn-bin` | — |  | [FILE] NGenomeSyn二进制文件路径｜NGenomeSyn binary path |
| `--use-syri` | — | store_true | [FLAG] 使用SyRI进行结构变异分析｜Use SyRI for structural variation analysis |
| `--syri-bin` | — |  | [FILE] SyRI二进制文件路径｜SyRI binary path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **NGenomeSyn**（必需，需在 PATH 或可用 --ngenomesyn-bin 指定）
- **minimap2**（默认比对器）或 **nucmer + show-coords**（MUMmer 比对器）
- **GetTwoGenomeSyn.pl**（PAF 转 LINK，随 NGenomeSyn 提供；缺失则退回用 PAF）
- **ImageMagick convert**（仅输出 PNG 需要）
- **syri**（仅 --use-syri 时必需）
- Python 3 + Biopython（封装脚本自身）

## 常见问题 | FAQ

**Q1：换参数重跑，比对没重算？**
比对环节按「.paf/.link 文件存在且非空」跳过。改了比对参数想重算，先删除输出目录里对应的 .paf/.link 文件。

**Q2：-c/--config 为什么不生效？**
当前实现里配置文件模式未接通，传入 --config 会报「暂不支持从配置文件读取模式」。请改用 --sample-map 提供样本清单。

**Q3：PNG 没生成？**
PNG 由 ImageMagick 的 convert 从 SVG 转出。缺 convert 时只出 SVG（日志有警告）。可只保留 SVG，或装 ImageMagick 后再转。

**Q4：想标结构变异，加 --use-syri 报「SyRI 未安装」？**
SyRI 必须在 PATH 中，或用 --syri-bin 指到实际路径。确认 SyRI 装好并加对路径。

**Q5：报「找不到 NGenomeSyn 二进制文件」？**
程序先查 PATH、再查常见安装位置。找不到时用 --ngenomesyn-bin 直接指定 NGenomeSyn 的完整路径。
