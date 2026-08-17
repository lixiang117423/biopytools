# 端粒识别 | Telomere Finder (tidk)

一句话理解：**自动找出基因组每条染色体「两端戴的帽子」（端粒重复序列）在哪、有多少，并画成分布图，用来判断组装是否真正「封顶」**。

## 功能概述 | Overview

- 基于 tidk（Telomere Identification toolKit）构建
- 五种分析模式：explore（探索）、find（按分类群查找）、search（自定义序列搜索）、plot（绘图）、pipeline（一键全流程，默认）
- pipeline 模式自动「探索 → 智能搜索 → 绘图」三连，无需知道该物种的端粒序列
- 内置分类群端粒数据库（昆虫/脊椎动物/植物等 50+ 类群），find 模式按分类群查
- 输出每窗口端粒重复次数统计（TSV）+ 端粒分布 SVG 图

## 快速开始 | Quick Start

```bash
biopytools find-telomere -g genome.fa -o results
```

默认走 pipeline 模式，自动完成端粒序列探索、定位和绘图，无需指定端粒序列。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 端粒(telomere) | 染色体两端「戴的帽子」，由一小段序列反复重复组成，保护染色体不被磨损 |
| 端粒重复序列 | 帽子上的「花纹」——如人/脊椎动物是 TTAGGG，植物多是 TTTAGGG，昆虫多是 TTAGG |
| 分类群(clade) | 生物学上的「家族门类」，如 Mammalia（哺乳纲）、Coleoptera（鞘翅目） |
| 窗口(window) | 把染色体切成一段段固定长度（默认 10000 bp）的小格，逐格数端粒重复次数 |
| 反向互补 | DNA 两条链互为「镜像翻译」；TTAGGG 的反向互补是 CCCTAA，搜索时两者常一起考虑 |
| SVG | 矢量图格式，放大不糊，适合看整条染色体的端粒分布 |

## 输入 | Input

一个基因组 FASTA 文件（.fa / .fasta / .fa.gz）：

```text
>Chr1
TTAGGGTTAGGGTTAGGG...ACGT...CCCTAACCCTAA
>Chr2
TTAGGGTTAGGG...
```

plot 模式还需要一个由 find/search 生成的窗口 TSV 文件（通过 -t/--tsv 传入）。

## 参数说明 | Parameters

### 分析模式 | Analysis mode

**通俗理解|In plain words:** -m 决定「这次只做哪一步」。默认 pipeline 一条龙；只想单独跑某一步时才切到对应模式。find 模式必须配 -c（分类群），search 模式必须配 -s（搜索序列），plot 模式必须配 -t（TSV 文件）。

### 输出参数 | Output parameters

**通俗理解|In plain words:** -o 是结果放哪，-p 是输出文件名的前缀（默认 telomere，生成 telomere_xxx.tsv）。一般不用动。

### 软件路径 | Software path

**通俗理解|In plain words:** --tidk-path 是 tidk 可执行文件的位置，默认指向 asm conda 环境。装好了就别动它。

### Explore 模式参数 | Explore mode parameters

**通俗理解|In plain words:** explore 是「在基因组里地毯式搜索，哪些短重复序列反复出现」。--explore-min/--explore-max 划定「找多长的重复单元」（默认 5-12 bp，端粒单元通常是 6 bp）；--explore-threshold 是「重复多少次才叫候选」；--explore-distance 是「离染色体末端多近才算端粒」。**全部用默认值即可，几乎不用动。**

### Find 模式参数 | Find mode parameters

**通俗理解|In plain words:** find 是「按物种家族查现成的端粒序列」。-c 告诉程序你的物种属于哪个类群（如 Mammalia、Aves、Poales）；-w 是切窗口的大小。窗口越大越省事、越小越精细，默认 10000 bp 通用。

### Search 模式参数 | Search mode parameters

**通俗理解|In plain words:** search 是「拿你指定的那条序列去基因组里搜」。-s 填端粒单元（如 TTTAGGG）；--format 决定输出 tsv 还是 bedgraph（给别的工具用的格式）。

### Plot 模式参数 | Plot mode parameters

**通俗理解|In plain words:** plot 是「把 find/search 的结果画成图」。-t 传窗口 TSV；--plot-height/--plot-width 控制图的大小，只在图太挤或太小时才调。

### 日志与其他 | Logging & others

**通俗理解|In plain words:** -v 打开详细日志；--log-file 把日志写到指定文件；--print-clades 打印程序内置支持的所有分类群及对应端粒序列（查表用，不跑分析）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** pipeline 模式像「先广撒网找出候选帽子花纹，再按优先级依次去验证哪个花纹真的出现在染色体两端，最后把最靠谱的结果画出来」。

```text
基因组 FASTA
    │
    ▼
步骤1/3: explore 探索端粒重复序列 → {prefix}_explore.tsv
    │
    ▼
步骤2/3: 智能搜索（按优先级依次尝试）
    ├─ 1. TTTAGGG（典型植物端粒）
    ├─ 2. 其他植物端粒（最多 5 个）
    ├─ 3. 探索结果的第一条候选
    └─ 首个「>=10 个端粒窗口」的序列即采纳
    │
    ▼
步骤3/3: plot 绘制端粒分布图 → {prefix}.svg
```

## 输出 | Output

```text
results/
├── telomere_explore.tsv                       # explore: 候选端粒重复序列
├── telomere_telomeric_repeat_windows.tsv      # 每窗口端粒重复次数（核心结果）
├── telomere_telomeric_repeat_counts.tsv       # 每条染色体的端粒重复总数
└── telomere.svg                              # 端粒分布图
```

注：文件名前缀随 -p/--prefix 改变，上面是默认前缀 telomere 的结果。

## 结果解读 | Interpreting Results

### 1. 端粒窗口文件（telomere_telomeric_repeat_windows.tsv）

**通俗理解|In plain words:** 这是核心结果——把每条染色体切成窗口，数每个窗口里「帽子花纹」出现了多少次。看染色体两端是否都有高计数的窗口。

- 第一列通常是染色体/序列名，后面是该窗口内端粒重复的出现次数；
- **好的组装**：每条染色体的两端各有一小段高计数区，中间几乎为 0；
- **端粒缺失**：某条染色体一端没有高计数区，说明这一端可能没组装完整。

### 2. 端粒分布图（telomere.svg）

**通俗理解|In plain words:** 把上面的数字画成图，横轴是染色体位置，纵轴是端粒重复次数。理想的图是「每条染色体两端各有一个尖峰」。

### 3. 候选端粒序列（telomere_explore.tsv）

**通俗理解|In plain words:** 探索步骤找到的「疑似帽子花纹」清单，告诉你这个物种的端粒重复序列大概长什么样。

## 参数选择建议 | Parameter Guidance

- **-m, --mode**：不确定就用默认 pipeline；已知物种端粒序列想省时间，用 find + 对应 -c；想精确控制用 search + -s
- **-c, --clade**：动物多为 Mammalia/Aves（TTAGGG），植物多为 Poales/Rosales 等（TTTAGGG），昆虫多为 TTAGG；先用 --print-clades 查表
- **-w, --window**：默认 10000 通用；端粒很短或想更精细看边界时可调小到 1000-5000
- **--explore-*** 参数：默认值即可，几乎不用动

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | — |  | 基因组文件｜Genome sequence file (FASTA format) |
| `--mode, -m` | `pipeline` | explore/find/search/plot/pipeline | 分析模式｜Analysis mode |
| `--output-dir, -o` | `./telomere_output` | Path | 输出目录｜Output directory |
| `--prefix, -p` | `telomere` |  | 输出前缀｜Output prefix |
| `--clade, -c` | — |  | 分类群名称｜Clade name |
| `--window, -w` | `10000` | int | 窗口大小｜Window size |
| `--search-string, -s` | — |  | 搜索字符串｜Search string |
| `--format` | `tsv` | tsv/bedgraph | 输出格式｜Output format |
| `--tsv, -t` | — |  | TSV文件｜TSV file |
| `--plot-height` | `200` | int | 图像高度｜Plot height |
| `--plot-width` | `1000` | int | 图像宽度｜Plot width |
| `--tidk-path` | `~/miniforge3/envs/asm/bin/tidk` |  | tidk软件路径｜tidk software path |
| `--explore-min` | `5` | int | 最小重复长度｜Min repeat length |
| `--explore-max` | `12` | int | 最大重复长度｜Max repeat length |
| `--explore-threshold` | `100` | int | 重复阈值｜Repeat threshold |
| `--explore-distance` | `0.01` | float | 距离比例｜Distance ratio |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--print-clades` | — |  | 打印支持的分类群列表｜Print supported clades list |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组序列文件｜Genome sequence file (FASTA format) |
| `-m, --mode` | `pipeline` | explore/find/search/plot/pipeline | 分析模式｜Analysis mode |
| `-o, --output-dir` | `./telomere_output` |  | 输出目录｜Output directory |
| `-p, --prefix` | `telomere` |  | 输出前缀｜Output prefix |
| `--explore-min` | `5` | int | 最小重复长度｜Minimum repeat length |
| `--explore-max` | `12` | int | 最大重复长度｜Maximum repeat length |
| `--explore-threshold` | `100` | int | 重复阈值｜Repeat threshold |
| `--explore-distance` | `0.01` | float | 距离比例｜Distance ratio |
| `-c, --clade` | — |  | 分类群名称｜Clade name |
| `-w, --window` | `10000` | int | 窗口大小｜Window size |
| `-s, --search-string` | — |  | 搜索字符串｜Search string |
| `--format` | `tsv` | tsv/bedgraph | 输出格式｜Output format |
| `-t, --tsv` | — |  | TSV文件路径｜TSV file path |
| `--plot-height` | `200` | int | 图像高度｜Plot height |
| `--plot-width` | `1000` | int | 图像宽度｜Plot width |
| `--plot-fontsize` | `12` | int | 字体大小｜Font size |
| `--plot-strokewidth` | `2` | int | 线条宽度｜Stroke width |
| `--tidk-path` | `~/miniforge3/envs/asm/bin/tidk` |  | tidk软件路径｜tidk software path |
| `-v, --verbose` | — | store_true | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--print-clades` | — | store_true | 打印支持的分类群列表｜Print supported clades list |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- tidk（默认路径 ~/miniforge3/envs/asm/bin/tidk，conda 环境 asm；程序会自动用 conda run -n asm --no-capture-output 包装调用）
- 无其他外部依赖

## 常见问题 | FAQ

**Q1：pipeline 模式会自动尝试哪些端粒序列？**
按优先级依次尝试：TTTAGGG（植物端粒）→ 其他植物端粒（最多 5 个）→ 探索结果的第一条候选。每试一个都搜一遍，首个能找到「>=10 个端粒窗口」的序列即被采纳并进入绘图。

**Q2：为什么 find 模式报错要 --clade？**
find 模式按分类群查内置数据库，必须用 -c 指定类群。先用 --print-clades 看支持哪些类群。

**Q3：search/plot 模式需要什么配套参数？**
search 必须配 -s（搜索序列）；plot 必须配 -t（一个由 find/search 生成的窗口 TSV 文件）。

**Q4：支持断点续传吗？**
不支持。每次运行都从头跑，重新执行即可覆盖旧输出。

**Q5：tidk 报找不到怎么办？**
检查 --tidk-path 是否指向正确的 tidk 可执行文件（默认在 asm conda 环境里）；路径不对时用 --tidk-path 指定。

