# Samplot 结构变异可视化 | Samplot SV Visualization

一句话理解：**把某个结构变异(SV)位点附近，多条 BAM 里的"读段如何跨越断点"画成图，直观验证这个 SV 是真的还是比对假象**。

## 功能概述 | Overview

- `plot` 子命令：单区域绘图
- `vcf` 子命令：批量绘制 VCF 中的结构变异
- 支持 BAM / CRAM，多样本并排对比
- 支持长读长展示（`--long-read`）

## 快速开始 | Quick Start

```bash
biopytools samplot plot -b sample.bam -c chr1 -s 1000 -e 5000 -t DEL
```

最小输入：一个 BAM + 染色体名 + 起止位置。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 结构变异(SV) | 大片段删除、重复、倒位、易位 |
| BND | 易位（断点连接） |
| 覆盖深度(depth) | 某位置 read 数；删除区会显著下降 |
| 配对 read | 一对两端 read；跨度异常提示 SV |
| 长读长(long-read) | 几千到几万 bp 的读段，能直接跨越断点 |

## 输入 | Input

- `plot`：`-b` BAM/CRAM（可多个）、`-c` 染色体、`-s` 起始、`-e` 结束
- `vcf`：`--vcf` 变异文件 + `-b` BAM/CRAM（可多个）

> CRAM 必须配合 `-r` 参考 FASTA 才能解码。

## 参数说明 | Parameters

### plot 子命令 | plot

**通俗理解|In plain words:** `-b` 给一个或多个 BAM（多样本并排对比），`-c/-s/-e` 指定看哪一段，`-t` 标注 SV 类型（DEL/DUP/INV/BND）。`-o` 指定输出文件名，不填则自动生成；`--output-dir` 输出目录。`-d`（默认 1）是"最多画多少正常配对 read"——调大画更多 read、更密，调小只画代表性 read。`-w` 窗口大小（不填自动）、`-z` 异常判定阈值（默认 4 倍标准差）。`-r` 参考（CRAM 必须）。`-n` 给每个样本起标题（与 `-b` 顺序对应）。`--coverage-only` 只画覆盖度、`--same-yaxis-scales` 统一各样本 Y 轴。`--samplot-path` 是 samplot 可执行文件路径。

### vcf 子命令 | vcf

**通俗理解|In plain words:** 从 VCF 里批量挑出 SV 逐一画图。`--vcf` 输入、`-b` BAM、`-d` 输出目录（默认 `samplot-out`）、`-O` 输出格式（png/pdf/eps/jpg，默认 png）。过滤参数：`--min-bp`（默认 20，最小 SV 长度）、`--max-mb`（最大 SV 长度）、`--min-entries`/`--max-entries`（默认 6-10，只画样本数在这个范围的 SV）、`--sample-ids` 指定样本、`--plot-all` 画所有样本、`--min-call-rate`/`--max-hets` 按 call rate 和杂合数过滤、`--downsample` 下采样、`-t` 线程、`--gff3` 叠加基因注释。**VCF 里 SV 很多时务必用过滤参数控制数量**。

## 分析流程 | Pipeline

```text
（单步，无多阶段流程）

plot:  BAM + 区域 → samplot plot → 单张图
vcf:   VCF + BAM → 按过滤条件挑 SV → samplot vcf → 输出目录多张图
```

## 输出 | Output

```text
# plot 子命令
{output_dir}/{自动命名或-o指定}.png        # 单区域 SV 图
{output_dir}/samplot_plot_{时间戳}.log      # 运行日志

# vcf 子命令
samplot-out/*.png（或其他格式）             # 每个 SV 一张图
samplot-out/samplot_vcf_{时间戳}.log        # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 单区域图（plot）

**通俗理解|In plain words:** 从上到下是每个样本的 read 堆叠，横轴是基因组坐标，纵轴是覆盖深度。

- 删除(DEL)：中间区域深度明显下降，两侧 read 配对跨度异常
- 重复(DUP)：中间区域深度升高
- 倒位(INV)：read 方向在该区反转
- 易位(BND)：read 突然跳到另一条染色体

### 2. 批量图（vcf）

**通俗理解|In plain words:** 每个 SV 一张小图，用于批量筛查。结合 `--gff3` 注释能看出 SV 是否落在基因上。

### 3. 好坏判据

- 深度信号清晰、read 支撑充分：SV 可信
- 深度忽高忽低、无规律：可能是比对噪声或低覆盖
- 图太密：用 `-d` 限制 read 数，或 vcf 模式用过滤参数

## 参数选择建议 | Parameter Guidance

- 快速看一个位点：`plot` + 默认参数
- 多样本对比：`plot -b a.bam -b b.bam -n A -n B`
- 批量筛查：`vcf --min-bp 50 --max-entries 10` 控制输出数量
- CRAM 输入：必须加 `-r reference.fa`
- 论文出图：`-O pdf` 或调 `--dpi 300`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --bams` | 必填 |  | BAM/CRAM文件路径｜BAM/CRAM file paths |
| `-c, --chrom` | 必填 |  | 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | 起始位置｜Start position |
| `-e, --end` | 必填 | int | 结束位置｜End position |
| `-t, --sv-type` | — |  | SV类型(DEL/DUP/INV/BND)｜SV type |
| `-o, --output-file` | — |  | 输出文件名｜Output file name |
| `--output-dir` | `.` |  | 输出目录｜Output directory |
| `-r, --reference` | — |  | 参考基因组(CRAM必须)｜Reference genome (required for CRAM) |
| `-d, --max-depth` | `1` | int | 最大正常pair数｜Max normal pairs to plot |
| `-w, --window` | — | int | 窗口大小｜Window size |
| `-z, --z` | `4` | int | 标准差倍数｜Number of stdevs from mean |
| `-H, --plot-height` | — | int | 图高｜Plot height |
| `-W, --plot-width` | `8` | int | 图宽｜Plot width |
| `--dpi` | `300` | int | DPI｜Dots per inch |
| `--long-read` | `1000` | int | 长读长最小长度｜Min length for long-read |
| `--coverage-only` | `False` |  | 仅显示覆盖度｜Show only coverage |
| `--same-yaxis-scales` | `False` |  | 统一Y轴｜Use same Y-axis scales |
| `-n, --titles` | — |  | 样本标题｜Sample titles |
| `--samplot-path` | `~/miniforge3/envs/viz/bin/samplot` |  | samplot路径｜samplot binary path |
| `--vcf` | 必填 |  | VCF文件路径｜VCF file path |
| `-d, --output-dir` | `samplot-out` |  | 输出目录｜Output directory |
| `-O, --output-type` | `png` | png/pdf/eps/jpg | 输出格式｜Output format |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |
| `--downsample` | `1` | int | 下采样数｜Downsample count |
| `--min-bp` | `20` | int | 最小SV长度(bp)｜Min SV length in bp |
| `--max-mb` | — | int | 最大SV长度(MB)｜Max SV length in MB |
| `--sample-ids` | — |  | 样本ID列表｜Sample ID list |
| `--plot-all` | `False` |  | 绘制所有样本｜Plot all samples |
| `--min-call-rate` | — | float | 最小call rate｜Min call rate |
| `--max-hets` | — | int | 最大杂合数｜Max heterozygotes |
| `--min-entries` | `6` | int | 最小样本数｜Min entries to plot |
| `--max-entries` | `10` | int | 最大样本数｜Max entries to plot |
| `--gff3` | — |  | GFF3注释文件｜GFF3 annotation file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --bams` | 必填 |  | [FILE] BAM/CRAM文件路径｜BAM/CRAM file paths |
| `-c, --chrom` | 必填 |  | [STR] 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | [INT] 起始位置｜Start position |
| `-e, --end` | 必填 | int | [INT] 结束位置｜End position |
| `-t, --sv-type` | — |  | [STR] SV类型(DEL/DUP/INV/BND)｜SV type |
| `-o, --output-file` | — |  | [FILE] 输出文件名｜Output file name |
| `--output-dir` | `.` |  | [DIR] 输出目录｜Output directory |
| `-r, --reference` | — |  | [FILE] 参考基因组(CRAM必须)｜Reference genome (required for CRAM) |
| `-d, --max-depth` | `1` | int | [INT] 最大正常pair数｜Max normal pairs to plot |
| `-w, --window` | — | int | [INT] 窗口大小｜Window size |
| `-z, --z` | `4` | int | [INT] 标准差倍数｜Number of stdevs from mean |
| `-H, --plot-height` | — | int | [INT] 图高｜Plot height |
| `-W, --plot-width` | `8` | int | [INT] 图宽｜Plot width |
| `--dpi` | `300` | int | [INT] DPI｜Dots per inch |
| `--long-read` | `1000` | int | [INT] 长读长最小长度｜Min length for long-read |
| `--coverage-only` | — | store_true | [FLAG] 仅显示覆盖度｜Show only coverage |
| `--same-yaxis-scales` | — | store_true | [FLAG] 统一Y轴｜Use same Y-axis scales |
| `-n, --titles` | — |  | [STR] 样本标题列表｜Sample title list |
| `--samplot-path` | `~/miniforge3/envs/viz/bin/samplot` |  | [FILE] samplot可执行文件路径｜samplot binary path |
| `--vcf` | 必填 |  | [FILE] VCF文件路径｜VCF file path |
| `-d, --output-dir` | `samplot-out` |  | [DIR] 输出目录｜Output directory |
| `-O, --output-type` | `png` | png/pdf/eps/jpg | [STR] 输出格式｜Output format |
| `-t, --threads` | `1` | int | [INT] 线程数｜Number of threads |
| `--downsample` | `1` | int | [INT] 下采样数｜Downsample count |
| `--min-bp` | `20` | int | [INT] 最小SV长度(bp)｜Min SV length in bp |
| `--max-mb` | — | int | [INT] 最大SV长度(MB)｜Max SV length in MB |
| `--sample-ids` | — |  | [STR] 样本ID列表｜Sample ID list |
| `--plot-all` | — | store_true | [FLAG] 绘制所有样本｜Plot all samples |
| `--min-call-rate` | — | float | [FLOAT] 最小call rate｜Min call rate |
| `--max-hets` | — | int | [INT] 最大杂合数｜Max heterozygotes |
| `--min-entries` | `6` | int | [INT] 最小样本数｜Min entries to plot |
| `--max-entries` | `10` | int | [INT] 最大样本数｜Max entries to plot |
| `--gff3` | — |  | [FILE] GFF3注释文件｜GFF3 annotation file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- samplot（conda 环境 `viz`，默认 `~/miniforge3/envs/viz/bin/samplot`）

## 常见问题 | FAQ

**Q1：CRAM 报错？**
CRAM 必须提供 `-r` 参考 FASTA 才能解码。

**Q2：vcf 模式一下生成太多图？**
用 `--min-bp`/`--max-mb` 过滤 SV 长度、`--min-entries`/`--max-entries` 过滤样本数，必要时 `--sample-ids` 只留目标样本。

**Q3：输出文件名是什么？**
plot 模式不指定 `-o` 时由 samplot 自动命名；指定 `-o` 则用该名。vcf 模式输出到 `-d` 目录，每个 SV 一个文件。

**Q4：-n 标题和 -b 的顺序？**
`-n` 逐个对应 `-b` 给的 BAM 顺序，数量要一致。

**Q5：samplot 找不到？**
检查 `~/miniforge3/envs/viz/bin/samplot` 或用 `--samplot-path` 指定路径。
