# MSA 可视化 | MSA Visualization (MAFFT + matplotlib)

一句话理解：**把多条 DNA/蛋白序列自动比对并对齐，画成"彩色方块矩阵"图，让保守位置和差异一眼可见**。

## 功能概述 | Overview

- 默认自动 MAFFT 比对 + 可视化；也可 `--skip-align` 直接用现成比对结果
- 自动检测 DNA / 蛋白并选用默认配色（Nucleotide / Zappo）
- 支持多种配色方案、区域截取、换行、consensus 显示、NJ 树排序
- 输出 png / jpg / svg / pdf

## 快速开始 | Quick Start

```bash
biopytools msaviz -i sequences.fa -o output.png
```

最小输入：一个 FASTA（多条序列），自动比对后出图。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| MSA(多序列比对) | 把多条序列对齐，让同源位置排在同一列 |
| 保守(conserved) | 某一列所有序列的碱基/氨基酸相同 |
| consensus | 每列出现最多的字符，代表"主流" |
| gap(-) | 为了对齐而插入的空位 |
| NJ 树 | 按序列相似度建的树，用于给序列排序 |
| 配色方案(color scheme) | 不同字符染不同颜色的规则 |

## 输入 | Input

- `-i/--infile`：FASTA 多序列文件；若未比对，程序先自动跑 MAFFT
- 若已比对好：加 `--skip-align` 直接可视化

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是输入序列，`-o` 是输出图文件（扩展名决定格式：png/jpg/svg/pdf）。输入可以是"没对齐的多条序列"（自动比对）或"已对齐的 MSA"（加 `--skip-align`）。

### 比对参数 | Alignment

**通俗理解|In plain words:** `--skip-align` 表示输入已是比对结果、跳过 MAFFT；`--mafft-path` 指定 MAFFT 路径；`--mafft-params` 是传给 MAFFT 的参数（默认 `--auto --preservecase`）；`--threads` 是 MAFFT 线程数；`--keep-alignment` 默认开，会把比对结果存成 `{输出名}.aligned.fa` 方便复用。**一般全默认即可**；序列特别多时调大 `--threads`。

### 格式参数 | Format

**通俗理解|In plain words:** `--format` 是输入 MSA 的文件格式（默认 fasta，也支持 phylip/clustal 等，仅 `--skip-align` 时相关）；`--color-scheme` 指定配色，不填则 DNA 用 Nucleotide、蛋白用 Zappo。**一般不用动**。

### 区域参数 | Region

**通俗理解|In plain words:** `--start`/`--end` 只展示比对的某一段（1-based）。序列太长、只想看局部时用，**默认显示全长**。

### 显示参数 | Display

**通俗理解|In plain words:** `--wrap-length`（默认 60）是每行显示多少字符，超过就换行——长序列一定要用它，否则图宽到没法看；`--wrap-space-size` 是换行间距。`--label-type` 控制标签用 ID 还是完整描述；`--show-label` 显示序列名；`--show-grid` 加网格线；`--show-count` 在右侧显示每行非 gap 字符数；`--show-consensus` 显示 consensus 条（`--consensus-color`/`--consensus-size` 控制其颜色与高度）；`--sort` 按 NJ 树把相似的序列排到一起；`--dpi` 是图像清晰度。**日常只动 wrap-length 和 sort 即可**。

## 分析流程 | Pipeline

```text
输入 FASTA
    │
    ▼
（可选）MAFFT 多序列比对 → .aligned.fa
    │
    ▼
载入比对 → 检测 DNA/蛋白 → 选配色
    │
    ▼
（可选）区域截取 / 换行 / NJ 树排序 / consensus
    │
    ▼
matplotlib 绘制彩色方块矩阵 → 保存图像
```

## 输出 | Output

```text
{outfile}                    # 可视化图（png/jpg/svg/pdf）
{outfile.stem}.aligned.fa    # 比对结果（--keep-alignment 默认保留，--no-keep-alignment 则删）
```

> 未保留比对结果时（`--no-keep-alignment`），中间比对文件写在输出目录 `tmp/` 下，用完即弃。

## 结果解读 | Interpreting Results

### 1. 序列矩阵

**通俗理解|In plain words:** 每一行是一条序列，每一列是一个比对位置，每个格子颜色代表该位置的字符种类。

- 整列颜色一致：该位置完全保守
- 列内颜色混杂：该位置变异较大
- 大片 `-`（gap）：该序列在此处有插入缺失

### 2. consensus 条（`--show-consensus` 时）

**通俗理解|In plain words:** 图下方的柱状条，柱越高表示该列序列一致性越高（绝大多数序列相同）。**连续的高柱 = 高度保守区**，低柱 = 高变区。

### 3. NJ 树排序（`--sort` 时）

**通俗理解|In plain words:** 相似的序列被排到相邻行，差异大的排到两端，便于肉眼分组。

## 参数选择建议 | Parameter Guidance

- 短序列（<100 bp）：默认即可
- 长序列：`--wrap-length 60` 换行（或 80/100），避免图过宽
- 想突出保守区：加 `--show-consensus`
- 想按亲缘排序：加 `--sort`
- 只画局部：`--start`/`--end` 截取
- 论文出图：`-o output.svg`（矢量）或 `-o output.pdf`，`--dpi 300`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--infile, -i` | 必填 |  | 输入序列文件或MSA文件｜Input sequences or MSA file |
| `--outfile, -o` | 必填 |  | 输出可视化文件｜Output visualization file (*.png｜*.jpg｜*.svg｜*.pdf) |
| `--skip-align` | — |  | 跳过MAFFT比对（输入已是比对结果）｜Skip MAFFT alignment (input is already aligned) |
| `--mafft-path` | `mafft` |  | MAFFT可执行文件路径｜MAFFT executable path (default: mafft) |
| `--mafft-params` | `--auto --preservecase` |  | MAFFT参数｜MAFFT parameters (default: --auto --preservecase) |
| `--threads` | `12` | int | MAFFT线程数｜MAFFT threads (default: 4) |
| `--keep-alignment/--no-keep-alignment` | `True` |  | 保留比对结果文件｜Keep alignment result file (default: True) |
| `--format` | `fasta` |  | MSA文件格式｜MSA file format (default: fasta) |
| `--color-scheme` | — |  | 颜色方案｜Color scheme (default: auto: DNA->Nucleotide / 蛋白｜protein->Zappo) |
| `--start` | `1` | int | 起始位置(1-based)｜Start position (1-based, default: 1) |
| `--end` | — | int | 结束位置(1-based)｜End position (1-based, default: MSA length) |
| `--wrap-length` | `60` | int | 换行长度｜Wrap length (default: 60) |
| `--wrap-space-size` | `3.0` | float | 换行间距｜Wrap space size (default: 3.0) |
| `--label-type` | `id` | id/description | 标签类型｜Label type (default: id) |
| `--show-label/--no-show-label` | `True` |  | 显示序列标签｜Show sequence labels (default: True) |
| `--show-grid` | — |  | 显示网格｜Show grid |
| `--show-count` | — |  | 显示字符统计｜Show sequence character count |
| `--show-consensus` | — |  | 显示consensus序列｜Show consensus sequence |
| `--consensus-color` | `#1f77b4` |  | Consensus颜色｜Consensus color (default: #1f77b4) |
| `--consensus-size` | `2.0` | float | Consensus大小｜Consensus size (default: 2.0) |
| `--sort` | — |  | 按NJ树排序｜Sort by NJ tree |
| `--dpi` | `300` | int | 图像DPI｜Figure DPI (default: 300) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --infile` | 必填 |  | 输入序列文件或MSA文件｜Input sequences or MSA file |
| `-o, --outfile` | 必填 |  | 输出可视化文件｜Output visualization file (*.png｜*.jpg｜*.svg｜*.pdf) |
| `--skip-align` | — | store_true | 跳过MAFFT比对（输入已是比对结果）｜Skip MAFFT alignment (input is already aligned) |
| `--mafft-path` | — |  | MAFFT可执行文件路径(默认域环境自动解析)｜MAFFT executable path (default: auto domain env) |
| `--mafft-params` | `--auto --preservecase` |  | MAFFT参数｜MAFFT parameters (default: --auto --preservecase) |
| `--threads` | `4` | int | MAFFT线程数｜MAFFT threads (default: 4) |
| `--keep-alignment` | `True` | store_true | 保留比对结果文件｜Keep alignment result file (default: True) |
| `--no-keep-alignment` | — | store_false | 不保留比对结果文件｜Do not keep alignment result file |
| `--format` | `fasta` |  | MSA文件格式｜MSA file format (default: fasta) |
| `--color-scheme` | — |  | 颜色方案｜Color scheme (default: auto: DNA->Nucleotide / 蛋白｜protein->Zappo) |
| `--start` | `1` | int | 起始位置(1-based)｜Start position (1-based, default: 1) |
| `--end` | — | int | 结束位置(1-based)｜End position (1-based, default: MSA length) |
| `--wrap-length` | `60` | int | 换行长度｜Wrap length (default: 60) |
| `--wrap-space-size` | `3.0` | float | 换行间距｜Wrap space size (default: 3.0) |
| `--label-type` | `id` | id/description | 标签类型｜Label type (default: id) |
| `--show-label` | `True` | store_true | 显示序列标签｜Show sequence labels (default: True) |
| `--no-show-label` | — | store_false | 不显示序列标签｜Do not show sequence labels |
| `--show-grid` | — | store_true | 显示网格｜Show grid |
| `--show-count` | — | store_true | 显示字符统计｜Show sequence character count |
| `--show-consensus` | — | store_true | 显示consensus序列｜Show consensus sequence |
| `--consensus-color` | `#1f77b4` |  | Consensus颜色｜Consensus color (default: #1f77b4) |
| `--consensus-size` | `2.0` | float | Consensus大小｜Consensus size (default: 2.0) |
| `--sort` | — | store_true | 按NJ树排序｜Sort by NJ tree |
| `--dpi` | `300` | int | 图像DPI｜Figure DPI (default: 300) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- matplotlib（绘图）
- BioPython（`Bio.AlignIO` 等，读比对）
- MAFFT（比对；默认路径 `mafft`，可用 `--mafft-path` 指定）

## 常见问题 | FAQ

**Q1：报 MAFFT 未找到？**
检查 mafft 是否在 PATH，或用 `--mafft-path` 指定路径；若输入已是比对结果，加 `--skip-align` 绕过。

**Q2：--skip-align 时报格式错误？**
输入必须是已比对的对齐文件，且 `--format` 要匹配实际格式（默认 fasta）。

**Q3：图太宽/太小？**
用 `--wrap-length` 控制换行（每行字符数），或 `--start`/`--end` 只画局部。

**Q4：start/end 越界报错？**
start/end 是 1-based 且在比对长度范围内（1 ≤ start ≤ end ≤ 比对长度）。

**Q5：DNA 和蛋白的配色怎么自动判定的？**
按序列中 ATGCUN 字符占比判断（>90% 判为核苷酸），也可用 `--color-scheme` 手动指定。
