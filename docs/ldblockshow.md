# 连锁不平衡热图分析 | LD Heatmap Analysis (LDBlockShow)

一句话理解：**把一段基因组里 SNP 两两之间「绑不绑定遗传」的关系画成一张热图**，让你一眼看出哪些变异是「抱团」的，常用来给 GWAS 显著位点做区域可视化。

## 功能概述 | Overview

- 封装 LDBlockShow + ShowLDSVG，从 VCF / Genotype / Plink 三种输入生成 LD 热图
- 支持 D' 和 R² 两种 LD 度量（默认 D'），也可两种都算
- 内置 5 种 LD block 检测方法（Gabriel、SolidSpine、BlockCut、FixBlock、NoBlock）
- 可在热图上叠加 GWAS P 值散点与 GFF3 基因注释
- 单区域 `-r` 或 BED 多区域批量 `-b` 出图
- 输出 SVG/PNG/PDF 图像，外加位点、block、配对 LD 值三份数据文件

## 快速开始 | Quick Start

```bash
biopytools ldblockshow -i variants.vcf.gz -o ld_result -r chr1:1000000-2000000
```

最小输入：一个 VCF（或 Genotype/Plink）+ 一个分析区域。输出落到 `ld_result/` 下，以区域命名（此处为 `chr1_1000000_2000000.*`）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 连锁不平衡(LD) | 染色体上离得近的 SNP 往往「一起遗传」，不独立；LD 就是衡量这种「绑定」有多强 |
| D' / R² | 两种给 LD 打分的尺子，取值都在 0（无关）到 1（完全绑定）之间；R² 更常用，D' 对罕见变异更敏感 |
| LD block | 一整段「几乎总是一起遗传」的区域，热图里像一块深色色块 |
| MAF | 「少数派等位基因」在群体里占多少；太低=这个位点区分度差 |
| 缺失率 | 多少样本在这个位点「没测出来」 |
| GWAS P 值 | 每个位点与性状关联的显著程度，越小越显著（图中取 -log10，越大越高） |

## 输入 | Input

### 变异输入（三选一）

- VCF（`-i`，支持 `.vcf` / `.vcf.gz`，最常用）
- SNP Genotype 格式（`--in-genotype`）
- Plink 前缀（`--in-plink`，`bed+bim+fam` 或 `ped+map`）

### 分析区域（`-r` 与 `-b` 二选一）

- `-r chr:start-end`：单个区域，**1-based 闭区间**（与 VCF 的 POS 对齐），如 `chr1:1000000-2000000`
- `-b regions.bed`：BED 文件（每行 `chrom start end [name]`），**0-based 半开区间**，等价于多个 `-r` 批量出图；有第 4 列 name 时用它做输出文件名前缀

### 可选注释文件

- GWAS（`--in-gwas`）：3 列 `chr pos pvalue`，染色体名须与 VCF 一致；P 值为 0/负/非数值的行会被自动丢弃
- GFF3（`--in-gff`）：标准 GFF3 注释，热图旁显示基因结构

## 分析流程 | Pipeline

```text
输入 VCF/Genotype/Plink + 区域
    │
    ▼
步骤1: LDBlockShow 计算两两 LD 值 + block 检测
    │   产出 .site.gz / .blocks.gz / .TriangleV.gz
    ▼
步骤2: ShowLDSVG 绘制 SVG（可选叠加 GWAS/GFF）
    │   产出 .svg / .png / .pdf
    ▼
输出
```

## 输出 | Output

```text
ld_result/
├── chr1_1000000_2000000.svg         # LD 热图（矢量，可无限放大）
├── chr1_1000000_2000000.png         # PNG 位图（默认输出）
├── chr1_1000000_2000000.pdf         # PDF（--out-pdf 时）
├── chr1_1000000_2000000.site.gz     # 过滤后保留的 SNP 位点列表
├── chr1_1000000_2000000.blocks.gz   # 检测到的 LD blocks
└── chr1_1000000_2000000.TriangleV.gz # 区域两两 SNP 的 LD 值（R²/D'）
```

## 结果解读 | Interpreting Results

- **热图颜色**：颜色越深（R²/D' 越接近 1）表示两个 SNP「绑定越紧」；越浅（接近 0）表示越独立。深色三角形/方形色块 = 一个 LD block
- **block 边界**：`blocks.gz` 里记录被判定为一个 block 的 SNP 范围，热图上通常有黑色边框标注
- **GWAS 轨道**：叠加 `--in-gwas` 后，热图上方/一侧出现散点，散点越高（-log10 P 越大）越显著；配合 `--cutline`（默认 5.0）画阈值线
- **好坏判据**：候选 GWAS 峰值位点若落在深色（高 LD）色块内，说明它附近有一整段「共分离」区域，功能解释更可信

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 绝大多数情况只用 `-i -o -r`（或 `-b`）三个参数就能出图；下面这些只有在特定需求下才动。

- **`--sele-var`（LD 度量）**：默认 1（D'）即可；多数文章用 R² 更有可比性，可设 2；想两种都看设 4
- **`--block-type`**：默认 1（Gabriel，最通用）一般不用动；只有明确要用自定义阈值才选 3（配 `--block-cut`），用预定义 block 文件选 4（配 `--fix-block`），不想画 block 选 5
- **SNP 过滤（MAF/Miss/HWE/Het）**：默认**关闭**——即默认不对任何 SNP 做过滤。注意：`biopytools ldblockshow` 命令行虽然暴露了 `--maf/--miss/--hwe/--het`，但 SNP 过滤的总开关 `--enable-snp-filter` 只在直接调用 `python -m biopytools.ldblockshow` 时可用，因此通过 `biopytools ldblockshow` 传这些值不会生效（详见 FAQ）
- **`--in-gwas` / `--in-gff`**：要给 GWAS 位点或基因结构做可视化时才加
- **`--out-pdf`**：需要 PDF 排版（论文插图）时加；PNG 默认已输出
- **`--no-show-ldist`**：区域内 SNP 极多、图太大时，给一个距离上限只显示近距离 LD，减小体积

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF变异文件路径｜VCF variant file path |
| `-o, --output-dir` | 必填 |  | 输出目录(自动创建)；每 region 产物落在 目录/<label>.* ｜Output directory (auto-created); per-region outputs land in dir/<label>.* |
| `-r, --region` | — |  | 分析区域，格式chr:start-end｜Analysis region, format chr:start-end |
| `-b, --bed` | — |  | 基因组BED文件(每行 chrom start end [name])，等价多个 -r 批量出图｜Genomic BED (cols: chrom start end [name]), equivalent to multiple -r |
| `--in-genotype` | — |  | SNP Genotype格式文件路径｜SNP Genotype format file path |
| `--in-plink` | — |  | Plink文件前缀(bed+bim+fam或ped+map)｜Plink file prefix |
| `--sele-var` | `1` | 1/2/3/4 | 选择LD度量统计量｜Select LD statistic (1=D', 2=R², 3/4=Both) |
| `--maf` | `0.05` | float | 最小次要等位基因频率｜Minimum minor allele frequency |
| `--miss` | `0.25` | float | 最大缺失率｜Maximum missing ratio |
| `--hwe` | `0.0` | float | Hardy-Weinberg平衡P值阈值｜Hardy-Weinberg equilibrium P-value threshold |
| `--het` | `1.0` | float | 最大杂合率｜Maximum heterozygosity ratio |
| `--enable-oth-var` | — |  | 允许indel/SV/CNV变异｜Allow bi-indel bi-sv bi-cnv variants |
| `--block-type` | `1` | 1/2/3/4/5 | Block检测方法 (1=Gabriel, 2=SolidSpine, 3=BlockCut, 4=FixBlock, 5=NoBlock) |
| `--block-cut` | `0.85:0.90` |  | BlockType3的cutoff｜Cutoff for BlockType3 |
| `--fix-block` | — |  | 固定block文件路径｜Fixed block file path |
| `--in-gwas` | — |  | GWAS文件，3列: chr pos pvalue，chr名须与VCF一致｜GWAS file, 3 cols: chr pos pvalue |
| `--in-gff` | — |  | GFF3注释文件路径｜GFF3 annotation file path |
| `--mer-min-snp-num` | `50` | int | 合并网格的最小SNP数｜Minimum SNP number to merge grids |
| `--no-out-png` | — |  | 不输出PNG格式图像｜Do not output PNG format image |
| `--out-pdf` | — |  | 输出PDF格式图像｜Output PDF format image |
| `--sub-pop` | — |  | 亚群样本文件路径｜Subgroup sample file path |
| `--tag-snp-cut` | `0.8` | float | TagSNP的LD cutoff｜LD cutoff for TagSNP |
| `--cutline` | `5.0` | float | GWAS P值显著性阈值(-log10)｜GWAS P-value significance cutoff (-log10) |
| `--point-size` | — | int | GWAS散点大小｜GWAS point size |
| `--top-site` | — |  | 指定GWAS峰值位点(chr:pos)｜Specify GWAS peak site (chr:pos) |
| `--no-log-p` | — |  | 不对P值取-log10｜Do not -log10 transform P-value |
| `--no-gene-name` | — |  | 不显示基因名｜Do not show gene names |
| `--show-num` | — |  | 在热图中显示R²/D'值｜Show R²/D' values in heatmap |
| `--spe-snp-name` | — |  | 特殊SNP名称文件(chr site Name)｜Special SNP name file |
| `--show-gwas-spe-snp` | — |  | 在GWAS图中显示特殊SNP名称｜Show special SNP names in GWAS plot |
| `--resize-h` | — | int | 图像高度，宽度按比例自动调整｜Image height, width auto-adjusted |
| `--no-show-ldist` | — | int | 超过此距离的SNP对不显示LD｜NoShow pairwise LD over this distance |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-o, --output-dir` | 必填 |  | 输出目录(自动创建)；每 region 产物落在 目录/<label>.* ｜Output directory (auto-created); per-region outputs land in dir/<label>.* |
| `-r, --region` | — |  | 单个分析区域，格式chr:start-end｜Single region, format chr:start-end |
| `-b, --bed` | — |  | 基因组BED文件(每行 chrom start end [name])，等价多个 -r 批量出图｜Genomic BED (cols: chrom start end [name]), equivalent to multiple -r |
| `-i, --vcf-file` | — |  | VCF变异文件路径｜VCF variant file path |
| `--in-genotype` | — |  | SNP Genotype格式文件路径｜SNP Genotype format file path |
| `--in-plink` | — |  | Plink文件前缀(bed+bim+fam或ped+map)｜Plink file prefix (bed+bim+fam or ped+map) |
| `--sele-var` | `1` | 1/2/3/4 | 选择LD度量统计量｜Select LD statistic: 1=D' (default), 2=R², 3/4=Both |
| `--enable-snp-filter` | `False` | store_true | 启用SNP过滤（默认关闭），使用-MAF/-Miss/-HWE/-Het参数过滤SNP｜Enable SNP filtering (default OFF), filter SNPs using -MAF/-Miss/-HWE/-Het parameters |
| `--maf` | `0.05` | float | 最小次要等位基因频率｜Minimum minor allele frequency (default: 0.05) |
| `--miss` | `0.25` | float | 最大缺失率｜Maximum missing ratio (default: 0.25) |
| `--hwe` | `0.0` | float | Hardy-Weinberg平衡P值阈值｜Hardy-Weinberg equilibrium P-value threshold (default: 0.0) |
| `--het` | `1.0` | float | 最大杂合率｜Maximum heterozygosity ratio (default: 1.0) |
| `--enable-oth-var` | — | store_true | 允许indel/SV/CNV变异｜Allow bi-indel bi-sv bi-cnv variants |
| `--block-type` | `1` | 1/2/3/4/5 | Block检测方法｜Block detection method: 1=Gabriel(default), 2=SolidSpine, 3=BlockCut, 4=FixBlock, 5=NoBlock |
| `--block-cut` | `0.85:0.90` |  | BlockType3的cutoff（格式：cutoff:ratio）｜Cutoff for BlockType3 (format: cutoff:ratio) |
| `--fix-block` | — |  | 固定block文件路径｜Fixed block file path (for BlockType=4) |
| `--in-gwas` | — |  | GWAS P值文件，制表符/空格分隔，3列: chr pos pvalue，染色体名须与VCF一致｜GWAS P-value file, tab/space delimited, 3 columns: chr pos pvalue, chr names must match VCF |
| `--in-gff` | — |  | GFF3注释文件路径｜GFF3 annotation file path |
| `--mer-min-snp-num` | `50` | int | 合并网格的最小SNP数｜Minimum SNP number to merge grids (default: 50) |
| `--no-out-png` | — | store_true | 不输出PNG格式图像｜Do not output PNG format image |
| `--out-pdf` | — | store_true | 输出PDF格式图像｜Output PDF format image |
| `--sub-pop` | — |  | 亚群样本文件路径｜Subgroup sample file path |
| `--tag-snp-cut` | `0.8` | float | TagSNP的LD cutoff｜LD cutoff for TagSNP (default: 0.80) |
| `--cutline` | `5.0` | float | GWAS P值显著性阈值(-log10)｜GWAS P-value significance cutoff (-log10) (default: 5) |
| `--point-size` | — | int | GWAS散点大小｜GWAS point size |
| `--top-site` | — |  | 指定GWAS峰值位点(chr:pos)｜Specify GWAS peak site (chr:pos) |
| `--no-log-p` | — | store_true | 不对P值取-log10｜Do not -log10 transform P-value |
| `--no-gene-name` | — | store_true | 不显示基因名｜Do not show gene names |
| `--show-num` | — | store_true | 在热图中显示R²/D'值｜Show R²/D' values in heatmap |
| `--spe-snp-name` | — |  | 特殊SNP名称文件(chr site Name)｜Special SNP name file (chr site Name) |
| `--show-gwas-spe-snp` | — | store_true | 在GWAS图中显示特殊SNP名称｜Show special SNP names in GWAS plot |
| `--resize-h` | — | int | 图像高度，宽度按比例自动调整｜Image height, width auto-adjusted |
| `--no-show-ldist` | — | int | 超过此距离的SNP对不显示LD｜NoShow pairwise LD over this distance |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- LDBlockShow（编译好的独立二进制，默认路径 `~/software/LDBlockShow/bin/LDBlockShow`，可用环境变量 `LDBLOCKSHOW_PATH` 或 `~/.config/biopytools/config.yml` 覆盖），需同目录附带 ShowLDSVG
- rsvg-convert（来自 librsvg，推荐安装；用于把大 SVG 可靠地渲染成 PNG，可用环境变量 `RSVG_CONVERT_PATH` 指定）
- ImageMagick `convert`（可选，PNG 渲染的兜底）

## 常见问题 | FAQ

**Q1：我设了 `--maf 0.01`，为什么结果没变化？**
SNP 过滤默认关闭（内部强制 `-MAF 0.0 -Miss 1.0`）。开关 `--enable-snp-filter` 目前只在直接调用 `python -m biopytools.ldblockshow` 时暴露，click 包装器 `biopytools ldblockshow` 未透传，所以这些过滤参数在包装器下不生效。

**Q2：换 `--sele-var` 或 `--block-type` 重跑，结果没变？**
LD 计算步骤按 `.site.gz` 是否存在做断点续传。改了 LD 度量或 block 方法后，需先删除旧的 `<label>.site.gz`（以及 `.blocks.gz`/`.TriangleV.gz`），否则会复用旧结果。

**Q3：PNG 是空白/打不开？**
ShowLDSVG 内置的 batik 对大 SVG（数百万 polygon）会渲染出空白 PNG。安装 librsvg（`conda install -c conda-forge librsvg`）或设 `RSVG_CONVERT_PATH` 指向 rsvg-convert；SVG 本身是正常的，也可手动转换。

**Q4：`-r` 和 `-b` 的坐标怎么不一样？**
`-r` 是 1-based 闭区间（和 VCF 的 POS 一致），BED 是 0-based 半开区间。同一段区域 BED 写 `chr1 999999 2000000`，等价于 `-r chr1:1000000-2000000`。
