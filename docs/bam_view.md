# BAM 比对可视化 | BAM Alignment Visualization (alignoth)

一句话理解：**把某个基因组区间内"每条测序读段怎么贴到参考序列上"画成一张图**，方便肉眼核对错配、插入缺失和变异位点。

## 功能概述 | Overview

- 用 alignoth 从 BAM + 参考 FASTA 生成指定区间的交互式比对视图
- 支持 html（内嵌 JS 可交互）/ json（Vega-Lite 规范）/ svg / pdf 四种输出
- 支持用 VCF、BED 或自定义区间高亮变异位点或目标区域
- 运行前校验 BAM 及索引（.bai）、VCF 索引（.csi/.tbi）、alignoth 是否存在

## 快速开始 | Quick Start

```bash
biopytools bam-view -b alignments.bam -r reference.fa -g chr1:1000-2000
```

最小输入：一个 BAM（需带 `.bai` 索引）+ 参考 FASTA + 一个区间（`chr:start-end`）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BAM | 压缩的比对结果文件，记录每条读段贴到参考的哪、怎么贴 |
| read（读段） | 测序仪读出的一段序列，像一条"碎纸片" |
| 比对(alignment) | 把 read 与参考序列"对位"的过程 |
| 错配(mismatch) | read 与参考在某个位置碱基不一致 |
| 覆盖深度(depth) | 某个位置叠了多少条 read，像"叠了多少层纸" |
| 索引(.bai) | BAM 的"目录"，让程序能快速定位到某个区间 |

## 输入 | Input

- `-b/--bam`：BAM 文件，**必须带同名 `.bai` 索引**（如 `alignments.bam.bai`）
- `-r/--reference`：参考序列 FASTA
- `-g/--region`：可视化区间，格式 `chr:start-end`（如 `chr1:1000-2000`）
- 可选：`-v/--vcf`（需 `.csi` 或 `.tbi` 索引）、`--bed`、`-H/--highlight`

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 三个参数缺一不可：`-b` 是要看的比对文件，`-r` 是它比对到的那条参考，`-g` 是"放大镜对准哪一段"。BAM 和参考必须是同一次比对的产物，否则图对不上。

### 软件配置 | Software

**通俗理解|In plain words:** `--alignoth-path` 是底层画图软件 alignoth 的路径，**默认路径一般不用改**，只有装到了别处才需要指定。

### 输出配置 | Output

**通俗理解|In plain words:** `-f` 选输出格式：`html`（默认，浏览器可交互缩放/悬停）、`json`（Vega-Lite 原始规范）、`svg`（矢量图，可编辑）、`pdf`（论文用）。`-o` 是输出目录。**日常查看用 html 即可**。

### 可视化参数 | Visualization

**通俗理解|In plain words:** `-d` 控制最多显示多少条 read——深度太高的位点全画会很慢很密，默认 500 已够看；`-w` 是图的最大宽度；`--mismatch-display-min-percent` 是"错配低于这个比例就不标色"的阈值，调小会标出更多错配，调大只看明显的错配。**一般不用动**。

### 高亮选项 | Highlight

**通俗理解|In plain words:** 想把某些位置"圈出来"就用这些参数：`-v` 高亮 VCF 里的变异位点，`--bed` 高亮 BED 区域，`-H name:start-end` 高亮自定义区间（可多次使用）。**不需要高亮就全部省略**。

### 其他选项 | Other

**通俗理解|In plain words:** `-x/--aux-tag` 显示 BAM 里的辅助标签（如 NM 错配数、MD 标签）；`--no-embed-js` 让 html 不内嵌 JS（体积更小但失去交互）；`--plot-all` 强制画出所有 read（无视深度限制）。**一般不用动**。

## 分析流程 | Pipeline

```text
BAM(+索引) + 参考FASTA + 区间
    │
    ▼
校验 BAM/.bai、参考、VCF 索引、alignoth
    │
    ▼
构建 alignoth 命令（按输出格式走不同管道）
    ├─ html → alignoth stdout 直接存为 {bam}.html
    ├─ json → alignoth stdout 存为 {bam}_vl.json
    ├─ svg  → alignoth | vl2vg > {bam}.svg
    └─ pdf  → alignoth | vl2vg | vg2pdf > {bam}.pdf
    │
    ▼
输出单文件 + bam_view.log
```

## 输出 | Output

```text
bam_view_output/
├── {bam前缀}.html            # html 格式（默认，交互式）
├── {bam前缀}_vl.json         # json 格式（Vega-Lite）
├── {bam前缀}.svg             # svg 格式（矢量）
├── {bam前缀}.pdf             # pdf 格式（论文）
└── bam_view.log              # 运行日志
```

> 实际只生成所选 `-f` 格式对应的一个文件；文件名取 BAM 文件名去掉 .bam/.sam/.cram 后缀。

## 结果解读 | Interpreting Results

### 1. html 交互式视图

**通俗理解|In plain words:** 打开 html，横轴是基因组坐标，每一行是一条 read，不同颜色标记不同碱基。悬停可看 read 名、比对详情。

- 某列大量颜色与参考不一致：该位点可能有变异或测序/比对错误
- 大片连续插入/缺失：read 上出现长空白或长色块
- 某位置叠了很多 read（行数骤增）：该区域覆盖深

### 2. 错配标记

**通俗理解|In plain words:** 错配用颜色高亮，直观显示哪些 read 与参考不同。若几乎每条 read 都在同一位置错配，往往是真变异；若零散分布，多为测序噪声。

### 3. 高亮区

**通俗理解|In plain words:** VCF/BED/`-H` 指定的位置会被额外标注，方便把变异和 read 支撑情况对照起来看。

## 参数选择建议 | Parameter Guidance

- 日常查看：默认 `html` 即可
- 要可编辑矢量图：`-f svg`（需安装 vega-cli 的 vl2vg）
- 论文插图：`-f pdf`（需 vega-cli 的 vl2vg + vg2pdf）
- 深度很高的位点：先保持默认 `-d 500`，不够再看 `--plot-all`
- 只想看明显错配：调大 `--mismatch-display-min-percent`（如 5）减少噪点

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--bam, -b` | 必填 |  | BAM文件路径｜BAM file path |
| `--reference, -r` | 必填 |  | 参考序列FASTA文件路径｜Reference FASTA file path |
| `--region, -g` | 必填 |  | 可视化区域(格式: chr:start-end)｜Visualization region (format: chr:start-end) |
| `--alignoth-path` | `~/miniforge3/envs/alignoth/bin/alignoth` |  | alignoth软件路径｜alignoth software path |
| `--output-dir, -o` | `./bam_view_output` |  | 输出目录｜Output directory |
| `--output-format, -f` | `html` | html/json/svg/pdf | 输出格式｜Output format |
| `--max-read-depth, -d` | `500` | int | 最大reads显示深度｜Maximum read depth to display |
| `--max-width, -w` | `1024` | int | 最大宽度｜Maximum width |
| `--mismatch-display-min-percent` | `1.0` | float | 显示错配的最小百分比｜Minimum percentage of mismatches to display |
| `--vcf, -v` | — |  | VCF文件路径(高亮变异位点)｜VCF file path (highlight variants) |
| `--bed` | — |  | BED文件路径(高亮区域)｜BED file path (highlight regions) |
| `--highlight, -H` | — |  | 高亮区间(可多次使用, 格式: name:start-end)｜Highlight interval (can be used multiple times, format: name:start-end) |
| `--aux-tag, -x` | — |  | 辅助标签(可多次使用)｜Auxiliary tag (can be used multiple times) |
| `--no-embed-js` | — |  | 不嵌入JavaScript(仅HTML格式)｜Do not embed JavaScript (HTML format only) |
| `--plot-all` | — |  | 绘制所有reads｜Plot all reads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --bam` | 必填 |  | BAM文件路径｜BAM file path |
| `-r, --reference` | 必填 |  | 参考序列FASTA文件路径｜Reference FASTA file path |
| `-g, --region` | 必填 |  | 可视化区域(格式: chr:start-end)｜Visualization region (format: chr:start-end) |
| `--alignoth-path` | `~/miniforge3/envs/alignoth/bin/alignoth` |  | alignoth软件路径｜alignoth software path |
| `-o, --output-dir` | `./bam_view_output` |  | 输出目录｜Output directory |
| `-f, --output-format` | `html` | html/json/svg/pdf | 输出格式｜Output format |
| `-d, --max-read-depth` | `500` | int | 最大reads显示深度｜Maximum read depth to display |
| `-w, --max-width` | `1024` | int | 最大宽度｜Maximum width |
| `--mismatch-display-min-percent` | `1.0` | float | 显示错配的最小百分比｜Minimum percentage of mismatches to display |
| `-v, --vcf` | — |  | VCF文件路径(高亮变异位点)｜VCF file path (highlight variants) |
| `--bed` | — |  | BED文件路径(高亮区域)｜BED file path (highlight regions) |
| `-H, --highlight` | — | append | 高亮区间(可多次使用, 格式: name:start-end)｜Highlight interval (can be used multiple times, format: name:start-end) |
| `-x, --aux-tag` | — | append | 辅助标签(可多次使用)｜Auxiliary tag (can be used multiple times) |
| `--no-embed-js` | — | store_true | 不嵌入JavaScript(仅HTML格式)｜Do not embed JavaScript (HTML format only) |
| `--plot-all` | — | store_true | 绘制所有reads｜Plot all reads |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- alignoth（conda 环境 `alignoth`，默认 `~/miniforge3/envs/alignoth/bin/alignoth`）
- 可选：vl2vg / vg2pdf（vega-cli，`npm install -g vega-cli vega-lite-cli`，仅 svg/pdf 需要）

## 常见问题 | FAQ

**Q1：报"BAM 索引文件不存在"？**
需要同名 `.bai` 索引，用 `samtools index alignments.bam` 生成。

**Q2：VCF 报"索引文件不存在"？**
VCF 需要 `.csi` 或 `.tbi` 索引，用 `bcftools index` 或 `tabix -p vcf` 生成。

**Q3：svg/pdf 生成失败、提示 vl2vg 未找到？**
svg/pdf 需要 vega-cli 的 vl2vg/vg2pdf，执行 `npm install -g vega-cli vega-lite-cli` 安装。html/json 不需要。

**Q4：CRAM 文件能用吗？**
能，但需要 `-r` 参考 FASTA（CRAM 依赖参考才能解码）；BAM 一般不需要。

**Q5：区间格式怎么写？**
`chr:start-end`，冒号分隔、连字符区间，如 `chr1:1000-2000`。染色体名要和 BAM 里的一致。
