# GenomeSyn2 比较基因组学可视化 | GenomeSyn2 Comparative Genomics Visualization

一句话理解：**把多个基因组、基因或 SNP 数据「画到同一张图里」做横向比较**——比共线性、比序列差异、比 SNP 密度——本质上是把成熟的 GenomeSyn2（Perl）工具做了统一入口，按你要做的事自动选对命令。

## 功能概述 | Overview

- 封装 GenomeSyn2（Perl 脚本，胡教授实验室出品）的比较基因组可视化能力，统一命令行入口
- 五种工作模式按参数自动切换：序列比对、VCF 算 SNP、血统解析图、共线性图、文件列表生成
- 支持 DNA 与蛋白质两类比对（mummer/minimap2 对基因组，blastp/mmseqs/diamond 对蛋白）
- 依赖 conda 环境 genomesyn2 提供的 perl 与 GenomeSyn2.pl，无需手动装 Perl 模块
- 无断点续传（每条命令独立执行，见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools genomesyn2 --align mummer --genome ./genome_dir/ --outdir ./output/ -t 12
```

最小输入：一个按数字排序命名的基因组文件目录，用 mummer 做两两比对并出图。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 比较基因组学 | 把几个基因组摆在一起找「同」与「不同」的学问，本工具负责把结果画出来 |
| 共线性图 | 一条条染色体排成行，用连线标出两两之间的同源段落，一眼看出重排/倒位 |
| SNP | 单个碱基位点上的差异（单核苷酸多态性），是群体间差异的最小单位 |
| bin | 把染色体切成固定大小的「格子」，每个格子里数有多少 SNP，得到密度分布 |
| 血统解析图 | 用 SNP 一致性 + 密度两类数据，画出基因组里哪些区段更像哪个祖先来源 |
| BLAST/MMseqs/Diamond | 三类做蛋白质比对的软件，速度与灵敏度各有取舍 |
| 配置文件（conf） | GenomeSyn2 绘图的「配方」文件，写清楚画什么、颜色、标注 |

## 输入 | Input

输入取决于你选哪种模式（见下方参数说明），常见组合：

- **序列比对模式**：`--genome` 基因组文件目录 +（蛋白比对时再加 `--gene` 基因注释目录）+ `--outdir`
- **VCF 模式**：`--vcf` 一个 VCF 文件 + `--bin`
- **血统解析模式**：`--identity`（SNP 一致性 BED）+ `--density`（SNP 密度 BED）
- **共线性绘图模式**：`--conf` 配置文件
- **文件列表生成模式**：`--type`（fa/prot/anno）+ `--path` 目录 + `--out` 输出文件名

基因组目录里的文件需按数字排序命名（如 01.fa、02.fa），便于程序确定顺序。

## 参数说明 | Parameters

### 比对模式 | Alignment mode

**通俗理解|In plain words:** 用 --align 选比对软件，DNA 用 mummer/minimap2，蛋白质用 blastp/mmseqs/diamond。选了蛋白质比对就**必须**再给 --gene（基因注释目录）。--genome 里的文件要按数字排序。

- --align：mummer / minimap2 / blastp / mmseqs / diamond
- --genome：基因组文件目录（文件名需按数字排序）
- --gene：基因注释文件目录（蛋白质比对必填）
- --outdir：输出目录

### VCF 模式 | VCF mode

**通俗理解|In plain words:** 给一个 VCF（记录了哪些位点有变异），程序按固定大小的 bin 统计每个格子的 SNP 密度与一致性。--bin 是格子大小，调大 = 图更粗但更平滑，调小 = 更细但噪点多。**一般默认 50000 就行**。

- --vcf：VCF 文件路径
- --bin：bin 大小（bp），默认 50000

### 血统解析绘图模式 | Ancestry plotting mode

**通俗理解|In plain words:** 直接给已经算好的 SNP 一致性文件（identity）和 SNP 密度文件（density），画出血统解析图。这两个文件通常由 VCF 模式先产出。

- --identity：SNP 一致性 BED 文件
- --density：SNP 密度 BED 文件

### 共线性绘图模式 | Synteny plotting mode

**通俗理解|In plain words:** 给一份 GenomeSyn2 的配置文件（conf）就出图。--anno 用于查看/生成注释配置选项。想精细控制图的每个细节就用这个模式。

- --conf：配置文件路径
- --anno：显示注释配置选项

### 文件列表生成模式 | File list generation mode

**通俗理解|In plain words:** 只生成一份「文件清单」（把某个目录下某类文件列成表），不画图。--type 决定列 fa/prot/anno 哪类文件，--path 指目录，--out 是输出文件名。

- --type：fa / prot / anno
- --path：文件路径（目录）
- --out：输出文件名

### 通用参数 | Common

**通俗理解|In plain words:** --threads 是并行线程数，比对类任务核多就给大点。

- -t, --threads：线程数，默认 12

## 分析流程 | Pipeline

```text
根据参数自动判定模式
    |
    v
  --align      -> 序列比对（DNA/蛋白）
  --vcf        -> 计算 SNP 密度与一致性
  --identity + --density -> 血统解析图
  --conf       -> 共线性绘图
  --type       -> 生成文件列表
    |
    v
conda run -n genomesyn2 perl GenomeSyn2.pl <参数>
```

## 输出 | Output

输出由 GenomeSyn2 脚本决定，落在 --outdir（比对模式）或当前目录（其他模式）；文件列表模式输出到 --out 指定的文件。常见产物：

- 序列比对模式：比对结果文件与共线性图（SVG/PNG 等，具体命名由 GenomeSyn2 决定）
- VCF 模式：SNP 一致性（identity）与密度（density）BED 文件，供后续血统解析绘图用
- 血统解析/共线性绘图模式：最终图片
- 文件列表模式：一份文件清单（TSV）
- 每次运行在输出目录写 genomesyn2.log 日志

## 结果解读 | Interpreting Results

- **共线性图**：连线规整平行 = 基因组间保守；大量交叉/断线 = 结构重排。
- **SNP 密度图**：密度高峰区 = 变异富集（可能提示选择压力或比对噪声）；密度低谷 = 保守区。
- **血统解析图**：区段颜色越单一 = 该区段越像某个特定祖先来源；颜色混杂 = 祖先来源混合。
- **好坏判据**：脚本返回码为 0 且输出文件生成即成功；日志出现「命令执行失败」需看 stderr。

## 参数选择建议 | Parameter Guidance

- **DNA 比对**：默认示例用 mummer（精细）；想快用 minimap2。
- **蛋白比对**：blastp 最经典稳妥，mmseqs/diamond 更快（大数据量首选）。
- **--bin**：基因组大、想看得细可降到 10000–20000；想画全局趋势用默认 50000 或更大。
- **--conf + --anno**：需要精细排版时，先用 --anno 看有哪些配置项，再写好 conf 绘图。
- **文件列表**：不确定 GenomeSyn2 要什么输入格式时，先用 --type 生成清单对照。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--align` | — | mummer/minimap2/blastp/mmseqs/diamond | 比对软件类型｜Alignment software type |
| `--genome` | — |  | 基因组文件目录｜Genome files directory |
| `--gene` | — |  | 基因注释文件目录｜Gene annotation files directory |
| `--outdir` | — |  | 输出目录｜Output directory |
| `--vcf` | — |  | VCF文件路径｜VCF file path for SNP analysis |
| `--bin` | `50000` | int | Bin大小(用于SNP分析)｜Bin size for SNP analysis |
| `--identity` | — |  | SNP一致性文件｜SNP identity BED file |
| `--density` | — |  | SNP密度文件｜SNP density BED file |
| `--conf` | — |  | 配置文件路径｜Configuration file path |
| `--anno` | — |  | 显示注释配置选项｜Show annotation configuration options |
| `--type` | — | fa/prot/anno | 文件类型(用于生成文件列表)｜File type for generating file list |
| `--path` | — |  | 文件路径｜File path for generating list |
| `--out` | — |  | 输出文件名｜Output file name |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--align` | — | mummer/minimap2/blastp/mmseqs/diamond | 比对软件类型｜Alignment software type (mummer/minimap2/blastp/mmseqs/diamond) |
| `--genome` | — |  | 基因组文件目录｜Genome files directory (文件名需按数字排序｜Files must be numbered) |
| `--gene` | — |  | 基因注释文件目录｜Gene annotation files directory (蛋白质比对需要｜required for protein alignment) |
| `--outdir` | — |  | 输出目录｜Output directory |
| `--vcf` | — |  | VCF文件路径｜VCF file path for SNP analysis |
| `--bin` | `50000` | int | Bin大小(用于SNP分析)｜Bin size for SNP analysis (default: 50000) |
| `--identity` | — |  | SNP一致性文件｜SNP identity BED file |
| `--density` | — |  | SNP密度文件｜SNP density BED file |
| `--conf` | — |  | 配置文件路径｜Configuration file path |
| `--anno` | — | store_true | 显示注释配置选项｜Show annotation configuration options |
| `--type` | — | fa/prot/anno | 文件类型(用于生成文件列表)｜File type for generating file list (fa/prot/anno) |
| `--path` | — |  | 文件路径｜File path for generating list |
| `--out` | — |  | 输出文件名｜Output file name |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--perl-path` | `~/miniforge3/envs/genomesyn2/bin/perl` |  |  |
| `--genomesyn2-pl` | `~/miniforge3/envs/genomesyn2/bin/GenomeSyn2.pl` |  |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **conda 环境 genomesyn2**（提供 perl 与 GenomeSyn2.pl，默认路径 ~/miniforge3/envs/genomesyn2/）
- **Perl 模块**：Bio::SeqIO（Bioperl）、SVG（缺失会告警但仍尝试运行）
- 按模式另需：mummer / minimap2 / blastp / mmseqs / diamond（比对软件）
- Python 3（封装脚本自身）

## 常见问题 | FAQ

**Q1：换了参数重跑，结果没更新？**
本工具没有断点续传，每次都会重新执行 GenomeSyn2.pl；若结果没变，通常是 GenomeSyn2 脚本自身对已有同名输出不覆盖，或你指定的 --outdir 没变、旧文件被复用。必要时换个输出目录。

**Q2：报「GenomeSyn2 Perl 脚本不存在 / Perl 不存在」？**
默认路径是 ~/miniforge3/envs/genomesyn2/bin/GenomeSyn2.pl 与同目录 perl。若你的环境名或安装位置不同，用隐藏参数 --perl-path / --genomesyn2-pl 指定（见脚本 help 被隐藏的两个参数）。

**Q3：蛋白比对报「需要 --gene 参数」？**
blastp/mmseqs/diamond 是蛋白比对，必须有基因注释目录（--gene）提供蛋白序列。补上 --gene 即可。

**Q4：Perl 模块告警会影响结果吗？**
程序会检查 Bio::SeqIO 与 SVG 模块，缺失只发 WARNING 不阻断；但相关功能（读序列/画 SVG）可能失败。建议在 genomesyn2 环境里装齐 Bioperl 与 SVG。

**Q5：文件列表模式（--type）是干什么的？**
它只是生成文件清单，不跑分析。当你需要给 GenomeSyn2 准备输入列表（哪些 fa/蛋白/注释文件）时用，产出 --out 指定的清单文件。
