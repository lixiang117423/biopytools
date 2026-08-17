# SubPhaser 异源多倍体亚基因组分离 | SubPhaser Subgenome Phasing

一句话理解：**基于「重复序列 k-mer」的差异，自动把异源多倍体基因组的染色体分成几个亚基因组并命名（如 Chr1A / Chr1B）**，解决「多倍体组装完成后，不知道哪些染色体属于哪个祖先」的问题。

## 功能概述 | Overview

- 封装 SubPhaser：基于重复 k-mer 的异源多倍体亚基因组自动分离与命名
- 自动模式（默认）：无需配置文件，自动按染色体大小配对、生成临时配置并运行
- 配置模式：提供 `-c` 亚基因组配置文件，精细控制分组
- 验证模式：提供两个父本基因组，用比对结果把亚基因组标签校正为 A/B
- 后处理自动生成 Chr1A/B/C 风格命名、重命名 FASTA、phasing 结果汇总与最终配置文件
- 支持禁用 LTR 分析、Circos 图、同源区块等耗时步骤

## 快速开始 | Quick Start

```bash
biopytools subphaser -i genome.fa --nsg 2 -o output/
```

最小输入：一个基因组 FASTA 和亚基因组数量（--nsg）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 异源多倍体 | 由两个不同物种杂交加倍形成的多倍体（如小麦、油菜） |
| 亚基因组 | 异源多倍体里来自某个祖先的那套染色体（A 亚基因组、B 亚基因组） |
| k-mer | 从序列里滑窗取出的长度为 k 的小片段 |
| 重复 k-mer | 在基因组里出现多次的 k-mer；不同亚基因组的重复序列组成有差异，可用来区分 |
| phasing（分相） | 把混在一起的多套染色体分开、各自还原 |
| LTR | 一类转座子（长末端重复），可用于估算亚基因组分化时间 |
| bootstrap | 统计里「重抽样」检验，用来判断差异是否显著 |
| Chr1A / Chr1B | 分相后的命名：第 1 号同源组里来自 A、B 两个亚基因组的两条染色体 |

## 输入 | Input

### 基因组

`-i` / `--genomes`，FASTA 格式，可指定多个。

### 亚基因组数量

`--nsg`，必须 >= 2。

### 可选输入

- `-c/--sg-cfgs`：亚基因组配置文件（不提供则进入自动模式）。
- `--parental-genomes`：两个父本基因组，启用验证模式（校正 A/B 标签）。

```text
# 自动模式（最简）
biopytools subphaser -i genome.fa --nsg 2 -o output/

# 验证模式（用父本校正A/B标签）
biopytools subphaser -i genome.fa --nsg 2 --parental-genomes parentA.fa parentB.fa
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 基因组（-i）和亚基因组数量（--nsg）必填。nsg 按物种实际亚基因组数填（如异源四倍体填 2，异源六倍体小麦填 3），填错会整体分组错误。

### 输入模式 | Input mode

**通俗理解|In plain words:** 不给 `-c` 就是自动模式（程序自己按染色体大小配对分列）；给了 `-c` 就是配置模式（用你准备的配置文件精细控制）。`--parental-genomes` 加两个父本进入验证模式，自动把分好的组校正为 A/B 标签。

### 输出与资源 | Output and resources

**通俗理解|In plain words:** `-o` 是输出目录，`--prefix` 是输出前缀（默认无）。`-t` 线程数按机器资源给（默认 24）。

### 染色体过滤 | Chromosome filtering

**通俗理解|In plain words:** `--min-chrom-size` 过滤太短的序列（默认 100 万 bp），短于它的 contig 不参与分相。若过滤后染色体数少于 nsg 会报错，此时调小该值。大基因组的小 contig 多时保持默认即可。

### K-mer 参数 | K-mer

**通俗理解|In plain words:** 这是 SubPhaser 的核心参数。`-k` 是 k-mer 长度（默认 15）；`-f/--min-fold` 是最小倍数差异（判断某 k-mer 是否在亚基因组间富集）；`-q/--min-freq` 是最小 k-mer 频率。这些是上游工具的默认值，一般不用动。

### 聚类与统计 | Clustering and statistics

**通俗理解|In plain words:** `--max-pval` 是统计检验的显著性阈值（默认 0.05），`--replicates` 是 bootstrap 重复次数（默认 1000），`--test-method` 选统计检验方法（默认 ttest_ind）。重复次数越多越稳但越慢，一般用默认。

### 步骤控制 | Step control

**通俗理解|In plain words:** LTR 分析（估算分化时间）、Circos 图、同源区块这几个步骤耗时较长，大基因组可用 `--disable-ltr` / `--disable-circos` / `--disable-blocks` 关掉提速；`--just-core` 只跑核心分相。只想快速出分组结果就用这些开关。

### LTR 与 Circos | LTR and Circos

**通俗理解|In plain words:** 仅当开启 LTR/Circos 分析时才相关。`--ltr-detectors` 选 LTR 检测工具，`--mu` 是替换率/年（用于分化时间估算）；`--window-size` 是 Circos 窗口，`--aligner` 是画同源区块的比对工具。默认即可。

### 高级参数 | Advanced

**通俗理解|In plain words:** `--sg-assigned` 提供已知分配文件可跳过聚类；`--target` 指定目标染色体；`--labels` 自定义标签；`--no-label` 不加标签前缀；`--custom-features` 附加自定义特征。这些都是进阶用法，正常流程用不到。

## 分析流程 | Pipeline

```text
输入基因组FASTA + nsg
  -> [自动模式] 解析染色体 -> 按大小配对 -> 生成临时配置(00_auto_config/)
  -> 运行SubPhaser(重复k-mer聚类分相)
  -> 解析chrom-subgenome.tsv -> 分亚基因组
  -> [验证模式] 用父本比对校正A/B标签
  -> 生成Chr1A/B/C命名 -> 重命名FASTA(phased_genome/)
  -> 写phasing_result.tsv + 最终配置文件phased_sg.config
```

## 输出 | Output

```text
subphaser_output/
├── 00_auto_config/
│   └── auto_sg.config               # 自动模式生成的临时配置
├── 00_pipeline_info/
│   └── software_versions.yml        # 软件版本信息
├── phased_genome/
│   └── phased_genome.fa             # 重命名后的基因组(Chr1A/Chr1B...)(核心)
├── phasing_result.tsv               # phasing结果汇总表
├── phased_sg.config                 # 最终配置文件
├── ...                              # SubPhaser自身产出的其他文件(图/报告)
└── 99_logs/
    └── subphaser.log
```

关键文件说明：

- **phased_genome.fa**：染色体已按 Chr1A/B/C 风格重命名的基因组，可直接用于下游。
- **phasing_result.tsv**：五列（原始名、染色体编号、亚基因组、新命名、长度），是分相结果的完整清单。
- **phased_sg.config**：按同源组组织的最终配置，可复用/手工微调。
- **软件版本**：记录 SubPhaser 版本、conda 环境与关键参数。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看 phasing 结果，核心是「每个亚基因组是否分到了合理数量的染色体、同源组是否配对正确」。

- **phasing_result.tsv**：每条染色体属于哪个亚基因组、被命名成 Chr 几 A/B/C。同一「ChrN」的不同亚基因组（如 Chr1A、Chr1B）应是同源关系。
- **染色体数量**：每个亚基因组的染色体数应接近 nsg 均分后的预期；某组明显偏多/偏少说明分相或 nsg 设置有误。
- **验证模式**：用父本校正后，A/B 标签对应该亚基因组与哪个父本更相似，比默认的 SG1/SG2 更有生物学意义。
- **日志**：99_logs/subphaser.log 记录了每条染色体到新命名的映射，可逐一核对。

## 参数选择建议 | Parameter Guidance

- **--nsg**：按物种实际亚基因组数填，务必准确（异源四倍体=2，异源六倍体=3）。
- **--min-chrom-size**：只想要染色体级别的分相结果时保持默认；若小 contig 也要参与，可调小，但别小于过滤后染色体数不足 nsg。
- **--disable-ltr / --disable-circos / --disable-blocks**：大基因组或只想快速出分相结果时关闭，显著提速。
- **--parental-genomes**：有父本参考时强烈建议加，让 A/B 标签有生物学依据。
- **--just-core**：只想验证分相是否成功时，先跑核心分相看结果，确认后再补全分析。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genomes` | 必填 |  | 基因组FASTA文件(可多个)｜Genome FASTA files |
| `--nsg` | 必填 | int | 亚基因组数量(>=2)｜Number of subgenomes (>=2) |
| `-c, --sg-cfgs` | — |  | 亚基因组配置文件(可选，不提供则自动模式)｜Subgenome config files (optional, auto mode if omitted) |
| `--parental-genomes` | — |  | 父本基因组(验证模式，需2个)｜Parental genomes for validation (requires 2) |
| `-o, --output-dir` | `./subphaser_output` |  | 输出目录｜Output directory |
| `--prefix` | — |  | 输出前缀｜Output prefix |
| `-t, --threads` | `24` | int | 线程数｜Number of threads |
| `--min-chrom-size` | `1000000` | int | 最小染色体长度(bp)，过滤小contigs｜Min chromosome size (bp), filter small contigs |
| `-k, --kmer-size` | `15` | int | K-mer大小｜K-mer size |
| `-f, --min-fold` | `2.0` | float | 最小倍数差异｜Minimum fold difference |
| `-q, --min-freq` | `200` | int | 最小k-mer频率｜Minimum k-mer frequency |
| `--max-pval` | `0.05` | float | 最大P值｜Maximum P-value |
| `--replicates` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates |
| `--test-method` | `ttest_ind` | ttest_ind/kruskal/wilcoxon/mannwhitneyu | 统计检验方法｜Statistical test method |
| `--disable-ltr` | — |  | 禁用LTR分析(大基因组耗时较长)｜Disable LTR analysis |
| `--disable-circos` | — |  | 禁用Circos图｜Disable Circos plot |
| `--disable-blocks` | — |  | 禁用同源区块｜Disable homologous blocks |
| `--just-core` | — |  | 仅运行核心phasing｜Only run core phasing |
| `--ltr-detectors` | — | ltr_finder/ltr_harvest | LTR检测工具｜LTR detection tools |
| `--mu` | `1.3e-08` | float | 替换率/年｜Substitution rate per year |
| `--window-size` | `1000000` | int | Circos窗口大小｜Circos window size (bp) |
| `--aligner` | `minimap2` | minimap2/unimap | 比对工具｜Aligner for homologous blocks |
| `--sg-assigned` | — |  | 已知亚基因组分配文件(跳过聚类)｜Pre-assigned subgenome file (skip clustering) |
| `--target` | — |  | 目标染色体文件｜Target chromosomes file |
| `--labels` | — |  | 基因组标签｜Genome labels |
| `--no-label` | — |  | 不添加标签前缀｜No label prefix |
| `--custom-features` | — |  | 自定义特征FASTA｜Custom feature FASTA files |
| `--figfmt` | `pdf` | pdf/png | 图片格式｜Figure format |
| `--overwrite` | — |  | 覆盖已有结果｜Overwrite existing results |
| `--cleanup` | — |  | 清理临时文件｜Clean up temporary files |
| `--conda-env` | `SubPhaser` |  | conda环境名称｜Conda environment name |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- SubPhaser（conda 环境，默认 `SubPhaser`，通过 `--conda-env` 可改；以 `conda run -n SubPhaser --no-capture-output subphaser` 调用）
- 验证模式额外需要 minimap2（在 PATH 中）

## 常见问题 | FAQ

**Q1：报「过滤后仅剩 N 条染色体 (< nsg)」？**
`--min-chrom-size` 把太多短序列过滤掉了。调小该值（或设 0 关闭过滤）再跑。

**Q2：subphaser 命令未找到？**
确认 conda 环境 `SubPhaser`（或用 `--conda-env` 指定的环境）已安装 SubPhaser。

**Q3：断点续传支持吗？**
本封装层不额外做断点续传，直接委托 SubPhaser 自身；用 `--overwrite` 可覆盖已有结果，`--cleanup` 清理临时文件。

**Q4：自动模式和配置模式选哪个？**
首次跑用自动模式（最省事）；自动配对不理想（如同源染色体大小差异大导致配错）时，手写 `-c` 配置文件精细控制。

**Q5：A/B 标签是随便定的吗？**
默认是 SG1/SG2 等，无生物学含义；加 `--parental-genomes` 用父本比对后，才会校正成对应父本的 A/B 标签。

