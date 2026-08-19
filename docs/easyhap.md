# easyhap 区域单倍型分析 | EasyHap Regional Haplotype Analysis

一句话人话：给一批有相（phased）基因型的样本，圈出一段基因组区域，看这段区域里样本们各自"是什么版本"（单倍型），哪些样本版本相同、哪些不同，还能按群体比较、画图——不建系统发育树也能快速看出群体内/群体间的单倍型构成。

## 功能概述 | Overview

- 封装上游 EasyHap v1.0（纯 Python，GPL-3.0），支持单倍体/二倍体/多倍体 phased VCF
- 两种分析模式：`inbred`（纯合为主的近交群体，基因型级单倍型）与 `hybrid`（杂合/异交群体，相位拷贝级单倍型）
- 可选 Fisher 精确检验过滤两群体间分化位点、按 Hamming 距离聚类合并冗余单倍型
- 输出 FASTA/NEXUS/PHYLIP 比对（可接下游建树/单倍型网络）与多种图（基因结构、热图、饼图、堆叠柱、性状箱线图）
- 逐区域运行、断点续传、单区域失败不影响其余区域

## 快速开始 | Quick Start

```bash
biopytools easyhap -i sample.phased.vcf.gz --group groups.tsv --region Chr1:1-10000 --plot --gff genes.gff3 -o easyhap_out/
```

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|---|---|
| 单倍型<br>haplotype | 一段染色体上连着一起遗传的一组变异组合，像一本书的"版本" |
| 相位<br>phasing | 知道杂合位点的等位基因分别来自父方还是母方那条染色体（分相以竖线分隔，未分相以斜杠分隔） |
| inbred 模式 | 样本多为纯合（自交/近交群体），把每个样本当作一种"基因型组合" |
| hybrid 模式 | 样本杂合（杂交/异交群体），把每条染色体拷贝单独当作一个单倍型 |
| 聚类<br>cluster | 把相似的单倍型（差异小于阈值）合并成一类，减少稀有类型噪音 |

## 输入 | Input

五类输入中，VCF、分组表、区域三者为必填，性状表和 GFF3/GTF 可选：

1. **VCF（必填，`-i/--input`）**——phased VCF/VCF.gz/BCF：基因型必须已分相（`0|1` 分隔，而非未分相的 `0/1`）；建议使用序列解析型 REF/ALT（PAV 用真实序列表示，不用 `<PAV>` 符号）；大数据集建议 bgzip 压缩并 tabix 建索引。
2. **分组表（必填，`--group`）**——TAB 分隔两列、无表头，格式为 `样本名<TAB>分组名`，样本名必须与 VCF 中一致：
   ```
   Sample1	GroupA
   Sample2	GroupB
   ```
3. **区域（必填其一）**——单区域直接传 `--region Chr1:1-1000`；批量区域用 `--region-file` 传 TAB 三列无表头文件（`chr<TAB>start<TAB>end`），支持 `.gz` 压缩、`#` 开头为注释行。
4. **性状表（可选，`--traits`）**——有表头，首列列名为 Accession，样本名与 VCF 一致；其余列即性状，可用 `--trait-cols` 指定参与分析的列：
   ```
   Accession	Height	Weight
   Sample1	120	30
   ```
5. **GFF3/GTF 注释（可选，`--gff`）**——染色体名必须与 VCF 一致，供基因结构图使用。

## 参数说明 | Parameters

参数按功能分为输入、分析、过滤、出图、通用五组（参数速查表见文末自动生成区块）。

**通俗理解|In plain words — 输入组:** `-i/--input`、`--group`、`--region/--region-file`、`-o/--output-dir` 决定了"分析什么、怎么分组、看哪段、结果放哪"，格式要求见「输入」节，照填即可。

**通俗理解|In plain words — 分析组:** `--mode` 根据群体类型选单倍型重建策略：纯合为主的近交群体（自交系、近交系）选 `inbred`（默认），把每个样本当作一种"基因型组合"；杂交/异交群体样本杂合多，选 `hybrid`，把每条染色体拷贝单独当一个单倍型。`--cluster-threshold` 是"相似度门槛"：调大=更容易把相似单倍型合并成同一类（类更少），调小=更严格（类更多），类比给书按内容相似度归类的松紧度，一般不用动，只在结果单倍型过多/过少时调整。`--hetero-policy` 控制 inbred 模式下杂合位点的编码方式（slash/iupac/missing）；`--vcf-backend` 选 VCF 读取后端，保持 `auto` 即可（没装 cyvcf2/pysam 自动回退纯 Python）。

**通俗理解|In plain words — 过滤组:** `--fisher-groups` 只在关注两个群体间的分化位点时用——给出两个分组名（逗号分隔），仅对这两个分组做 Fisher 精确检验，其余样本不参与该过滤；`--fisher-alpha` 是 P 值阈值，越小越严格、保留位点越少；`--fisher-adjust bh` 是位点很多时的多重检验校正（默认 none 不做校正）。

**通俗理解|In plain words — 出图组:** 只有开了 `--plot` 才生成图形；`--gff` 提供基因注释才有基因结构图，`--traits` 提供性状表才有性状箱线图；`--plot-min-count` 只是"显示门槛"——图中类别计数小于它的不单独显示，**不是过滤**，不会删除任何单倍型数据。

**通俗理解|In plain words — 通用组:** 默认断点续传——已完成的区域（`{label}.easyhap.done` 存在）自动跳过，`--force` 可强制重跑；`--log-level`/`--log-file` 控制日志级别与日志位置（默认 `99_logs/easyhap.log`）。

## 分析流程 | Pipeline

```
读区域清单 → 逐区域独立运行 → 记录版本信息
```

1. **读区域清单**：`--region` 单区域或 `--region-file` 批量区域展开为区域列表；
2. **逐区域独立运行**：已完成区域自动跳过（done-marker），`--force` 强制重跑；
3. **每区域调用上游 EasyHap analyze**：读取区域内变异 →（开启时）Fisher 精确检验过滤两群体分化位点 → 按 `--mode` 重建单倍型 → 按 Hamming 距离聚类合并冗余单倍型 → 导出 FASTA/NEXUS/PHYLIP 比对 →（`--plot` 时）生成各图；
4. **版本记录**：写 `00_pipeline_info/software_versions.yml`。

逐区域独立运行，单区域失败不影响其余区域——某区域无变异或上游报错时仅标记该区域失败，其余区域照常完成，退出码 1 但成功产出保留。

## 输出 | Output

```
easyhap_out/
├── 00_pipeline_info/software_versions.yml   # 工具版本+参数记录
├── 01_haplotypes/                           # 每区域一组文件(前缀 {chr}_{start}_{end})
│   ├── Chr1_1_10000.HapSummary.tsv          # 单倍型总表(每行一个单倍型+聚类ID+样本+计数)
│   ├── Chr1_1_10000.HapGroup.tsv            # 样本→单倍型/聚类/分组/性状
│   ├── Chr1_1_10000.Haplotype.fa/.nex/.phy  # 单倍型比对(供下游建树)
│   ├── Chr1_1_10000.Haplotype_sample.fa     # 样本拷贝级序列
│   ├── Chr1_1_10000.AlleleStateMap.tsv      # 原始等位基因→编码状态映射
│   ├── Chr1_1_10000.ProcessedVariants.tsv / .SampleGenotypeTokens.tsv  # 处理表(--no-processed 可关)
│   ├── Chr1_1_10000.FisherFilter.tsv        # Fisher 过滤统计(开启时)
│   └── Chr1_1_10000.{GeneHaplotype,HaplotypeHeatmap,GroupPie,GroupStackedBar,TraitBoxplot}.pdf  # --plot 时
└── 99_logs/easyhap.log                      # 运行日志
```

关键文件说明：

- `HapSummary.tsv`：单倍型总表（行=单倍型，列=位点，Accession=样本清单，Number=计数），回答"有哪些版本、每个版本多少人"
- `HapGroup.tsv`：样本维度表（样本→单倍型/聚类/分组/性状），供下游统计
- `Haplotype.fa/.nex/.phy`：非冗余单倍型比对，可直接喂 IQ-TREE 等下游建树/单倍型网络
- `Haplotype_sample.fa`：样本（拷贝）级序列
- `AlleleStateMap.tsv`：原始 token→编码状态映射，看懂 HapSummary 编码用
- `FisherFilter.tsv`：各位点两群体 2x2 表 + pvalue/padj/keep（开启 Fisher 时）
- 图文件按 `--plot` 参数组合出现：基因结构（GeneHaplotype）、热图（HaplotypeHeatmap）、饼图（GroupPie）、堆叠柱（GroupStackedBar）、性状箱线图（TraitBoxplot）

## 结果解读 | Interpreting Results

- **HapSummary 行数 = 单倍型数**：单倍型过多（接近样本数）提示区域变异过杂或聚类阈值太严，可调大 `--cluster-threshold` 或缩小区域
- **HapGroup 中 ClusterID 相同的样本 = 相似单倍型**；群体饼图中某单倍型集中在单一群体 = 潜在分化信号（结合 FisherFilter 表确认）
- **性状箱线图**：组间显著性（p<0.05，`***`）提示该单倍型与性状关联
- **inbred 模式**下含 `/` 的 token（如 `A/G`）= 杂合位点，按 `--hetero-policy` 策略保留

## 参数选择建议 | Parameter Guidance

| 场景<br>Scenario | 推荐参数<br>Recommended flags |
|---|---|
| 纯合自交群体 | `--mode inbred --hetero-policy slash` |
| 杂合/异交群体 | `--mode hybrid` |
| 关注两群体分化 | `--fisher-groups A,B --fisher-alpha 0.05 --fisher-adjust bh` |
| 发文章出图 | `--plot --plot-format pdf --plot-hap-level cluster --cluster-threshold 0.5` |
| 大 VCF | 先 bgzip+tabix 建索引，再 `--vcf-backend auto` |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | phased VCF/VCF.gz/BCF 文件｜Phased VCF/VCF.gz/BCF file |
| `--group` | 必填 |  | 样本分组表(样本<TAB>分组, 无表头)｜Sample group table |
| `--region` | — |  | 单区域 chr:start-end｜Single region chr:start-end |
| `--region-file` | — |  | 批量区域文件(TAB 三列)｜Batch region file |
| `-o, --output-dir` | `./easyhap_output` |  | 输出目录(默认./easyhap_output)｜Output directory |
| `--mode` | `inbred` | inbred/hybrid | 单倍型重建策略(默认inbred)｜Mode (default inbred) |
| `--hetero-policy` | `slash` | slash/iupac/missing | inbred 模式杂合编码(默认slash)｜Hetero policy (default slash) |
| `--cluster-threshold` | `0.15` | float | 聚类阈值(默认0.15)｜Cluster threshold (default 0.15) |
| `--vcf-backend` | `auto` | auto/cyvcf2/pysam/plain | VCF 读取后端(默认auto)｜VCF backend (default auto) |
| `--no-processed` | `False` |  | 不写处理表｜Skip processed tables |
| `--fisher-groups` | — |  | Fisher 过滤两分组名｜Fisher filter groups |
| `--fisher-alpha` | — | float | Fisher 过滤 P 值阈值｜Fisher filter p cutoff |
| `--fisher-adjust` | `none` | none/bh | 多重检验校正(默认none)｜Fisher adjust (default none) |
| `--plot` | `False` |  | 生成图形输出｜Generate plots |
| `--gff` | — |  | GFF3/GTF 注释｜GFF3/GTF annotation |
| `--traits` | — |  | 性状表｜Trait table |
| `--trait-cols` | — |  | 性状列名(逗号分隔)｜Trait columns |
| `--plot-format` | `pdf` |  | 图格式(默认pdf)｜Plot formats (default pdf) |
| `--plot-hap-level` | `hap` | hap/cluster | 按单倍型/聚类出图(默认hap)｜Plot level (default hap) |
| `--plot-min-count` | `1` | int | 图中最小类别计数(默认1)｜Min count in plots (default 1) |
| `--force` | `False` |  | 强制重跑已完成区域｜Force re-run |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |
| `--log-file` | — |  | 日志文件(默认99_logs/easyhap.log)｜Log file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | phased VCF/VCF.gz/BCF 文件｜Phased VCF/VCF.gz/BCF file |
| `--group` | 必填 |  | 样本分组表(样本<TAB>分组, 无表头)｜Sample group table (sample<TAB>group, no header) |
| `--region` | — |  | 单区域 chr:start-end｜Single region chr:start-end |
| `--region-file` | — |  | 批量区域文件(TAB 三列 chr start end)｜Batch region file |
| `-o, --output-dir` | `./easyhap_output` |  | 输出目录(默认./easyhap_output)｜Output directory (default ./easyhap_output) |
| `--mode` | `inbred` | inbred/hybrid | 单倍型重建策略(默认inbred)｜Haplotype reconstruction mode (default inbred) |
| `--hetero-policy` | `slash` | slash/iupac/missing | inbred 模式杂合编码(默认slash)｜Hetero encoding in inbred mode (default slash) |
| `--cluster-threshold` | `0.15` | float | 单倍型聚类阈值(默认0.15)｜Clustering threshold (default 0.15) |
| `--vcf-backend` | `auto` | auto/cyvcf2/pysam/plain | VCF 读取后端(默认auto)｜VCF reader backend (default auto) |
| `--no-processed` | — | store_true | 不写 ProcessedVariants/SampleGenotypeTokens 表｜Skip processed tables |
| `--fisher-groups` | — |  | Fisher 过滤两分组名(逗号分隔)｜Two group names for Fisher filter |
| `--fisher-alpha` | — | float | Fisher 过滤 P 值阈值(0,1]｜Fisher filter p cutoff (0,1] |
| `--fisher-adjust` | `none` | none/bh | 多重检验校正(默认none)｜Multiple-testing correction (default none) |
| `--plot` | — | store_true | 生成图形输出｜Generate plots |
| `--gff` | — |  | GFF3/GTF 注释(基因结构图)｜GFF3/GTF annotation for gene plot |
| `--traits` | — |  | 性状表(有表头)｜Trait table (with header) |
| `--trait-cols` | — |  | 性状列名(逗号分隔)｜Comma-separated trait columns |
| `--plot-format` | `pdf` |  | 图格式(逗号分隔 pdf/svg/png, 默认pdf)｜Plot formats (default pdf) |
| `--plot-hap-level` | `hap` | hap/cluster | 按单倍型/聚类出图(默认hap)｜Plot by hap or cluster (default hap) |
| `--plot-min-count` | `1` | int | 图中最小类别计数(默认1)｜Min class count in plots (default 1) |
| `--force` | — | store_true | 强制重跑已完成区域｜Force re-run completed regions |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |
| `--log-file` | — |  | 日志文件(默认99_logs/easyhap.log)｜Log file (default 99_logs/easyhap.log) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- 上游 EasyHap v1.0（pip 安装于 pop 域环境）；pandas/numpy/matplotlib/scipy（pop 环境自带）
- cyvcf2/pysam 为上游可选加速后端（pop 环境未装，`--vcf-backend auto` 自动回退 plain）

## 常见问题 | FAQ

- **Q: 报错 "No variants found in region"?** 该区域在 VCF 中无变异位点，上游直接报错；本模块会标记该区域失败并继续其余区域（退出码 1 但成功产出保留）。
- **Q: 重跑会覆盖已有结果吗?** 不会——已完成的区域（`{label}.easyhap.done` 存在）自动跳过；要重跑用 `--force`。
- **Q: 为什么默认 inbred 模式?** 上游默认即 inbred；杂交/异交群体请用 `--mode hybrid`。
- **Q: pop 环境没有 cyvcf2 会怎样?** `--vcf-backend auto` 自动回退纯 Python plain 后端，功能一致只是大文件略慢。
