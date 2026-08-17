# OrthoFinder 泛基因组分析 | OrthoFinder Pangenome Analysis

一句话理解：把多个物种或个体的蛋白质序列放在一起，先自动找出"同一个祖先基因的后代"(直系同源群)，再把它们分成**核心基因、软核心基因、非必需基因、私有基因**四类，回答"哪些是大家共有的家底、哪些是可有可无的"。

## 功能概述 | Overview

- 以 OrthoFinder 为引擎，先做直系同源群(orthogroup)推断，再做泛基因组四分类
- 自动分类：core(所有基因组都有)、softcore(最多缺 N 个)、dispensable(非必需)、private(仅一个基因组有)
- 额外识别单拷贝核心基因(每个基因组恰好 1 份)，可按同源群或按基因组提取序列
- 内置稀释曲线(rarefaction)分析，估计"再加样本还能发现多少新基因"
- 可选生成分类饼图、频率分布图、稀释曲线图与综合文字报告
- 断点续传：自动复用已完成的 OrthoFinder 结果，跳过重复计算

## 快速开始 | Quick Start

```bash
biopytools orthofinder -i protein_sequences/ -o pangenome_results
```

`-i` 指向一个只放 FASTA 文件的目录（默认按蛋白质序列处理，每个文件是一个基因组/样本）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 直系同源群(orthogroup) | 一群"来自同一个祖先基因"的基因，像一大家人各有各的分支 |
| 泛基因组(pangenome) | 一个物种里所有个体基因的"总和集合"，包括大家共有的和个别人独有的 |
| 核心基因(core) | 每个基因组都有的基因，像"基本器官"，缺了就活不了 |
| 软核心基因(softcore) | 几乎所有基因组都有、只缺极少数的基因，重要但非绝对必需 |
| 非必需基因(dispensable) | 只在部分基因组出现的基因，提供环境适应性等"锦上添花"功能 |
| 私有基因(private) | 只在一个基因组里出现的基因，往往是这个个体/菌株的独门本领 |
| 单拷贝基因(single copy) | 每个基因组里恰好只有一份的核心基因，常用作建树或比较的"标准尺子" |
| MCL inflation | 给同源群"分家"的松紧度；越大分得越细(群越多越小)，越小越粗 |
| 稀释曲线(rarefaction) | 模拟"只测 1 个、2 个、3 个样本能发现多少基因"，看泛基因组是否已经"饱和" |

## 输入 | Input

输入为一个**目录**，目录里每个 FASTA 文件代表一个基因组/样本的序列集合。

- 默认按**蛋白质序列**处理（`.fa`、`.faa`、`.fasta`、`.pep` 等扩展名均会被识别）
- 若输入是 DNA 序列，加 `-d/--dna` 切换
- 文件名(去扩展名)即基因组名，会出现在所有输出表里
- 序列 ID 中若含非法字符，程序会先尝试清洗后再验证

```text
protein_sequences/
├── strain_A.faa
├── strain_B.faa
└── strain_C.faa
```

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 给"材料"(蛋白质序列目录)，`-o` 给"成果存放地"；`-n` 给这次分析起个名字，不填就用目录名自动生成。

### 分类阈值 | Classification thresholds

**通俗理解|In plain words:** 决定"缺几个基因组还算软核心"的边界线。`--softcore-threshold` 是实际生效的判定线：缺失数小于等于它就归为软核心；缺失更多且多个基因组共有，则归为非必需；只有一个基因组有的归私有。`--dispensable-threshold` 目前只影响摘要文字里的定义描述，不参与实际判定(见 FAQ)。**默认 1，一般不用动。**

### OrthoFinder 引擎参数 | OrthoFinder engine

**通俗理解|In plain words:** 这些参数透传给 OrthoFinder 本体。`-s` 选"拿什么搜索引擎找同源"(blast 最经典最慢，diamond/mmseqs 快很多)；`--mcl-inflation` 控制同源群"分家松紧"；`-t` 是线程数。换搜索程序能显著改变速度，但会轻微影响结果。

### 分析开关 | Analysis toggles

**通俗理解|In plain words:** 控制"除了基础分类还要不要做额外分析"。单拷贝提取、稀释曲线、画图默认都开；不需要哪块就用对应 `--disable-*`/`--no-plots` 关掉以提速。`--basic-only`/--generate-trees` 控制 OrthoFinder 是否额外构建物种树(默认只跑同源群推断，不建树)。

### 断点续传 | Resume

**通俗理解|In plain words:** 默认自动复用已完成的 OrthoFinder 结果，重跑不会从头再来。想彻底重算加 `--force`；已有结果但想跳过 OrthoFinder 直接做分类用 `--skip-orthofinder`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入蛋白质序列文件目录｜Input protein sequence files directory |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--project-name, -n` | — | str | 项目名称｜Project name |
| `--softcore-threshold` | `1` | int | Softcore基因缺失阈值｜Softcore gene missing threshold |
| `--dispensable-threshold` | `1` | int | Dispensable基因缺失阈值｜Dispensable gene missing threshold |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--search, --search-program, -s` | `blast` | blast/diamond/diamond_ultra_sens/mmseqs/blast_nucl | 序列搜索程序｜Sequence search program |
| `--mcl-inflation` | `1.2` | float | MCL膨胀参数｜MCL inflation parameter |
| `--dna, -d` | — |  | 输入序列为DNA｜Input sequences are DNA |
| `--basic-only` | — |  | 仅执行基础分析｜Perform basic analysis only |
| `--generate-trees` | — |  | 生成系统发育树｜Generate phylogenetic trees |
| `--msa-program` | `mafft` | mafft/muscle | 多序列比对程序｜Multiple sequence alignment program |
| `--tree-program` | `fasttree` | fasttree/fasttree_fastest/raxml/raxml-ng/iqtree | 系统发育树推断程序｜Phylogenetic tree inference program |
| `--orthofinder-path` | `orthofinder` |  | OrthoFinder程序路径｜OrthoFinder program path |
| `--force` | — |  | 强制重新分析覆盖已有结果｜Force reanalysis overwriting existing results |
| `--skip-orthofinder` | — |  | 跳过OrthoFinder步骤直接分类｜Skip OrthoFinder step and go directly to classification |
| `--disable-rarefaction` | — |  | 禁用稀释曲线分析｜Disable rarefaction curve analysis |
| `--disable-single-copy` | — |  | 禁用单拷贝基因分析｜Disable single copy gene analysis |
| `--no-plots` | — |  | 不生成图表｜Do not generate plots |
| `--rarefaction-iterations` | `100` | int | 稀释分析迭代次数｜Rarefaction analysis iterations |
| `--single-copy-format` | `both` | by_orthogroup/by_genome/both | 单拷贝基因序列输出格式｜Single copy gene sequence output format |
| `--plot-format` | `png` | png/pdf/svg | 图表格式｜Plot format |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入蛋白质序列文件目录｜Input protein sequence files directory |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-n, --project-name` | — |  | 项目名称｜Project name |
| `--softcore-threshold` | `1` | int | Softcore基因缺失阈值｜Softcore gene missing threshold |
| `--dispensable-threshold` | `1` | int | Dispensable基因缺失阈值｜Dispensable gene missing threshold |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-s, --search, --search-program` | `blast` | blast/diamond/diamond_ultra_sens/mmseqs/blast_nucl | 序列搜索程序｜Sequence search program |
| `--mcl-inflation` | `1.2` | float | MCL inflation参数｜MCL inflation parameter |
| `--disable-single-copy` | — | store_true | 禁用单拷贝基因分析｜Disable single copy gene analysis |
| `--extract-sequences` | `True` | store_true | 提取单拷贝基因序列｜Extract single copy gene sequences |
| `--single-copy-format` | `both` | by_orthogroup/by_genome/both | 单拷贝基因序列输出格式｜Single copy gene sequence output format |
| `--disable-rarefaction` | — | store_true | 禁用稀释曲线分析｜Disable rarefaction curve analysis |
| `--rarefaction-iterations` | `100` | int | 稀释分析迭代次数｜Rarefaction analysis iterations |
| `-d, --dna` | — | store_true | 输入序列为DNA序列｜Input sequences are DNA |
| `--basic-only` | — | store_true | 仅进行基础分析｜Perform basic analysis only |
| `--generate-trees` | — | store_true | 生成系统发育树｜Generate phylogenetic trees |
| `--msa-program` | `mafft` | mafft/muscle | 多序列比对程序｜Multiple sequence alignment program |
| `--tree-program` | `fasttree` | fasttree/fasttree_fastest/raxml/raxml-ng/iqtree | 系统发育树构建程序｜Phylogenetic tree inference program |
| `--force` | — | store_true | 强制重新分析覆盖已有结果｜Force reanalysis overwriting existing results |
| `--skip-orthofinder` | — | store_true | 跳过OrthoFinder步骤直接进行分类｜Skip OrthoFinder step and go directly to classification |
| `--no-plots` | — | store_true | 不生成图表｜Do not generate plots |
| `--plot-format` | `png` | png/pdf/svg | 图表格式｜Plot format |
| `--orthofinder-path` | `~/miniforge3/envs/annot/bin/orthofinder` |  | OrthoFinder程序路径｜OrthoFinder program path |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先"认亲戚"(OrthoFinder 找同源群)，再"排家谱"(按出现频率分四类)，最后"画曲线/提序列/写报告"。

```text
输入目录(蛋白质 FASTA)
    │
    ▼
步骤1: 验证并清洗输入序列
    │
    ▼
步骤2: 运行 OrthoFinder 推断直系同源群 (可断点续传)
    │
    ▼
步骤3: 加载 Orthogroups.tsv / GeneCount.tsv
    │
    ▼
步骤4: 泛基因组四分类(core/softcore/dispensable/private)
    │
    ├─ 步骤5: 稀释曲线分析(可选,默认开)
    ├─ 步骤6: 单拷贝基因提取(可选,默认开)
    ├─ 步骤7: 可视化(可选,默认开)
    └─ 步骤8: 生成综合报告
```

## 输出 | Output

```text
pangenome_results/
├── 00_pipeline_info/
│   ├── software_versions.yml          # 软件版本与运行参数
│   └── pipeline_params.yml            # 流程参数快照
├── 01_orthofinder_analysis/
│   ├── orthofinder_results/           # OrthoFinder 原始结果(含 Orthogroups.tsv 等)
│   └── orthofinder_working_dir/       # OrthoFinder 工作目录
├── 02_pangenome_classification/
│   ├── pangenome_gene_families.txt    # 分类结果:Category/Orthogroup_ID/Gene_ID/Genome_Name
│   ├── pangenome_classification_summary.txt  # 分类统计摘要
│   ├── gene_families_detailed.tsv     # 基因家族详细表格(含是否单拷贝)
│   └── frequency_distribution_table.tsv     # 各类别的出现频率分布
├── 03_rarefaction_analysis/
│   ├── rarefaction_curve_detailed.tsv # 稀释曲线逐次抽样明细
│   └── rarefaction_curve_summary.tsv  # 稀释曲线各样本数均值/标准差
├── 04_single_copy_genes/
│   ├── by_orthogroup/                 # 每个单拷贝同源群一个 FASTA
│   ├── by_genome/                     # 每个基因组一个 FASTA
│   └── single_copy_genes_summary.txt  # 单拷贝基因汇总
├── 05_visualization/
│   ├── pangenome_classification_pie.png    # 分类饼图
│   ├── frequency_distribution.png          # 频率分布图
│   └── pangenome_rarefaction_curve.png     # 稀释曲线图
├── 06_reports/
│   └── pangenome_analysis_report.txt  # 综合分析报告
└── 99_logs/
    └── pangenome_analysis.log         # 运行日志
```

关键文件说明：

- **`pangenome_gene_families.txt`**：一行一个"同源群-基因-基因组"三元组，第一列是分类标签(core/softcore/dispensable/private)，是做下游功能分析的主表
- **`pangenome_classification_summary.txt`**：各类别的同源群数、基因数及占比，快速看泛基因组构成
- **`gene_families_detailed.tsv`**：与上面等价但多一列 `Is_Single_Copy`，方便程序化处理
- **`rarefaction_curve_summary.tsv`**：稀释曲线均值数据，画泛基因组"饱和曲线"用
- **`pangenome_analysis_report.txt`**：一页式综合报告，含结果解读指南

## 结果解读 | Interpreting Results

### 1. 分类构成（`pangenome_classification_summary.txt`）

**通俗理解|In plain words:** 看四类基因各占多少，就知道这个物种的泛基因组是"封闭"还是"开放"。

- **核心基因占比高**：泛基因组趋于封闭，群体基因集合稳定
- **非必需/私有基因占比高**：泛基因组开放，群体内基因多样性大（常见于细菌、真菌等）
- 单拷贝基因数是核心基因的子集，可作为跨基因组比较的"标准件"

### 2. 稀释曲线（`rarefaction_curve_summary.tsv`）

**通俗理解|In plain words:** Pan 曲线往上走=还在不断发现新基因；Core 曲线往下掉=共有基因随样本增多而减少。

- Pan 曲线趋平：再加样本也发现不了多少新基因，泛基因组已饱和
- Pan 曲线仍上升：还有大量"未探测"的私有/非必需基因，可考虑加样本
- Core 曲线下降后趋稳：核心基因组大小基本确定

### 3. 单拷贝基因（`04_single_copy_genes/`）

**通俗理解|In plain words:** 这些是"每个基因组里都恰好一份"的基因，最适合做系统发育树或多序列比对的标准输入。

- `by_orthogroup/`：每个同源群一个文件，序列名带基因组前缀，适合逐群比对建树
- `by_genome/`：每个基因组一个文件，适合做物种树输入或后续串联比对

## 参数选择建议 | Parameter Guidance

- **`-s/--search`**：基因组多(>50)或追求速度用 `diamond` 或 `mmseqs`；小规模或想复现经典结果用默认 `blast`
- **`--mcl-inflation`**：默认 1.2 通用；想分得更细(群更多)调大，想合并得更粗调小，范围 0.1~10
- **`--softcore-threshold`**：样本多时(>50)可放宽到 2~3，让"软核心"容纳偶发缺失；样本少(≤10)保持 1
- **`--rarefaction-iterations`**：默认 100；样本多时降采样组合数巨大，可降到 20~50 提速
- **`-d/--dna`**：只有输入是 DNA 序列时才加；蛋白质输入不要加
- **`--force`**：换了输入序列或想完全重算时才用，否则默认续跑更省时间

## 依赖 | Dependencies

- Python 3（pandas、numpy、matplotlib、pyyaml、biopython）
- OrthoFinder（conda 环境 `annot`，默认路径 `~/miniforge3/envs/annot/bin/orthofinder`）
- 若 `--generate-trees`：还需 mafft/muscle 与 fasttree/raxml/iqtree 等建树工具

## 常见问题 | FAQ

**Q1：重跑为什么没有重新算 OrthoFinder？**
默认开启断点续传：只要检测到有效的 `Orthogroups.tsv` + `GeneCount.tsv` 就直接复用。想完全重算加 `--force`。

**Q2：`--dispensable-threshold` 改了为什么分类没变？**
实际分类判定只用 `--softcore-threshold`（缺失数 ≤ 它=软核心；缺失更多且多个基因组共有=非必需）。`--dispensable-threshold` 目前只影响摘要/报告里的定义文字，不影响判定结果。

**Q3：输入是 DNA 序列怎么办？**
加 `-d/--dna`，程序会把序列类型记为 dna 并给 OrthoFinder 传 `-d`。

**Q4：单拷贝序列没有生成？**
确认没有用 `--disable-single-copy`，且输入序列 ID 与 OrthoFinder 结果里的基因 ID 能对上；ID 不一致会导致提取为空。

**Q5：输入目录里放了非 FASTA 文件会怎样？**
验证阶段会扫描 FASTA 扩展名；目录里找不到任何 FASTA 会直接报错退出。

**Q6：生成的图想换成 PDF/SVG？**
用 `--plot-format pdf` 或 `--plot-format svg`（默认 png）。
