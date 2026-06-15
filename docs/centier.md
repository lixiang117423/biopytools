# CentIER着丝粒鉴定模块|CentIER Centromere Identification Module

## 模块简介|Module Introduction

CentIER模块是对[CentIER](https://github.com/xxxCUIxxxx/CentIER)软件的封装，用于T2T(Telomere-to-Telomere)基因组组装的着丝粒识别和注释。

The CentIER module is a wrapper for the [CentIER](https://github.com/xxxCUIxxxx/CentIER) software, designed for centromere identification and annotation in T2T (Telomere-to-Telomere) assembled genomes.

## 功能特点|Features

- **多特征融合|Multi-feature Integration**: 结合K-mer频率分析、串联重复序列、LTR反转座子进行着丝粒预测|Combines k-mer frequency analysis, tandem repeats, and LTR retrotransposons for centromere prediction
- **Hi-C数据支持|Hi-C Data Support**: 可选使用Hi-C相互作用矩阵进行验证|Optional Hi-C interaction matrix for validation
- **高准确率|High Accuracy**: 对植物基因组着丝粒识别准确率>90%|>90% accuracy for plant genome centromere identification
- **多着丝粒支持|Multiple Centromeres Support**: 支持识别多着丝粒染色体|Supports identification of metapolycentric chromosomes
- **可视化输出|Visualization Output**: 自动生成着丝粒区域可视化图表|Automatically generates centromere region visualization

## 安装依赖|Installation

### 1. 安装CentIER软件|Install CentIER Software

```bash
# 克隆CentIER仓库|Clone CentIER repository
cd ~/software
git clone https://github.com/xxxCUIxxxx/CentIER.git
```

CentIER已自带以下工具|CentIER includes the following bundled tools:
- `hmmsearch` - HMMER序列搜索|HMMER sequence search
- `ltr_finder` - LTR反转座子查找|LTR retrotransposon finder
- `trf409.linux64` - 串联重复序列查找|Tandem repeat finder
- HMM数据库|HMM databases: REXdb.hmm, Ty3_gypsy.hmm

### 2. 可选外部工具|Optional External Tools

以下工具为可选，如未安装可跳过|The following tools are optional and can be skipped if not installed:

- **genometools (gt)**: 用于LTRharvest|For LTRharvest
- **LTR_retriever**: 用于LTR质量优化|For LTR quality refinement

如需安装|If needed:
```bash
# genometools可跳过（CentIER有内置的ltr_finder）|genometools can be skipped (CentIER has built-in ltr_finder)
# LTR_retriever安装|LTR_retriever installation
conda install -c bioconda ltr_retriever
```

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 最简单的用法|Simplest usage
biopytools centier -i genome.fa -o output_dir/

# 使用GFF注释|With GFF annotation
biopytools centier -i genome.fa -g annotation.gff3 -o output_dir/

# 使用Hi-C数据验证|With Hi-C data validation
biopytools centier -i genome.fa \
    --matrix1 hic_100k.mnd \
    --matrix2 hic_200k.mnd \
    --bed1 hic_100k.bed \
    --bed2 hic_200k.bed \
    -o output_dir/
```

### 参数说明|Parameters

| 参数|Parameter | 简写|Short | 类型|Type | 默认值|Default | 说明|Description |
|---------|---------|---------|---------|-------------|---------|
| `--input` | `-i` | 必需|required | str | - | 基因组FASTA文件|Genome FASTA file |
| `--output-dir` | `-o` | 可选|optional | str | `./centier_output` | 输出目录|Output directory |
| `--centier-path` | - | 可选|optional | str | `~/software/CentIER/CentIER-main` | CentIER软件路径|CentIER software path |
| `--gff` | - | 可选|optional | str | - | GFF/GTF注释文件|GFF/GTF annotation file |
| `--matrix1` | - | 可选|optional | str | - | Hi-C矩阵(100k)|Hi-C matrix at 100k resolution |
| `--matrix2` | - | 可选|optional | str | - | Hi-C矩阵(200k)|Hi-C matrix at 200k resolution |
| `--bed1` | - | 可选|optional | str | - | Hi-C BED(对应matrix1)|Hi-C BED for matrix1 |
| `--bed2` | - | 可选|optional | str | - | Hi-C BED(对应matrix2)|Hi-C BED for matrix2 |
| `--kmer-size` | `-k` | 可选|optional | int | `21` | K-mer大小|K-mer size |
| `--center-tolerance` | `-c` | 可选|optional | int | `15` | 中心容差|Center tolerance |
| `--step-len` | - | 可选|optional | int | `10000` | 分析步长|Analysis step length |
| `--mul-cents` | - | 可选|optional | flag | `False` | 保留所有潜在着丝粒|Retain all potential centromeres |
| `--mingap` | - | 可选|optional | int | `2` | 最小Gap值|Minimum gap value |
| `--signal-threshold` | - | 可选|optional | float | `0.7` | Hi-C信号阈值|Hi-C signal threshold |
| `--summary` | - | 可选|optional | flag | `False` | 输出结果摘要|Output result summary |

## 输出文件|Output Files

### 主要输出文件|Main Output Files

运行完成后，输出目录将包含以下文件|After completion, the output directory will contain:

| 文件名|Filename | 说明|Description |
|---------|---------|---------|
| `{prefix}_centromere_range.txt` | 着丝粒区域坐标|Centromere region coordinates |
| `{prefix}_all_centromere_seq.txt` | 着丝粒区域序列|Centromere region sequences |
| `{prefix}_monomer_seq.txt` | 单体序列|Monomer sequences |
| `{prefix}_monomer_in_centromere.txt` | 着丝粒内单体位置|Monomer positions in centromeres |
| `{prefix}_ltr_position.txt` | LTR反转座子位置|LTR retrotransposon positions |
| `{prefix}_LTR_statistics.txt` | LTR统计信息|LTR statistics |
| `{prefix}_draw_cen.svg` | 着丝粒可视化图表|Centromere visualization |
| `centier_summary.json` | 分析结果摘要（使用--summary时）|Analysis summary (when using --summary) |

### 输出目录结构|Output Directory Structure

```
centier_output/
├── 99_logs/
│   └── centier.log              # 运行日志|Run log
├── genome_centromere_range.txt  # 着丝粒区域|Centromere regions
├── genome_all_centromere_seq.txt  # 着丝粒序列|Centromere sequences
├── genome_ltr_position.txt      # LTR位置|LTR positions
├── genome_draw_cen.svg          # 可视化|Visualization
└── centier_summary.json         # 结果摘要|Result summary
```

## 分析原理|Analysis Principles

CentIER使用多种特征组合进行着丝粒预测|CentIER uses a combination of multiple features for centromere prediction:

### 1. K-mer频率分析|K-mer Frequency Analysis

- 滑动窗口计算基因组K-mer种类数|Calculates k-mer diversity in sliding windows
- 着丝粒区域K-mer种类显著降低|Centromeric regions show significantly reduced k-mer diversity
- 使用动态阈值识别候选区域|Uses dynamic thresholds to identify candidate regions

### 2. 串联重复序列|Tandem Repeats

- 使用TRF (Tandem Repeats Finder)识别串联重复|Uses TRF to identify tandem repeats
- 着丝粒区域富含长串联重复序列|Centromeres are rich in long tandem repeats
- 分析重复单元大小和拷贝数|Analyzes repeat unit size and copy number

### 3. LTR反转座子|LTR Retrotransposons

- 使用LTR_Finder和LTRharvest识别LTR|Uses LTR_Finder and LTRharvest to identify LTRs
- 着丝粒区域富集特定类型LTR（如Ty3_gypsy）|Centromeres are enriched for specific LTR types (e.g., Ty3_gypsy)
- 使用HMMER进行LTR分类|Uses HMMER for LTR classification

### 4. Hi-C相互作用矩阵（可选）|Hi-C Interaction Matrix (Optional)

- 着丝粒区域在Hi-C矩阵中表现为信号减弱区|Centromeres appear as signal-depleted regions in Hi-C matrices
- 用于验证和精细化着丝粒边界|Used for validation and refinement of centromere boundaries

## 输出文件格式详解|Output File Format Details

### 1. centromere_range.txt

着丝粒区域坐标文件|Centromere region coordinates file:

```
chr1    1000000    2500000
chr2    500000     1800000
...
```

格式|Format: `染色体\t起始位置\t终止位置|Chromosome\tStart\tEnd`

### 2. monomer_seq.txt

单体序列文件|Monomer sequences file:

```
chr1    ATCGATCGATCGATCG...
chr2    GCTAGCTAGCTAGCTA...
...
```

格式|Format: `染色体\t单体序列|Chromosome\tMonomer Sequence`

### 3. LTR_statistics.txt

LTR统计文件|LTR statistics file:

```
ID    Ty1_copia    Ty3_gypsy    unknown    ...
chr1    5    12    3    ...
chr2    3    8    1    ...
...
```

## 示例工作流|Example Workflow

### 完整工作流（带Hi-C验证）|Complete Workflow (with Hi-C validation)

```bash
# 1. 准备输入文件|Prepare input files
# - 基因组: genome.fa
# - 注释: annotation.gff3 (可选)
# - Hi-C数据: hic_100k.mnd, hic_200k.mnd, hic_100k.bed, hic_200k.bed

# 2. 运行CentIER|Run CentIER
biopytools centier -i genome.fa \
    -g annotation.gff3 \
    --matrix1 hic_100k.mnd \
    --matrix2 hic_200k.mnd \
    --bed1 hic_100k.bed \
    --bed2 hic_200k.bed \
    -o centier_results/ \
    --summary

# 3. 查看结果|View results
cat centier_results/genome_centromere_range.txt
cat centier_results/centier_summary.json

# 4. 可视化|Visualization
# 在浏览器中打开SVG文件|Open SVG file in browser
firefox centier_results/genome_draw_cen.svg
```

### 快速工作流（无Hi-C）|Quick Workflow (without Hi-C)

```bash
# 仅使用基因组序列|Using genome sequence only
biopytools centier -i genome.fa -o centier_results/ --summary
```

## 常见问题|FAQ

### Q1: Hi-C数据格式要求是什么？|What are the Hi-C data format requirements?

A: CentIER需要Juicer工具输出的mnd格式矩阵文件和对应的BED文件：
- `matrix1`: 100kb分辨率的Hi-C矩阵|Hi-C matrix at 100kb resolution
- `matrix2`: 200kb分辨率的Hi-C矩阵|Hi-C matrix at 200kb resolution
- `bed1`: 对应matrix1的BED文件|BED file corresponding to matrix1
- `bed2`: 对应matrix2的BED文件|BED file corresponding to matrix2

### Q2: 是否必须提供GFF注释文件？|Is GFF annotation file mandatory?

A: 不是。GFF注释文件是可选的，用于提高预测精度。如果没有注释文件，CentIER仍可正常工作。

### Q3: 什么是mul_cents参数？|What is the mul_cents parameter?

A: 某些物种的染色体可能具有多个着丝粒区域（多着丝粒染色体）。默认情况下，CentIER只保留最可能的着丝粒。使用`--mul-cents`参数会保留所有潜在的着丝粒区域。

### Q4: 如何解读draw_cen.svg可视化图？|How to interpret the draw_cen.svg visualization?

A: SVG文件显示了每条染色体上的着丝粒区域：
- 水平线：染色体|Horizontal lines: chromosomes
- 红色区域：串联重复序列（单体）|Red regions: tandem repeats (monomers)
- 蓝色区域：LTR反转座子|Blue regions: LTR retrotransposons
- 不同颜色的LTR单体代表不同类型|Different colored LTR monomers represent different types

### Q5: 运行时间大概多久？|What is the approximate runtime?

A: 运行时间取决于基因组大小：
- 小基因组(<500 Mb): 1-2小时
- 中等基因组(500 Mb - 1.5 Gb): 2-6小时
- 大基因组(>1.5 Gb): 6-12小时

## 引用|Citation

如果使用CentIER，请引用以下论文|If you use CentIER, please cite:

> Xu, D. et al. CentIER: accurate centromere identification for plant genome. *Plant Communications* 101046 (2024) doi:10.1016/j.xplc.2024.101046

## 相关模块|Related Modules

- **assembly_qc**: 基因组组装质量评估|Genome assembly quality evaluation
- **busco**: 基因组完整性评估|Genome completeness assessment
- **merqury_qv**: 基于K-mer的QV评估|K-mer based QV evaluation

## 更新日志|Changelog

### v1.0.0 (2026-03-13)
- 初始版本发布|Initial release
- 封装CentIER v3.0.1|Wrapper for CentIER v3.0.1
- 支持Hi-C数据验证|Support for Hi-C data validation
- 添加Click CLI包装器|Added Click CLI wrapper
