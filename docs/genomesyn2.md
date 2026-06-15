# GenomeSyn2 比较基因组学可视化工具

**GenomeSyn2 比较基因组学可视化工具 | GenomeSyn2 Comparative Genomics Visualization Tool**

## 功能概述 | Overview

GenomeSyn2 是一个专业的比较基因组学可视化平台，支持大规模基因组共线性分析、结构变异可视化和祖先血统解析。提供基因组、染色体和基因三个尺度的可视化功能，适用于多基因组比较、泛基因组分析和种群遗传研究。

## 主要特性 | Key Features

- **多序列比对支持**: 支持MUMmer、Minimap2、BLASTp、MMseqs2、DIAMOND等多种比对工具
- **多种可视化模式**: 基因组共线性图、局部基因结构视图、祖先血统解析视图
- **灵活的注释系统**: 支持BED、GFF3格式，多种可视化样式
- **完整分析流程**: 从序列比对到可视化的端到端解决方案
- **高度可配置**: 自定义软件路径、输出设置和绘图参数

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基因组比对分析
biopytools genomesyn2 \
    --align mummer \
    --genome ./genome_dir/ \
    --outdir ./output/ \
    --threads 12

# 蛋白质比对分析
biopytools genomesyn2 \
    --align blastp \
    --genome ./genome_dir/ \
    --gene ./gene_dir/ \
    --outdir ./output/

# 从VCF计算SNP密度和一致性
biopytools genomesyn2 \
    --vcf variants.vcf \
    --bin 50000

# 绘制祖先血统解析图
biopytools genomesyn2 \
    --identity SNP_identity.bed \
    --density SNP_density.bed

# 绘制共线性图
biopytools genomesyn2 --conf total.conf
```

## 参数说明 | Parameters

### 比对模式参数 | Alignment Mode Parameters

| 参数 | 类型 | 描述 |
|------|------|------|
| `--align` | choice | 比对软件类型 (mummer/minimap2/blastp/mmseqs/diamond) |
| `--genome` | path | 基因组文件目录（文件名需按数字排序，如1.genome1.fa, 2.genome2.fa） |
| `--gene` | path | 基因注释文件目录（蛋白质比对时必需） |
| `--outdir` | path | 输出目录 |

**文件命名规范**：
- 基因组文件：`1.MH63RS3.fasta`, `2.T.mark.fasta`, `3.Y.mark.fasta`
- 基因注释文件：`1.MH63.gene.gff3`, `2.T.gene.gff3`, `3.Y.gene.gff3`

### VCF模式参数 | VCF Mode Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--vcf` | - | VCF文件路径（需包含染色体长度信息） |
| `--bin` | 50000 | SNP分析的bin大小 |
| `--identity` | - | SNP一致性BED文件（已预计算） |
| `--density` | - | SNP密度BED文件（已预计算） |

### 绘图模式参数 | Plotting Mode Parameters

| 参数 | 描述 |
|------|------|
| `--conf` | 配置文件路径 |
| `--anno` | 显示注释配置选项 |

### 文件生成模式参数 | File Generation Mode Parameters

| 参数 | 描述 |
|------|------|
| `--type` | 文件类型 (fa/prot/anno) |
| `--path` | 文件路径 |
| `--out` | 输出文件名 |

### 通用参数 | Common Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | 12 | 线程数 |

## 使用场景 | Usage Scenarios

### 场景1: 多基因组共线性分析 | Scenario 1: Multi-Genome Synteny Analysis

```bash
# 使用MUMmer进行全基因组比对
biopytools genomesyn2 \
    --align mummer \
    --genome ./genome_data/ \
    --outdir ./mummer_results/ \
    --threads 24

# 生成配置文件模板
biopytools genomesyn2 --conf ? > anno.conf

# 编辑配置文件后绘制共线性图
biopytools genomesyn2 --conf anno.conf
```

### 场景2: 蛋白质序列比较 | Scenario 2: Protein Sequence Comparison

```bash
# 使用DIAMOND进行快速蛋白质比对
biopytools genomesyn2 \
    --align diamond \
    --genome ./genome_data/ \
    --gene ./gene_data/ \
    --outdir ./diamond_results/ \
    --threads 24
```

### 场景3: 祖先血统解析 | Scenario 3: Ancestry Deconvolution

```bash
# 步骤1: 从VCF计算SNP统计
biopytools genomesyn2 \
    --vcf parents.progeny.snps.genotype.vcf \
    --bin 50000

# 步骤2: 使用生成的文件绘制血统解析图
biopytools genomesyn2 \
    --identity ./SNP_identity.50Kb.bed \
    --density ./SNP_density.50Kb.bed
```

### 场景4: 局部基因结构可视化 | Scenario 4: Local Gene Structure Visualization

在配置文件中设置：
```ini
[show_region]
region = MH63:Chr10:24,850,000-24,885,000
gene_list = gene.info.tsv
```

### 场景5: 生成文件信息表 | Scenario 5: Generate File Info Table

```bash
# 生成基因组文件列表
biopytools genomesyn2 \
    --type fa \
    --path ./genome_data/ \
    --out genome.info.tsv

# 生成基因注释文件列表
biopytools genomesyn2 \
    --type anno \
    --path ./annotation_data/ \
    --out anno.info.tsv
```

## 配置文件格式 | Configuration File Format

配置文件采用INI格式，包含以下部分：

### [genome_info] - 基因组信息
```ini
gonomes_filetype = bed    # 文件类型 (fasta/bed)
gonomes_list = chr_length.info.tsv  # 染色体长度信息文件
```

### [synteny_info] - 共线性信息
```ini
line_type = curve         # 连接线样式 (curve/line)
synteny_list = synteny.info.tsv  # 共线性信息文件
```

### [save_info] - 保存设置
```ini
figure_type = pdf         # 图形格式 (svg/pdf/png)
savefig1 = GenomeSyn2.figure1.pdf
savefig2 = GenomeSyn2.figure2.pdf
```

### [centromere_info] - 着丝粒信息
```ini
centromere_list = centromere.info.tsv
```

### [telomere_info] - 端粒信息
```ini
telomere_list = telomere.info.tsv
telomere_color = #441680
opacity = 100%
```

### [anno_info] - 注释信息
```ini
anno_number = [1,2,3]
anno_name = [PAV,SNP,GC Content]
anno_color = ['#5FB6DE','#0000FF','#000000']
anno_type = [rectangle,barplot,lineplot]
anno_position = [top,top,bottom]
anno_height = [5,5,5]
min_max_value = [normal,auto,0.4:0.5]
anno_window = [none,none,100000]
opacity = [50%,100%,100%]
file_type = [bed,bed,gff3]
filter_type = [none,none,none]
anno_list = [PAV.info.tsv,SNP.info.tsv,GC.info.tsv]
```

## 输出结果 | Output Results

### 比对模式输出 | Alignment Mode Output

```
outdir/
├── fa_bed/                          # 处理后的基因组文件
├── mummer/ (或minimap2/等)           # 比对结果
│   ├── 1.genome1_vs_genome2.delta.filter.tsv
│   └── ...
├── chr_length.info.tsv              # 染色体长度信息
├── genomes.info.tsv                 # 基因组信息
├── synteny.info.tsv                 # 共线性信息
├── total.conf                       # 配置文件
├── GenomeSyn2.figure1.pdf           # 单染色体共线性图
└── GenomeSyn2.figure2.pdf           # 多染色体共线性图
```

### VCF模式输出 | VCF Mode Output

```
./SNP_identity.50Kb.bed              # SNP一致性统计
./SNP_density.50Kb.bed               # SNP密度统计
./GenomeSyn2.50Kb.pdf                # 血统解析图
```

## VCF文件格式要求 | VCF File Format Requirements

VCF文件必须包含：

1. **染色体长度信息**（必需）：
```
##contig=<ID=Chr01,length=45027022>
##contig=<ID=Chr02,length=37301368>
```

2. **样本颜色信息**（可选）：
```
##color=<Sample=Sample1,color="#39A5D6">
##color=<Sample=Sample2,color="#43A98C">
```

3. **SNP基因型数据**：
```
#CHROM  POS    ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
Chr01   5741   .   G    T    192.885  PASS  *    *      0/0      1/1
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Perl** 5.32+ (带BioPerl和SVG模块)
- **Python** 3.8+
- **比对工具** (根据需要选择)：
  - MUMmer4
  - Minimap2
  - BLAST+
  - DIAMOND
  - MMseqs2

### Python包依赖 | Python Package Dependencies

```bash
pip install click
```

## 故障排除 | Troubleshooting

### Q1: "Can't locate Bio/SeqIO.pm" 错误

确保安装了BioPerl：
```bash
# 通过conda安装
conda install -c bioconda perl-bioperl-core
```

### Q2: 比对软件未找到

检查软件是否在PATH中或使用完整路径：
```bash
which nucmer  # 检查MUMmer
which minimap2  # 检查Minimap2
```

### Q3: VCF文件格式错误

确保VCF文件头部包含染色体长度信息：
```bash
grep "##contig=" your.vcf
```

### Q4: 文件命名不符合要求

基因组文件必须按数字前缀命名：
```
正确|Correct:
├── 1.genome1.fa
├── 2.genome2.fa
└── 3.genome3.fa

错误|Wrong:
├── genome1.fa
├── genome2.fa
└── genome3.fa
```

## 相关资源 | Related Resources

- [GenomeSyn2官方文档](https://github.com/banzhou59/GenomeSyn2)
- [BioConda基因组_syn2](https://anaconda.org/bioconda/genomesyn2)
- [比较基因组学最佳实践](https://www.ncbi.nlm.nih.gov/pmc/articles/PMCPMC3098275/)

## 许可证 | License

本项目采用MIT许可证。

注意：GenomeSyn2软件本身遵循其原始许可证。

## 引用信息 | Citation

如果在研究中使用此工具，请引用GenomeSyn2相关文献：

```
Zhou, Z., Zhao, H., Chai, Y., Zhao, R., Qian, Y., Zhong, Y., Shao, Y., Chen, L., Song, J., 2026.
GenomeSyn-II: a comparative genomics framework integrating synteny visualization.
J. Genet. Genomics.
https://doi.org/10.1016/j.jgg.2026.01.011
```
