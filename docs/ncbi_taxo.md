# NCBI分类学注释工具|NCBI Taxonomy Annotation Tool

**从BLAST结果获取NCBI分类学注释和统计信息 | Get NCBI Taxonomy Annotations and Statistics from BLAST Results**

## 功能概述 | Overview

NCBI分类学注释工具是一个自动化流程，从BLAST比对结果或accession列表出发，获取完整的NCBI分类学注释信息，并生成详细的统计报告。支持多层级统计（属、种等），可分析BLAST hit次数分布和唯一accession多样性。

## 主要特性 | Key Features

- **📊 灵活输入**: 支持BLAST结果文件或accession ID列表
- **🔍 自动检测**: 自动识别输入文件类型
- **📈 多维统计**: 同时统计BLAST hit次数和唯一accession数量
- **🎯 多层级分析**: 支持Kingdom到Species的多层级统计
- **📄 多种输出**: TXT/CSV格式输出，包含详细统计表格
- **⚙️ 高度可配置**: 自定义数据库路径、统计层级、输出格式等

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 从BLAST结果获取分类学注释和统计
biopytools ncbi-taxo -i blast_results.txt -o output_prefix

# 从accession列表获取
biopytools ncbi-taxo -i accessions.txt -o output --input-type accession

# 自定义统计层级（只统计属）
biopytools ncbi-taxo -i blast.txt -o result --stats-by genus

# 输出CSV格式统计
biopytools ncbi-taxo -i blast.txt -o result --stats-output csv
```

### 高级用法 | Advanced Usage

```bash
# 完整参数示例
biopytools ncbi-taxo \
    -i blast_results.txt \
    -o output_prefix \
    --input-type blast \
    --blast-column 2 \
    --taxid-db ~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz \
    --stats-by genus species \
    --stats-target both \
    --stats-output txt \
    --threads 4

# 禁用覆盖度统计（默认开启）
biopytools ncbi-taxo \
    -i blast_results.txt \
    -o output_without_coverage \
    --disable-coverage
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入文件（BLAST结果或accession列表） | `-i blast.txt` |
| `-o, --output-prefix` | 输出文件前缀 | `-o result_prefix` |

### 输入类型配置 | Input Type Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--input-type` | `auto` | 输入文件类型 (auto/blast/accession) |
| `--blast-column` | `2` | BLAST结果中accession所在列（从1开始）|

### 数据库配置 | Database Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--taxid-db` | `~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz` | TaxID数据库路径 |

### 分类学配置 | Taxonomy Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--lineage-format` | `{k};{p};{c};{o};{f};{g};{s}` | 分类层级格式字符串 |
| `--no-full-lineage` | `False` | 不保留完整lineage |

### 统计配置 | Statistics Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--stats-by` | `genus species` | 统计层级（可多选） |
| `--stats-target` | `both` | 统计对象 (blast_hits/unique_accessions/both) |
| `--stats-output` | `txt` | 统计输出格式 (txt/csv) |
| `--disable-coverage` | `False` | 禁用覆盖度统计（默认开启） |
| `--query-fasta` | `None` | Query序列文件（用于计算真实长度，推荐提供） |

### 性能配置 | Performance Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `4` | 线程数 |

## 输入文件格式 | Input File Formats

### BLAST结果文件 | BLAST Results File

标准12列BLAST输出格式（制表符分隔）：

```
ptg000021l	NC_009385.1	99.828	42994	6	42	3842	46784	42977	1	0.0	78923
ptg000021l	PP108240.1	94.946	19191	803	71	39772	58916	16765	35834	0.0	29911
```

**默认提取第2列的accession**。

### Accession列表文件 | Accession List File

纯文本文件，每行一个accession：

```
AB284571.1
AB284572.1
AB649206.1
AY003908.1
```

## 输出文件 | Output Files

### 输出文件列表

```
output_prefix.accessions.txt            # 唯一accession列表
output_prefix.acc2taxid.txt             # accession到taxid的映射
output_prefix.taxonomy.txt              # 完整分类学注释
output_prefix.statistics.txt            # 基本统计报告（TXT格式）
output_prefix.coverage_detail.txt       # 覆盖度详细统计（每个query的覆盖度分布）
output_prefix.coverage_summary.txt      # 覆盖度汇总统计（所有query的汇总）
output_prefix.ncbi_taxo.log             # 处理日志
```

### 覆盖度文件说明 | Coverage Files Description

#### coverage_detail.txt - 详细覆盖度统计
显示每个query序列在不同分类群的覆盖度百分比分布：
```
Query: ptg000021l
--------------------------------------------------------------------------------
  Phytophthora sojae                                            60.50%
  Pythium insidiosum                                            35.20%
  Other (< 1% each)                                              4.30%

Query: ptg000022l
--------------------------------------------------------------------------------
  Phytophthora infestans                                         95.00%
  Other (< 1% each)                                              5.00%
```

#### coverage_summary.txt - 汇总覆盖度统计
显示所有query的汇总统计信息：
```
GENUS 汇总统计
Taxon                                        Queries    Avg Coverage    Max Coverage
Phytophthora                                     150          45.23%          95.00%
Pythium                                           80          32.10%          85.50%
Other (< 1% avg)                                  20           2.15%           8.30%

整体统计
平均每query总覆盖度: 105.50%
最大query总覆盖度: 250.00%
最小query总覆盖度: 45.20%
```

### taxonomy.txt 格式

```txt
AB284571.1	135480	Eukaryota;Stramenopiles;Oomycota;Peronosporomycetes;Pythiales;Lagenidiaceae;Salilagenidium;Salilagenidium callinectes
AB284572.1	135484	Eukaryota;Stramenopiles;Oomycota;Peronosporomycetes;Pythiales;Lagenidiaceae;Salilagenidium;Salilagenidium thermophilum
AY003908.1	4787	Eukaryota;Stramenopiles;Oomycota;Peronosporomycetes;Peronosporales;Peronosporaceae;Phytophthora;Phytophthora infestans
```

### statistics.txt 格式（汇总表格）

```txt
================================================================================
NCBI分类学注释统计结果|NCBI Taxonomy Annotation Statistics
================================================================================

总体统计|Overall Statistics:
--------------------------------------------------------------------------------
总accession数量|Total unique accessions: 1643
总BLAST hit次数|Total BLAST hits: 1.57M
平均每accession hit次数|Average hits per accession: 956.78

================================================================================

GENUS 统计|GENUS Statistics
--------------------------------------------------------------------------------
Level      Level_Name                                         Count   Percentage
--------------------------------------------------------------------------------
Genus      Phytophthora                                        743         45.2%
Genus      Pythium                                            527         32.1%
Genus      Hyphochytrium                                      126          7.7%
Genus      Halophytophthora                                    98          6.0%
...

SPECIES 统计|SPECIES Statistics
--------------------------------------------------------------------------------
Level      Level_Name                                         Count   Percentage
--------------------------------------------------------------------------------
Species    Phytophthora sojae                                521         31.7%
Species    Phytophthora infestans                            156          9.5%
Species    Pythium insidiosum                                134          8.2%
...
```

## 使用示例 | Usage Examples

### 示例1：大豆疫霉菌BLAST结果分析

```bash
biopytools ncbi-taxo \
    -i 41_medium_quality.nt.blast.txt \
    -o soybean_blight_taxonomy \
    --stats-by genus species \
    --stats-target both
```

**结果解读**：
- 生成1643个唯一accession的分类学注释
- 统计157万次BLAST hit在各属、种的分布
- 识别出主要病原体：Phytophthora sojae (31.7%), P. infestans (9.5%)

### 示例2：只统计属层级

```bash
biopytools ncbi-taxo \
    -i blast_results.txt \
    -o genus_only \
    --stats-by genus
```

### 示例3：输出CSV格式用于后续分析

```bash
biopytools ncbi-taxo \
    -i blast.txt \
    -o analysis \
    --stats-output csv
```

生成 `analysis.statistics.csv`，可用Excel/R/Python进一步分析。

### 示例4：从已有accession列表获取注释

```bash
# 先提取accession
cut -f2 blast.txt | sort -u > my_accessions.txt

# 获取分类学注释
biopytools ncbi-taxo \
    -i my_accessions.txt \
    -o my_annotation \
    --input-type accession
```

### 示例5：覆盖度统计（默认开启）

```bash
# 分析query序列的比对覆盖度（默认开启，推荐提供query序列文件）
biopytools ncbi-taxo \
    -i blast_results.txt \
    -o coverage_analysis \
    --query-fasta query_sequences.fa \
    --stats-by genus species

# 如果不需要覆盖度统计，使用 --disable-coverage
biopytools ncbi-taxo \
    -i blast_results.txt \
    -o no_coverage \
    --disable-coverage
```

**输出内容**（启用覆盖度时）：
- 基本统计（accession数量、hit次数）
- 属/种统计（基于唯一accessions和blast hits）
- **覆盖度详细统计** (`coverage_detail.txt`)：每个query在不同taxon的覆盖分布
- **覆盖度汇总统计** (`coverage_summary.txt`)：所有query的总体覆盖度（各taxon占比）

**为什么需要 --query-fasta**：
- 如果不提供query序列文件，工具使用BLAST结果中的q.end来估算query长度
- 这可能导致不准确，特别是当query序列的某些区域没有比对时
- 提供FASTA文件可以获取query的真实长度，计算更准确的覆盖度

## 统计说明 | Statistics Description

### 统计层级 | Statistics Levels

| 层级 | 描述 | 示例 |
|------|------|------|
| `kingdom` | 界 | Eukaryota |
| `phylum` | 门 | Stramenopiles |
| `class` | 纲 | Oomycota |
| `order` | 目 | Peronosporales |
| `family` | 科 | Peronosporaceae |
| `genus` | 属 | Phytophthora |
| `species` | 种 | Phytophthora sojae |

### 统计对象 | Statistics Targets

- **`blast_hits`**: 统计BLAST hit次数分布，反映序列覆盖度
- **`unique_accessions`**: 统计唯一accession数量，反映物种多样性
- **`both`**: 同时输出两种统计（默认）

### 覆盖度统计 | Coverage Statistics

**默认开启**。工具会自动计算每个query序列比对到各accession的覆盖度百分比。

如需禁用覆盖度统计（可加快处理速度），使用 `--disable-coverage` 参数。

**重要**：推荐使用 `--query-fasta` 参数提供query序列文件，以获取准确的序列长度。

**计算方式**：
- **每个query的覆盖度**：`(该query比对到某taxon的碱基数 / 该query的总长度) × 100`
- **总体覆盖度统计**：所有query序列的比对区域中，属于每个taxon的百分比
  - 例如：所有query比对总长度1M bp，其中Phytophthora占600K bp（60%）
  - **注意**：由于同一query可能比对到多个taxon，各taxon百分比之和可能超过100%

**输出内容**：
1. **总体覆盖度**：各taxon在所有比对区域中的占比
2. **详细query覆盖度**：每个query在不同taxon的覆盖分布

**适用场景**：
- 评估query序列与不同物种的比对完整性
- 识别query序列的主要比对目标
- 检测query是否被多个物种覆盖（可能提示嵌合体或污染）

## 依赖要求 | Dependencies

### 必需工具 | Required Tools

- **taxonkit**: NCBI分类学工具包
  ```bash
  conda install -c bioconda taxonkit
  ```

- **zgrep**: gzip文件搜索工具（通常系统自带）

### 数据库要求 | Database Requirements

需要NCBI accession2taxid数据库：
- 默认路径：`~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz`
- 下载地址：https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: 提示"taxonkit: command not found"**
```bash
# 安装taxonkit
conda install -c bioconda taxonkit
```

**Q: 提示"TaxID数据库不存在"**
```bash
# 检查数据库路径
ls -lh ~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz

# 或指定正确的路径
biopytools ncbi-taxo -i blast.txt -o result --taxid-db /path/to/nucl_gb.accession2taxid.gz
```

**Q: BLAST结果中没有第2列**
```bash
# 使用--blast-column指定正确的列
biopytools ncbi-taxo -i blast.txt -o result --blast-column 3
```

**Q: 很多accession找不到taxid**
- 检查accession格式是否正确（应包含版本号，如NC_123456.1）
- 确认使用的数据库类型正确（核酸数据库nucl_gb）

## 注意事项 | Important Notes

1. **Accession版本号**: accession应包含版本号（如NC_123456.1），不带版本号可能导致查找失败
2. **数据库大小**: nucl_gb.accession2taxid.gz约2.4GB，需要足够磁盘空间
3. **处理时间**: zgrep查找taxid是耗时步骤，大量accession（>10万）可能需要较长时间
4. **内存使用**: 默认配置可处理约10万accession，更多accession建议分批处理

## 引用 | Citation

如果在研究中使用此工具，请引用相关的NCBI Taxonomy数据库和taxonkit工具：

```
NCBI Resource Coordinators. (2024).
Database resources of the National Center for Biotechnology Information.
Nucleic Acids Research, 52(D1), D18–D28.

Shen, W., Ren, H., & Xue, Z. (2023).
TaxonKit: a cross-platform and efficient NCBI taxonomy toolkit.
Bioinformatics, 39(12), btad689.
```
