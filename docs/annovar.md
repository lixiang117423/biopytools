# 📝 ANNOVAR 变异注释分析模块

**专业的基因变异功能注释工具 | Professional Gene Variant Functional Annotation Tool**

## 📖 功能概述 | Overview

ANNOVAR 变异注释分析模块是一个强大的基因变异功能注释工具，基于ANNOVAR软件构建，提供从GFF3注释文件处理到最终变异功能注释的完整流程。支持自动化的多步骤处理、灵活的质量控制和详细的注释结果输出，适用于各种基因组变异分析研究。

## ✨ 主要特性 | Key Features

- **🔄 完整注释流程**: GFF3转换→序列提取→VCF处理→变异注释四步骤自动化
- **🎯 灵活步骤控制**: 支持单独运行任意步骤或完整流程
- **🛡️ 智能格式处理**: 自动GFF3格式清理和修复功能
- **🔍 质量控制过滤**: 可配置的VCF质量阈值过滤
- **📂 多格式支持**: 支持标准GFF3、FASTA、VCF格式文件
- **⚙️ 高度可配置**: 自定义软件路径、数据库路径和输出设置
- **📊 详细日志记录**: 完整的处理过程日志和错误追踪
- **🚀 高效处理**: 优化的处理流程，支持大规模基因组数据

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本变异注释分析
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b OV \
    -o annotation_results

# 使用自定义质量阈值
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b OV \
    --qual-threshold 30 \
    -o high_quality_results
```

### 高级用法 | Advanced Usage

```bash
# 自定义软件和数据库路径的完整分析
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b OV \
    -a /path/to/annovar \
    -d /path/to/annovar_db \
    --enable-vcf-filter \
    --qual-threshold 25 \
    -o custom_annotation

# 只运行特定步骤
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b OV \
    --step 1 \
    -o gff_conversion_only
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --gff3` | GFF3注释文件路径 | `-g annotation.gff3` |
| `-f, --genome` | 基因组序列文件路径 | `-f genome.fa` |
| `-v, --vcf` | VCF变异文件路径 | `-v variants.vcf` |
| `-b, --build-ver` | 基因组构建版本标识符 | `-b OV` 或 `-b KY131` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-a, --annovar-path` | `/share/org/YZWL/yzwl_lixg/software/annovar/annovar` | 🛠️ ANNOVAR软件安装路径 |
| `-d, --database-path` | `./database` | 💾 ANNOVAR数据库路径 |
| `-o, --output-dir` | `./annovar_output` | 📁 输出目录路径 |

### 质量控制配置 | Quality Control Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-q, --qual-threshold` | `20` | 🎯 VCF质量过滤阈值 |
| `--enable-vcf-filter` | `False` | 🔍 启用VCF过滤步骤（默认跳过） |

### 处理控制选项 | Processing Control Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-s, --step` | `全部` | 🎯 只运行指定步骤 (1/2/3/4) |
| `--skip-gff-cleaning` | `False` | ⏭️ 跳过GFF3文件的格式清理 |
| `--skip-gff-fix` | `False` | ⏭️ 跳过GFF3文件的自动修复 |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 🔄 GFF3转换 | 将GFF3注释文件转换为ANNOVAR兼容格式 |
| **2** | 🧬 序列提取 | 从基因组序列中提取基因和转录本序列 |
| **3** | 🔍 VCF处理 | VCF文件格式化和质量过滤 |
| **4** | 📝 变异注释 | 基于基因注释进行变异功能注释 |

## 📁 输入文件格式 | Input File Formats

### GFF3注释文件 | GFF3 Annotation File

标准GFF3格式注释文件，包含基因结构信息：

```gff3
##gff-version 3
##sequence-region chromosome1 1 50000000
chromosome1	RefSeq	gene	1000	5000	.	+	.	ID=gene1;Name=GENE1
chromosome1	RefSeq	mRNA	1000	5000	.	+	.	ID=transcript1;Parent=gene1
chromosome1	RefSeq	exon	1000	1500	.	+	.	ID=exon1;Parent=transcript1
chromosome1	RefSeq	CDS	1200	1400	.	+	0	ID=cds1;Parent=transcript1
```

**文件要求**:
- 标准GFF3格式（版本3）
- 包含基因、转录本、外显子、CDS信息
- 正确的层级关系（Parent-Child）

### 基因组序列文件 | Genome Sequence File

标准FASTA格式的基因组序列：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### VCF变异文件 | VCF Variant File

标准VCF格式的变异调用结果：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chromosome1,length=50000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chromosome1	1250	.	A	G	45.2	PASS	DP=30	GT:DP	1/1:30
chromosome1	2340	.	T	C	62.8	PASS	DP=25	GT:DP	0/1:25
```

## 💡 使用示例 | Usage Examples

### 示例1：完整注释流程 | Example 1: Complete Annotation Pipeline

```bash
# 对植物基因组变异进行完整注释
biopytools annovar \
    -g Arabidopsis.gff3 \
    -f Arabidopsis_genome.fa \
    -v SNP_variants.vcf \
    -b AT \
    -o arabidopsis_annotation \
    --qual-threshold 25
```

### 示例2：微生物基因组注释 | Example 2: Microbial Genome Annotation

```bash
# 微生物基因组变异注释
biopytools annovar \
    -g bacteria.gff3 \
    -f bacteria_genome.fasta \
    -v mutations.vcf \
    -b BACT001 \
    -a /opt/annovar \
    -d /data/annovar_db \
    --enable-vcf-filter \
    --qual-threshold 30 \
    -o bacterial_variants
```

### 示例3：分步骤处理 | Example 3: Step-by-Step Processing

```bash
# 第一步：GFF3格式转换
biopytools annovar \
    -g raw_annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b GEN01 \
    --step 1 \
    -o step1_gff_conversion

# 第二步：序列提取
biopytools annovar \
    -g clean_annotation.gff3 \
    -f genome.fa \
    -v variants.vcf \
    -b GEN01 \
    --step 2 \
    --skip-gff-cleaning \
    -o step2_seq_extract

# 第三步：VCF处理
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v raw_variants.vcf \
    -b GEN01 \
    --step 3 \
    --enable-vcf-filter \
    --qual-threshold 20 \
    -o step3_vcf_process

# 第四步：变异注释
biopytools annovar \
    -g annotation.gff3 \
    -f genome.fa \
    -v clean_variants.vcf \
    -b GEN01 \
    --step 4 \
    -o step4_annotation
```

### 示例4：跳过清理步骤的快速处理 | Example 4: Fast Processing with Skip Cleaning

```bash
# 使用已清理的GFF3文件进行快速注释
biopytools annovar \
    -g clean_annotation.gff3 \
    -f reference_genome.fa \
    -v filtered_variants.vcf \
    -b REF01 \
    --skip-gff-cleaning \
    --skip-gff-fix \
    -o fast_annotation
```

### 示例5：高质量严格过滤分析 | Example 5: High-Quality Strict Filtering Analysis

```bash
# 高质量变异的严格过滤和注释
biopytools annovar \
    -g comprehensive.gff3 \
    -f genome_v2.fa \
    -v high_conf_variants.vcf \
    -b HG002 \
    -a /usr/local/annovar \
    -d /database/annovar_humandb \
    --enable-vcf-filter \
    --qual-threshold 50 \
    -o high_quality_annotation
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
annovar_output/
├── 1_gff_conversion/           # GFF3转换结果
│   ├── cleaned_annotation.gff3
│   ├── gene_structure.txt
│   └── conversion.log
├── 2_sequence_extraction/      # 序列提取结果
│   ├── gene_sequences.fa
│   ├── transcript_sequences.fa
│   ├── protein_sequences.fa
│   └── extraction.log
├── 3_vcf_processing/          # VCF处理结果
│   ├── filtered_variants.vcf
│   ├── annovar_input.avinput
│   ├── vcf_statistics.txt
│   └── processing.log
├── 4_annotation_results/      # 变异注释结果
│   ├── variants.variant_function        # 基因区域功能
│   ├── variants.exonic_variant_function # 外显子变异功能
│   ├── variants.gene_based.txt         # 基于基因的注释
│   ├── annotation_summary.txt          # 注释结果汇总
│   └── annotation.log
├── database/                  # 本地数据库文件
│   ├── build_ver_gene.txt
│   ├── build_ver_geneMrna.fa
│   └── build_ver_refGene.txt
└── logs/                     # 完整日志记录
    ├── annovar_analysis.log
    ├── error.log
    └── step_summary.log
```

### 关键输出文件说明 | Key Output Files Description

- **variant_function**: 变异位点的基因组区域功能分类
- **exonic_variant_function**: 外显子区域变异的详细功能注释
- **gene_based.txt**: 基于基因模型的变异影响预测
- **annotation_summary.txt**: 所有变异的注释结果汇总表
- **filtered_variants.vcf**: 质量过滤后的VCF文件
- **annovar_input.avinput**: ANNOVAR标准输入格式文件

### 注释结果分类 | Annotation Result Categories

| 功能分类 | 描述 | 重要性 |
|----------|------|--------|
| **exonic** | 外显子区域变异 | 🔴 高 |
| **splicing** | 剪切位点变异 | 🔴 高 |
| **ncRNA** | 非编码RNA区域 | 🟡 中 |
| **UTR5/UTR3** | 5'/3'非翻译区 | 🟡 中 |
| **intronic** | 内含子区域 | 🟢 低 |
| **intergenic** | 基因间区域 | 🟢 低 |

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **ANNOVAR** (版本 2020-06-08 或更新)
  - 下载地址: http://annovar.openbioinformatics.org/
  - 需要注册和许可证
- **Perl** (版本 5.10+)
  - ANNOVAR的运行环境
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理
  - `subprocess` - 系统调用

### 安装依赖软件 | Installing Dependencies

```bash
# 安装ANNOVAR (需要注册下载)
# 1. 从官网下载: http://annovar.openbioinformatics.org/
# 2. 解压并配置环境变量
tar -xzf annovar.latest.tar.gz
export PATH=$PATH:/path/to/annovar

# 安装Perl (如果尚未安装)
# Ubuntu/Debian
sudo apt-get install perl

# CentOS/RHEL  
sudo yum install perl

# 安装Python包
pip install click pathlib
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大基因组推荐16GB以上）
- **存储**: 预留基因组文件大小5倍的磁盘空间
- **网络**: 如需下载ANNOVAR数据库，建议稳定网络

## ⚠️ 注意事项 | Important Notes

1. **许可证要求**: ANNOVAR需要学术许可证，请确保合规使用
2. **基因组版本**: build-ver参数需与基因组版本严格对应
3. **文件格式**: 确保输入文件符合标准格式要求
4. **路径设置**: 软件路径和数据库路径需正确配置
5. **内存使用**: 大基因组文件可能需要大量内存

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "perl: command not found" 错误**
```bash
# 安装Perl
sudo apt-get install perl  # Ubuntu/Debian
sudo yum install perl      # CentOS/RHEL
```

**Q: "ANNOVAR not found" 错误**
```bash
# 检查ANNOVAR安装
which annotate_variation.pl
# 如未找到，请检查--annovar-path参数
biopytools annovar ... --annovar-path /correct/path/to/annovar
```

**Q: GFF3格式错误**
```bash
# 启用GFF3清理和修复
biopytools annovar ... # 默认会自动清理

# 或手动检查GFF3文件格式
grep "##gff-version 3" annotation.gff3
```

**Q: 内存不足错误**
```bash
# 检查系统内存
free -h

# 考虑分步处理大文件
biopytools annovar ... --step 1  # 先转换GFF3
biopytools annovar ... --step 2  # 再提取序列
```

**Q: VCF文件质量问题**
```bash
# 启用VCF过滤并降低阈值
biopytools annovar ... --enable-vcf-filter --qual-threshold 10

# 或跳过VCF过滤步骤进行测试
biopytools annovar ... --step 4  # 只做注释
```

**Q: 找不到基因序列**
```bash
# 检查基因组文件和GFF3文件的染色体名称是否一致
grep "^>" genome.fa | head -5
grep "^#" annotation.gff3 | grep sequence-region
```

## 📊 结果解读指南 | Result Interpretation Guide

### 变异功能优先级 | Variant Function Priority

1. **高优先级** 🔴
   - 外显子变异 (frameshift, stopgain, stoploss)
   - 剪切位点变异 (splicing)

2. **中优先级** 🟡
   - 同义/非同义变异 (synonymous/nonsynonymous)
   - UTR区变异
   - ncRNA变异

3. **低优先级** 🟢
   - 内含子变异
   - 基因间区变异

### 注释结果筛选建议 | Annotation Result Filtering Suggestions

```bash
# 筛选高影响变异
grep -E "(exonic|splicing)" variants.variant_function

# 筛选蛋白质改变变异
grep -E "(frameshift|stopgain|stoploss|nonsynonymous)" variants.exonic_variant_function
```

## 📚 相关资源 | Related Resources

- [ANNOVAR官方文档](http://annovar.openbioinformatics.org/en/latest/)
- [GFF3格式规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [VCF格式规范](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [变异功能注释最佳实践](https://www.nature.com/articles/nrg3933)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: ANNOVAR软件本身需要单独的学术许可证，请访问官网获取。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用ANNOVAR相关文献：

```
Wang, K., Li, M., & Hakonarson, H. (2010). 
ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. 
Nucleic acids research, 38(16), e164-e164.
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | VCF变异文件｜VCF variant file path |
| `--gff3, -g` | 必填 |  | GFF3注释文件｜GFF3 annotation file path |
| `--genome, -G` | 必填 |  | 基因组序列文件｜Genome sequence file path |
| `--build-ver, -b` | 必填 |  | 基因组版本(例如: OV, KY131)｜Genome build version identifier (e.g., OV, KY131) |
| `--annovar-path, -a` | `~/software/annovar/annovar` |  | ANNOVAR软件路径｜ANNOVAR software installation path |
| `--output-dir, -o` | `./annovar_output` | Path | 输出目录(同时作为ANNOVAR数据库目录)｜Output directory (also the ANNOVAR database dir) |
| `--qual-threshold, -q` | `20` | int | VCF质量阈值｜VCF quality filtering threshold |
| `--step, -s` | — | IntRange | 运行指定步骤｜Run only specified step: 1: GFF3转换｜GFF3 conversion 2: 提取序列｜Extract sequences 3: VCF处理｜VCF processing 4: 变异注释｜Variant annotation |
| `--skip-gff-cleaning` | — |  | 跳过GFF3清理｜Skip GFF3 file format cleaning |
| `--skip-gff-fix` | — |  | 跳过GFF3修复｜Skip automatic GFF3 file fixes |
| `--enable-vcf-filter` | — |  | 启用VCF过滤(默认跳过)｜Enable VCF filtering step (skipped by default) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF变异文件路径｜VCF variant file path |
| `-g, --gff3` | 必填 |  | GFF3注释文件路径｜GFF3 annotation file path |
| `-G, --genome` | 必填 |  | 基因组序列文件路径｜Genome sequence file path |
| `-b, --build-ver` | 必填 |  | 基因组构建版本标识符(如: OV, KY131) - 不应包含路径分隔符｜Genome build version identifier (e.g., OV, KY131) - should not contain path separators |
| `-a, --annovar-path` | `~/software/annovar/annovar` |  | ANNOVAR软件安装路径｜ANNOVAR software installation path |
| `-o, --output-dir` | `./annovar_output` |  | 输出目录(同时作为ANNOVAR数据库目录,refGene在此生成)｜Output directory (also the ANNOVAR database dir; refGene is built here) |
| `-q, --qual-threshold` | `20` | int | VCF质量过滤阈值(仅在启用VCF过滤时生效)｜VCF quality filtering threshold (only effective when VCF filtering is enabled) |
| `-s, --step` | — | 1/2/3/4 | 只运行指定步骤｜Run only specified step (1: gff3转换｜gff3 conversion, 2: 提取序列｜extract sequences, 3: VCF处理｜VCF processing, 4: 注释｜annotation) |
| `--skip-gff-cleaning` | — | store_true | 跳过GFF3文件的格式清理(attributes清理和坐标修复)｜Skip GFF3 file format cleaning (attributes cleaning and coordinate fixing) |
| `--skip-gff-fix` | — | store_true | 跳过GFF3文件的自动修复(CDS phase等问题)｜Skip automatic GFF3 file fixes (CDS phase and other issues) |
| `--skip-vcf-filter` | `True` | store_true | 跳过VCF过滤步骤，直接使用输入的VCF文件(默认启用)｜Skip VCF filtering step, use input VCF file directly (enabled by default) |
| `--enable-vcf-filter` | — | store_true | 启用VCF过滤步骤(使用bcftools)｜Enable VCF filtering step (using bcftools) |

<!-- END PARAMS:auto -->
