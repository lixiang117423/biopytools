# 互作转录组分析工具

**双物种转录组同时分析工具 | Dual-species Transcriptome Simultaneous Analysis Tool**

## 功能概述 | Overview

互作转录组分析模块是一个专业的双物种转录组分析工具，能够同时对两个相互作用的物种进行转录组测序分析。该工具支持双参考基因组比对、智能reads物种分类、双物种表达定量和表达矩阵生成，适用于病原体-宿主互作、植物-微生物共生等多种研究场景。

## 主要特性 | Key Features

- **双基因组索引构建**: 分别构建两个物种的HISAT2索引，支持剪接位点和外显子信息
- **智能物种分类**: 基于MAPQ的reads物种分类策略，确保分类准确性
- **双物种定量**: 分别对两个物种的reads进行StringTie定量分析
- **表达矩阵生成**: 生成两个物种的表达矩阵文件，便于后续分析
- **灵活样本解析**: 支持目录扫描和样本信息文件两种输入方式
- **完整日志记录**: 详细的处理过程日志和分类统计信息

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
biopytools dual-rnaseq \
    --species1-name host \
    --species1-genome host.fa \
    --species1-gtf host.gtf \
    --species2-name pathogen \
    --species2-genome pathogen.fa \
    --species2-gtf pathogen.gtf \
    -i /data/fastq \
    -o dual_rnaseq_results
```

### 高级用法 | Advanced Usage

```bash
# 使用自定义参数
biopytools dual-rnaseq \
    --species1-name host \
    --species1-genome host.fa \
    --species1-gtf host.gtf \
    --species2-name pathogen \
    --species2-genome pathogen.fa \
    --species2-gtf pathogen.gtf \
    -i /data/fastq \
    -o results \
    -p "*.R1.fastq.gz" \
    -t 24 \
    --min-mapq 30
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `--species1-name` | 物种1名称 | `--species1-name host` |
| `--species1-genome` | 物种1基因组FASTA文件 | `--species1-genome host.fa` |
| `--species1-gtf` | 物种1 GTF注释文件 | `--species1-gtf host.gtf` |
| `--species2-name` | 物种2名称 | `--species2-name pathogen` |
| `--species2-genome` | 物种2基因组FASTA文件 | `--species2-genome pathogen.fa` |
| `--species2-gtf` | 物种2 GTF注释文件 | `--species2-gtf pathogen.gtf` |
| `-i, --input` | 输入FASTQ目录或样本信息文件 | `-i /data/fastq` |
| `-o, --output-dir` | 输出目录 | `-o results` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pattern` | `*_1.clean.fq.gz` | FASTQ文件命名模式 |
| `-t, --threads` | `12` | 线程数 |
| `--min-mapq` | `20` | 最小mapping quality值 |
| `--no-unique-only` | `False` | 不禁用非唯一比对 |

## 输入文件格式 | Input File Formats

### 基因组序列文件 | Genome Sequence File

标准FASTA格式的基因组序列：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### GTF注释文件 | GTF Annotation File

标准GTF格式注释文件：

```gtf
##gff-version 3
chromosome1	RefSeq	gene	1000	5000	.	+	.	ID=gene1;Name=GENE1
chromosome1	RefSeq	mRNA	1000	5000	.	+	.	ID=transcript1;Parent=gene1
chromosome1	RefSeq	exon	1000	1500	.	+	.	ID=exon1;Parent=transcript1
chromosome1	RefSeq	CDS	1200	1400	.	+	0	ID=cds1;Parent=transcript1
```

### FASTQ文件 | FASTQ Files

支持配对的FASTQ文件，文件名需包含read标识符（如R1/R2, _1/_2等）。

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
dual_rnaseq_output/
├── 01.index/                           # 索引文件
│   ├── host.hisat2.index.*.ht2
│   └── pathogen.hisat2.index.*.ht2
├── 02.classification/                  # Reads分类结果
│   └── sample1/
│       ├── sample1.host.bam            # 宿主特异性reads
│       ├── sample1.pathogen.bam        # 病原体特异性reads
│       ├── sample1.ambiguous.bam       # 无法确定的reads
│       └── sample1.unassigned.bam      # 未比对上的reads
├── 03.quantification/                  # 定量结果
│   ├── host/
│   │   ├── sample1.gtf
│   │   └── sample1.fpkm.txt
│   └── pathogen/
│       ├── sample1.gtf
│       └── sample1.fpkm.txt
├── 04.expression_matrix/               # 表达矩阵
│   ├── host_matrix.txt                 # 宿主表达矩阵
│   └── pathogen_matrix.txt             # 病原体表达矩阵
└── dual_rnaseq_summary.txt             # 分析总结报告
```

### 关键输出文件说明 | Key Output Files Description

- **species_matrix.txt**: 包含gene_id, transcript_id, cov, FPKM, TPM, sample列
- **dual_rnaseq_summary.txt**: 完整的分析流程总结报告
- **classification BAM文件**: 分类后的BAM文件，可用于后续可视化分析

## 分析流程 | Analysis Pipeline

1. **索引构建**: 分别为两个物种构建HISAT2索引
2. **样本解析**: 自动识别或从文件读取样本信息
3. **双基因组比对**: 每个样本的reads同时比对到两个参考基因组
4. **物种分类**: 根据MAPQ值将reads分类到不同物种
5. **定量分析**: 分别对两个物种的reads进行StringTie定量
6. **矩阵生成**: 合并生成双物种表达矩阵

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **HISAT2** (版本 2.2.0+)
- **StringTie** (版本 2.0+)
- **Samtools** (版本 1.10+)
- **Python** (版本 3.7+)
- **Python包**:
  - `pysam` - BAM文件处理
  - `pandas` - 数据处理
  - `click` - 命令行界面

### 安装依赖 | Installing Dependencies

```bash
# 安装HISAT2
# Ubuntu/Debian
sudo apt-get install hisat2

# CentOS/RHEL
sudo yum install hisat2

# 安装StringTie
sudo apt-get install stringtie

# 安装Samtools
sudo apt-get install samtools

# 安装Python包
pip install pysam pandas click
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐8核以上）
- **RAM**: 最少16GB（大基因组推荐32GB以上）
- **存储**: 预留基因组文件大小10倍的磁盘空间

## 使用示例 | Usage Examples

### 示例1：病原体-宿主互作分析

```bash
biopytools dual-rnaseq \
    --species1-name mouse \
    --species1-genome mouse.fa \
    --species1-gtf mouse.gtf \
    --species2-name salmonella \
    --species2-genome salmonella.fa \
    --species2-gtf salmonella.gtf \
    -i /data/mouse_salmonella_rnaseq \
    -o infection_results
```

### 示例2：植物-微生物共生分析

```bash
biopytools dual-rnaseq \
    --species1-name plant \
    --species1-genome plant.fa \
    --species1-gtf plant.gtf \
    --species2-name bacteria \
    --species2-genome bacteria.fa \
    --species2-gtf bacteria.gtf \
    -i /data/plant_bacteria_rnaseq \
    -o symbiosis_results \
    -t 24
```

### 示例3：使用样本信息文件

```bash
# 创建样本信息文件 samples.txt
cat > samples.txt << EOF
sample1	/data/sample1.R1.fq.gz	/data/sample1.R2.fq.gz
sample2	/data/sample2.R1.fq.gz	/data/sample2.R2.fq.gz
EOF

biopytools dual-rnaseq \
    --species1-name host \
    --species1-genome host.fa \
    --species1-gtf host.gtf \
    --species2-name pathogen \
    --species2-genome pathogen.fa \
    --species2-gtf pathogen.gtf \
    -i samples.txt \
    -o results
```

## 注意事项 | Important Notes

1. **基因组版本**: 确保基因组FASTA文件和GTF注释文件的版本一致
2. **文件格式**: GTF文件必须包含基因、转录本、外显子等特征信息
3. **内存使用**: 双基因组比对会消耗较多内存，建议预留足够的系统资源
4. **分类准确度**: min_mapq参数影响分类准确度，默认值为20，可根据数据质量调整
5. **样本命名**: FASTQ文件名需要包含read标识符以便识别配对关系

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "hisat2-build: command not found" 错误**
```bash
# 安装HISAT2
sudo apt-get install hisat2  # Ubuntu/Debian
sudo yum install hisat2      # CentOS/RHEL
```

**Q: "pysam: no module named 'pysam'" 错误**
```bash
# 安装pysam
pip install pysam
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools dual-rnaseq ... -t 4
```

**Q: 分类准确度低**
```bash
# 提高min_mapq阈值
biopytools dual-rnaseq ... --min-mapq 30
```

## 结果解读指南 | Result Interpretation Guide

### 分类统计说明

每次样本分类后会输出统计信息：
- **species1**: 分配到物种1的reads数量和比例
- **species2**: 分配到物种2的reads数量和比例
- **ambiguous**: 无法确定物种的reads数量和比例（两个基因组都有高质量比对）
- **unassigned**: 未比对或低质量比对的reads数量和比例

### 表达矩阵格式

表达矩阵包含以下列：
- **gene_id**: 基因ID
- **transcript_id**: 转录本ID
- **cov**: 覆盖度
- **FPKM**: Fragments Per Kilobase Million
- **TPM**: Transcripts Per Million
- **sample**: 样本名称

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用：

```
Dual RNA-seq Analysis Tool in BioPyTools
Version 1.0.0 (2026)
https://github.com/yourrepo/biopytools
```

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件
