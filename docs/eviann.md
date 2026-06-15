# 📝 EviAnn 基因组注释分析模块

**基于证据的真核生物基因组注释工具 | Evidence-based Eukaryotic Genome Annotation Tool**

## 📖 功能概述 | Overview

EviAnn 基因组注释分析模块是一个基于证据的真核生物基因组注释工具，使用RNA-seq数据和/或相关物种的蛋白质比对进行基因预测。EviAnn输出符合NCBI注释规范的GFF3格式，可直接提交到GenBank。支持蛋白编码基因和长链非编码RNA注释，自动识别假基因，可选功能注释。

## ✨ 主要特性 | Key Features

- **🎯 基于证据的注释**: 使用RNA-seq和蛋白质比对证据进行基因预测
- **📊 NCBI规范兼容**: 输出完全符合NCBI注释规范，可直接提交
- **🚀 快速高效**: 哺乳动物基因组注释少于1小时（24核服务器）
- **🧬 多数据类型支持**: 支持RNA-seq、Iso-seq、Nanopore等数据
- **💡 自动假基因识别**: 自动标记加工假基因
- **🔬 功能注释**: 可选UniProt-SwissProt功能注释
- **🧭 lncRNA注释**: 长链非编码RNA注释
- **⚡ 断点续传**: 自动保存进度，支持中断恢复
- **🧬 线粒体支持**: 支持线粒体和叶绿体注释

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用RNA-seq数据进行注释
biopytools eviann \
    -g genome.fa \
    -r rnaseq_list.txt \
    -t 24

# 使用转录本和蛋白质进行注释
biopytools eviann \
    -g genome.fa \
    -e transcripts.fa \
    -p proteins.fa \
    -t 24

# 完整注释（含功能注释）
biopytools eviann \
    -g genome.fa \
    -r rnaseq_list.txt \
    -p related_species_proteins.fa \
    --functional \
    -t 24
```

### 高级用法 | Advanced Usage

```bash
# 自定义参数的注释
biopytools eviann \
    -g genome.fa \
    -r rnaseq_list.txt \
    -p proteins.fa \
    -t 24 \
    -m 50000 \
    -d 2 \
    --lncrna-tpm 2.0 \
    --partial

# 线粒体基因组注释
biopytools eviann \
    -g mitochondria.fa \
    -e transcripts.fa \
    --mito-contigs mito_contigs.txt \
    -t 12
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 基因组FASTA文件 | `-g genome.fa` |

### 数据输入参数 | Data Input Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-r, --rnaseq` | RNA-seq文件列表 | `-r rnaseq.txt` |
| `-e, --transcripts` | 转录本FASTA文件 | `-e transcripts.fa` |
| `-p, --proteins` | 蛋白质FASTA文件 | `-p proteins.fa` |

**注意**: 必须提供 `-r` 或 `-e` 至少一个参数

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-s, --uniprot` | `自动下载` | UniProt-SwissProt FASTA文件 |
| `-t, --threads` | `1` | 线程数 |
| `-m, --max-intron` | `自动计算` | 最大内含子长度 (bp) |
| `-d, --ploidy` | `2` | 基因组倍性 |
| `-c, --cds-gff` | - | 现有CDS的GFF文件 |
| `--lncrna-tpm` | `1.0` | lncRNA最小TPM阈值 |
| `--partial` | `False` | 包含部分CDS |
| `--functional` | `False` | 执行功能注释 |
| `--mito-contigs` | - | 线粒体contig列表文件 |
| `--extra-gff` | - | 额外的GFF特征文件 |
| `--debug` | `False` | 调试模式（保留中间文件） |
| `--verbose` | `False` | 详细输出 |

## 📁 输入文件格式 | Input File Formats

### RNA-seq文件列表格式 | RNA-seq File List Format

`rnaseq.txt` 文件格式：

```
/path/file1_R1.fastq /path/file1_R2.fastq fastq
/path/file2.bam bam
/path/file3.fasta isoseq
/path/file4_R1.fastq /path/file4_R2.fastq /path/file4_longreads.fasta mix
```

**格式说明**:
- 每行代表一个实验的数据
- 支持的数据类型标签：
  - `fastq`: Illumina RNA-seq (fastq格式)
  - `fasta`: Illumina RNA-seq (fasta格式)
  - `bam`: 已比对的Illumina reads
  - `bam_isoseq`: 已比对的PacBio Iso-seq
  - `isoseq`: PacBio Iso-seq (fasta/fastq)
  - `mix`: Illumina + long reads混合
  - `bam_mix`: BAM格式的混合数据

### 基因组文件 | Genome File

标准FASTA格式的基因组序列：

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

## 💡 使用示例 | Usage Examples

### 示例1：基本RNA-seq注释 | Example 1: Basic RNA-seq Annotation

```bash
# 使用Illumina RNA-seq数据进行注释
biopytools eviann \
    -g Arabidopsis_thaliana.fa \
    -r rnaseq_data.txt \
    -t 24
```

### 示例2：结合RNA-seq和蛋白质 | Example 2: Combined RNA-seq and Proteins

```bash
# 使用RNA-seq和相关物种蛋白质
biopytools eviann \
    -g plant_genome.fa \
    -r rnaseq_list.txt \
    -p related_species_proteins.fa \
    -t 24
```

### 示例3：仅使用转录本数据 | Example 3: Transcripts Only

```bash
# 仅使用组装的转录本（无RNA-seq）
biopytools eviann \
    -g genome.fa \
    -e assembled_transcripts.fa \
    -p proteins.fa \
    -t 24
```

### 示例4：完整注释（含功能注释）| Example 4: Full Annotation with Functional

```bash
# 包含功能注释的完整流程
biopytools eviann \
    -g genome.fa \
    -r rnaseq.txt \
    -p proteins.fa \
    --functional \
    -t 24 \
    --verbose
```

### 示例5：线粒体基因组注释 | Example 5: Mitochondrial Genome Annotation

```bash
# 线粒体基因组注释（不同遗传密码）
biopytools eviann \
    -g mitochondria.fa \
    -e transcripts.fa \
    --mito-contigs mito_contigs.txt \
    -t 12
```

### 示例6：自定义参数注释 | Example 6: Custom Parameters

```bash
# 自定义内含子长度和lncRNA阈值
biopytools eviann \
    -g genome.fa \
    -r rnaseq.txt \
    -m 100000 \
    --lncrna-tpm 5.0 \
    --partial \
    -t 24
```

## 📊 输出结果 | Output Results

### 输出文件 | Output Files

EviAnn使用输入基因组文件名作为前缀，例如输入为`genome.fa`，则输出：

```
genome.fa.pseudo_label.gff      # 最终GFF3注释文件
genome.fa.proteins.fasta       # 蛋白质序列
genome.fa.transcripts.fasta    # 转录本序列
```

### GFF3输出格式 | GFF3 Output Format

EviAnn输出的GFF3文件包含以下属性：

**蛋白编码基因属性**:
- `ID`: 转录本ID
- `Parent`: 父特征ID
- `EvidenceProteinID`: 证据蛋白质ID
- `EvidenceTranscriptID`: 证据转录本ID
- `StartCodon`: 起始密码子
- `StopCodon`: 终止密码子
- `Class`: 匹配等级 (=, k, c, etc.)
- `Evidence`: 证据类型 (complete/protein_only/transcript_only)
- `pseudo=true`: 假基因标记（如适用）

**lncRNA属性**:
- `ID`: 转录本ID
- `Parent`: 父特征ID
- `EvidenceTranscriptID`: 证据转录本ID

### GFF3示例 | GFF3 Example

```gff3
chromosome1	EviAnn	gene	29462	43759	.	-	.	ID=XLOC_000048;geneID=XLOC_000048;type=protein_coding
chromosome1	EviAnn	mRNA	32745	43754	.	-	.	ID=XLOC_000048-mRNA-1;Parent=XLOC_000048;EvidenceProteinID=XP_001352289.2;EvidenceTranscriptID=MSTRG_00000148:4:7.70;StartCodon=atg;StopCodon=TGA;Class==;Evidence=complete
chromosome1	EviAnn	exon	32745	33125	.	-	.	Parent=XLOC_000048-mRNA-1
chromosome1	EviAnn	CDS	33018	33125	.	-	0	Parent=XLOC_000048-mRNA-1
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

**外部依赖**（需预安装）:
- **minimap2**: 序列比对
- **HISAT2**: RNA-seq比对

**内置依赖**（EviAnn自带）:
- StringTie 3.0.0
- gffread 0.12.7/0.12.6
- BLAST+ 2.8.1+
- TransDecoder 5.7.1
- samtools 1.15.1
- ufasta 1.0
- miniprot v0.15-r270

### 软件安装 | Software Installation

```bash
# 方法1: Bioconda安装（推荐）
conda create -n eviann_v.2.0.5
conda activate eviann_v.2.0.5
conda install -c bioconda eviann

# 方法2: 从GitHub安装
wget https://github.com/alekseyzimin/EviAnn_release/releases/download/v2.0.6/EviAnn-2.0.6.tar.gz
tar xvzf EviAnn-2.0.6.tar.gz
cd EviAnn-2.0.6
./install.sh
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐24核以上）
- **RAM**: 最少8GB（大基因组推荐64GB以上）
- **存储**: 预留基因组文件大小10倍的磁盘空间
- **网络**: 如需下载UniProt数据库，建议稳定网络

## ⚠️ 注意事项 | Important Notes

1. **数据要求**: 必须提供RNA-seq数据(`-r`)或转录本数据(`-e`)至少一种
2. **基因组大小**: 支持最大32Gbp的基因组
3. **内含子长度**: 默认自动计算，可根据需要手动设置
4. **断点续传**: EviAnn自动保存进度，中断后重新运行相同命令即可恢复
5. **染色体命名**: 确保基因组FASTA和输入文件的染色体名称一致
6. **线粒体注释**: 使用`--mito-contigs`指定线粒体contig以使用正确的遗传密码
7. **功能注释**: 使用`--functional`启用UniProt功能注释（会下载数据库）

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "必须提供RNA-seq或转录本数据" 错误**
```bash
# 解决方案：至少提供-r或-e参数
biopytools eviann -g genome.fa -r rnaseq.txt  # 正确
biopytools eviann -g genome.fa -e transcripts.fa  # 正确
biopytools eviann -g genome.fa  # 错误
```

**Q: "EviAnn未找到" 错误**
```bash
# 检查EviAnn安装
which eviann.sh
# 或使用conda安装
conda install -c bioconda eviann
```

**Q: 内存不足错误**
```bash
# 解决方案：减少线程数或使用更小的数据集
biopytools eviann -g genome.fa -r rnaseq.txt -t 12  # 减少线程
```

**Q: 内含子长度设置不合理**
```bash
# 解决方案：手动设置合理的内含子长度
biopytools eviann -g genome.fa -r rnaseq.txt -m 50000  # 设置为50kb
```

## 📚 相关资源 | Related Resources

- [EviAnn GitHub](https://github.com/alekseyzimin/EviAnn_release)
- [EviAnn预印本](https://www.biorxiv.org/content/10.1101/2025.05.07.652745v2)
- [NCBI注释规范](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/)
- [GFF3格式规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用EviAnn，请引用：

> EviAnn manuscript is under review. The preprint is available at: https://www.biorxiv.org/content/10.1101/2025.05.07.652745v2

**致谢 | Acknowledgments**:

开发得到NSF grant IOS-2432298，NIH grants R01-HG006677和R35-GM130151的支持。
