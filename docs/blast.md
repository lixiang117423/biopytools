# 🔍 BLAST比对分析模块

**高效的序列比对和相似性搜索工具 | Efficient Sequence Alignment and Similarity Search Tool**

## 📖 功能概述 | Overview

BLAST比对分析模块是一个强大的序列比对工具，支持多种BLAST算法，提供灵活的输入方式和丰富的参数配置，适用于各种生物信息学序列分析需求。

## ✨ 主要特性 | Key Features

- **🔬 多种BLAST算法**: 支持 blastn, blastp, blastx, tblastn, tblastx
- **📁 灵活输入**: 支持单文件、目录批处理、样品映射文件
- **⚡ 高性能**: 多线程并行处理，最大化计算效率
- **🎯 智能过滤**: 基于E-value、相似度、覆盖度的质量控制
- **🤖 自动化**: 自动样品检测和数据库构建
- **📊 结果分析**: 详细的比对结果统计和质量评估

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本BLAST比对
biopytools blast -i input.fa -t target.fa -o blast_results

# 批量处理多个文件
biopytools blast -i /path/to/sequences/ -t target.fa -o batch_results

# 使用样品映射文件
biopytools blast -s sample_map.txt -t target.fa -o mapped_results
```

### 高级用法 | Advanced Usage

```bash
# 蛋白质序列比对
biopytools blast -i proteins.fa -t target_proteins.fa \
    --blast-type blastp --target-db-type prot \
    -e 1e-10 --min-identity 80 -o protein_results

# 高通量分析
biopytools blast -i sequences/ -t targets.fa \
    --threads 32 --max-target-seqs 20 \
    --min-coverage 70 --high-quality-evalue 1e-15 \
    -o htp_results
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-t, --target-file` | 目标基因序列文件 | `-t reference.fa` |

### 输入选项 | Input Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | - | 输入文件或目录路径 |
| `-s, --sample-map-file` | - | 样品映射文件（格式：文件路径<TAB>样品名称） |
| `--input-suffix` | `*.fa` | 输入文件后缀模式（目录输入时使用） |
| `--auto-detect-samples` | `True` | 自动检测样品名称 |
| `--sample-name-pattern` | `([^/]+?)(?:\.fa\|\.fasta\|\.fna)?$` | 样品名称提取正则表达式 |

### BLAST配置 | BLAST Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--blast-type` | `blastn` | BLAST程序类型 (blastn/blastp/blastx/tblastn/tblastx) |
| `--target-db-type` | `nucl` | 目标数据库类型 (nucl/prot) |
| `-e, --evalue` | `1e-5` | E-value阈值 |
| `--max-target-seqs` | `10` | 最大目标序列数 |
| `--word-size` | `11` | 词大小（适用于blastn/tblastx） |

### 质量控制 | Quality Control

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-identity` | `70.0` | 最小序列相似度 (%) |
| `--min-coverage` | `50.0` | 最小覆盖度 (%) |
| `--high-quality-evalue` | `1e-10` | 高质量比对E-value阈值 |

### 性能配置 | Performance Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-j, --threads` | `88` | 线程数 |

### 工具路径 | Tool Paths

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--makeblastdb-path` | `makeblastdb` | makeblastdb程序路径 |
| `--blastn-path` | `blastn` | blastn程序路径 |
| `--blastp-path` | `blastp` | blastp程序路径 |
| `--blastx-path` | `blastx` | blastx程序路径 |
| `--tblastn-path` | `tblastn` | tblastn程序路径 |
| `--tblastx-path` | `tblastx` | tblastx程序路径 |

## 📁 输入文件格式 | Input File Formats

### 序列文件 | Sequence Files

支持标准FASTA格式文件：

```fasta
>sequence_1
ATCGATCGATCGATCGATCG
>sequence_2
GCTAGCTAGCTAGCTAGCTA
```

**支持的文件扩展名**: `.fa`, `.fasta`, `.fna`, `.fas`

### 样品映射文件 | Sample Mapping File

制表符分隔的文本文件，格式：

```
/path/to/sample1.fa    Sample_001
/path/to/sample2.fa    Sample_002
/path/to/sample3.fa    Sample_003
```

## 💡 使用示例 | Usage Examples

### 示例1：基本核酸序列比对 | Example 1: Basic Nucleotide Alignment

```bash
# 将查询序列与参考基因组比对
biopytools blast \
    -i query_sequences.fa \
    -t reference_genome.fa \
    -o nucleotide_blast_results \
    --blast-type blastn \
    -e 1e-6 \
    --min-identity 85 \
    --threads 16
```

### 示例2：蛋白质序列比对 | Example 2: Protein Sequence Alignment

```bash
# 蛋白质序列数据库搜索
biopytools blast \
    -i protein_queries.fa \
    -t protein_database.fa \
    -o protein_blast_results \
    --blast-type blastp \
    --target-db-type prot \
    -e 1e-8 \
    --min-identity 50 \
    --min-coverage 60 \
    --max-target-seqs 20
```

### 示例3：批量处理 | Example 3: Batch Processing

```bash
# 处理目录中的所有FASTA文件
biopytools blast \
    -i /data/sequences/ \
    -t targets.fa \
    -o batch_results \
    --input-suffix "*.fasta" \
    --threads 32 \
    --auto-detect-samples
```

### 示例4：使用样品映射 | Example 4: Using Sample Mapping

```bash
# 创建样品映射文件
echo -e "/data/sample1.fa\tPatient_A" > samples.map
echo -e "/data/sample2.fa\tPatient_B" >> samples.map

# 运行分析
biopytools blast \
    -s samples.map \
    -t disease_genes.fa \
    -o patient_analysis \
    --blast-type blastn \
    -e 1e-10 \
    --high-quality-evalue 1e-15
```

### 示例5：跨物种比对 | Example 5: Cross-species Alignment

```bash
# 翻译搜索（核酸vs蛋白质）
biopytools blast \
    -i genomic_sequences.fa \
    -t protein_families.fa \
    -o cross_species_results \
    --blast-type blastx \
    --target-db-type prot \
    -e 1e-5 \
    --min-identity 40 \
    --min-coverage 30
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
blast_results/
├── blast_database/          # BLAST数据库文件
│   ├── target.fa
│   ├── target.fa.nhr
│   ├── target.fa.nin
│   └── target.fa.nsq
├── blast_results/           # 原始BLAST结果
│   ├── sample1_blast.txt
│   ├── sample2_blast.txt
│   └── ...
├── summary/                 # 汇总分析结果
│   ├── blast_summary.txt
│   ├── quality_statistics.txt
│   └── best_hits_summary.txt
└── logs/                   # 运行日志
    ├── blast_analysis.log
    └── sample_processing.log
```

### 结果文件说明 | Result Files Description

- **blast_results/**: 每个样品的详细BLAST比对结果
- **summary/blast_summary.txt**: 所有样品的比对结果汇总
- **summary/quality_statistics.txt**: 比对质量统计信息
- **summary/best_hits_summary.txt**: 最佳匹配结果汇总
- **logs/**: 详细的分析日志和错误信息

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **NCBI BLAST+** (版本 2.10.0+)
  - `makeblastdb`
  - `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`
- **Python** (版本 3.7+)
- **Click** (命令行界面)

### 安装BLAST+ | Installing BLAST+

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# CentOS/RHEL
sudo yum install ncbi-blast+

# macOS (使用Homebrew)
brew install blast

# Conda
conda install -c bioconda blast
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐16核以上）
- **RAM**: 最少8GB（大数据集推荐32GB以上）
- **存储**: SSD硬盘（提升I/O性能）
- **网络**: 如需下载数据库，建议高速网络连接

## ⚠️ 注意事项 | Important Notes

1. **内存使用**: 大型数据库和查询文件可能需要大量内存
2. **磁盘空间**: 确保有足够空间存储中间文件和结果
3. **参数调优**: 根据具体应用场景调整E-value和相似度阈值
4. **数据质量**: 输入序列质量直接影响比对结果准确性
5. **版本兼容**: 确保BLAST+版本与工具兼容

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "makeblastdb not found" 错误**
```bash
# 检查BLAST+安装
which makeblastdb
# 如未安装，请按上述方法安装BLAST+
```

**Q: 内存不足错误**
```bash
# 减少线程数或处理较小的数据集
biopytools blast -i input.fa -t target.fa --threads 4
```

**Q: 没有找到匹配结果**
```bash
# 放宽筛选条件
biopytools blast -i input.fa -t target.fa \
    -e 1e-3 --min-identity 50 --min-coverage 30
```

**Q: 处理速度过慢**
```bash
# 增加线程数，减少最大目标序列数
biopytools blast -i input.fa -t target.fa \
    --threads 32 --max-target-seqs 5
```

## 📚 相关资源 | Related Resources

- [NCBI BLAST+ 用户手册](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [BLAST算法原理](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)
- [序列比对最佳实践](https://www.ncbi.nlm.nih.gov/books/NBK279684/)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件