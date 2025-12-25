# 🧬 BLAST序列比对分析模块 v2.0

**标准化BLAST比对分析工具 | Standardized BLAST Sequence Alignment Analysis Tool**

## 📖 功能概述 | Overview

BLAST序列比对分析模块是一个基于BLAST+套件的标准化序列比对分析工具，严格遵循BioPyTools开发规范。该模块支持多种BLAST算法、批量序列处理、标准化的日志输出和参数命名，适用于基因功能注释、同源性分析和序列比较研究。

## ✨ 主要特性 | Key Features

- **🔧 标准化参数**: 严格遵循BioPyTools参数命名规范
- **⚡ 懒加载设计**: CLI采用懒加载，显著提升`-h`查看帮助文档的速度
- **📊 标准化日志**: 采用标准日志格式，INFO输出到stdout，WARNING/ERROR输出到stderr
- **🔍 多算法支持**: 支持blastn、blastp、blastx、tblastn、tblastx等5种BLAST算法
- **⚡ 高性能处理**: 支持多线程并行处理
- **🛡️ 质量控制**: 完善的参数验证和错误处理
- **🔄 批量处理**: 支持单文件和目录批量处理
- **📈 详细统计**: 自动生成详细的分析统计报告

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本BLAST分析
biopytools blast -i query.fa -r database.fa -o results/

# 批量目录分析
biopytools blast -i sequences_dir/ -r database.fa -o results/ -t 8

# 蛋白质比对
biopytools blast -i proteins.fa -r protein_db.fa -o results/ --blast-type blastp

# 高严格性过滤
biopytools blast -i query.fa -r database.fa -o results/ --min-identity 90 --min-coverage 80
```

### 详细用法 | Advanced Usage

```bash
# 详细输出模式
biopytools blast -i input.fa -r db.fa -o results/ -v

# 静默模式
biopytools blast -i input.fa -r db.fa -o results/ --quiet

# 模拟运行
biopytools blast -i input.fa -r db.fa -o results/ --dry-run

# 强制覆盖
biopytools blast -i input.fa -r db.fa -o results/ -f

# 指定日志文件
biopytools blast -i input.fa -r db.fa -o results/ --log-file blast.log
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 📁 输入文件或目录路径 | `-i query.fa` |
| `-r, --reference` | 🎯 目标基因序列文件 | `-r database.fa` |

### 常用可选参数 | Common Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `./blast_output` | 📂 输出目录 |
| `-p, --prefix` | `blast_output` | 📝 输出文件前缀 |
| `-t, --threads` | `4` | 🧵 线程数 |
| `-q, --quality` | `1e-5` | 📊 E-value阈值 |
| `-m, --memory` | `8G` | 💾 内存限制 |

### 样本信息参数 | Sample Information Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sample-id` | None | 🏷️ 样本ID |
| `--sample-name` | None | 📛 样本名称 |
| `--read-group` | None | 📄 Read Group信息 |

### 质控参数 | Quality Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-quality` | `20` | 🎯 最小质量值 |
| `--min-length` | `50` | 📏 最小序列长度 |
| `--min-depth` | `10` | 📊 最小测序深度 |
| `--max-depth` | `1000` | 📈 最大测序深度 |
| `--mapping-quality` | `20` | 🎯 最小mapping质量 |

### 日志控制参数 | Logging Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-v, --verbose` | `0` | 详细输出级别（-v=INFO, -vv=DEBUG） |
| `--quiet` | `False` | 静默模式（仅ERROR） |
| `--log-file` | None | 📄 日志文件路径 |

### 执行控制参数 | Execution Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-f, --force` | `False` | 强制覆盖已存在文件 |
| `--dry-run` | `False` | 模拟运行不执行 |
| `--keep-intermediate` | `False` | 保留中间文件 |

### BLAST特定参数 | BLAST-Specific Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--blast-type` | `blastn` | ⚡ BLAST程序类型 |
| `--max-target-seqs` | `10` | 🔢 最大目标序列数 |
| `--min-identity` | `70.0` | 🎯 最小序列相似度 (%) |
| `--min-coverage` | `50.0` | 📐 最小覆盖度 (%) |
| `--target-db-type` | `nucl` | 💾 目标数据库类型 |
| `--high-quality-evalue` | `1e-10` | ⭐ 高质量比对E-value阈值 |

## 🧬 BLAST算法选择指南 | BLAST Algorithm Selection Guide

### blastn (核酸-核酸)
- **查询序列**: DNA/RNA
- **数据库**: DNA/RNA
- **用途**: 基因定位、转录本比对、SNP分析
- **示例**: `biopytools blast -i genome.fa -r ref_genome.fa --blast-type blastn`

### blastp (蛋白质-蛋白质)
- **查询序列**: 蛋白质
- **数据库**: 蛋白质
- **用途**: 蛋白质功能预测、结构域分析
- **示例**: `biopytools blast -i proteins.fa -r protein_db.fa --blast-type blastp --target-db-type prot`

### blastx (核酸-蛋白质)
- **查询序列**: DNA/RNA (6框翻译)
- **数据库**: 蛋白质
- **用途**: 编码基因功能注释、ORF预测
- **示例**: `biopytools blast -i transcripts.fa -r protein_db.fa --blast-type blastx --target-db-type prot`

### tblastn (蛋白质-核酸)
- **查询序列**: 蛋白质
- **数据库**: DNA/RNA (6框翻译)
- **用途**: 基因发现、假基因识别
- **示例**: `biopytools blast -i proteins.fa -r genome.fa --blast-type tblastn`

### tblastx (核酸-核酸翻译)
- **查询序列**: DNA/RNA (6框翻译)
- **数据库**: DNA/RNA (6框翻译)
- **用途**: 进化分析、假基因比较
- **示例**: `biopytools blast -i sequences.fa -r genome.fa --blast-type tblastx`

## 📊 标准化日志输出 | Standardized Logging Output

### stdout输出内容 | stdout Content

```
[2025-12-19 10:30:15] INFO: Pipeline started
[2025-12-19 10:30:15] INFO: Program: BLAST Analysis
[2025-12-19 10:30:15] INFO: Version: 2.0.0
[2025-12-19 10:30:15] INFO: ============================================================
[2025-12-19 10:30:15] INFO: STEP 1: Creating BLAST Database
[2025-12-19 10:30:15] INFO: ============================================================
[2025-12-19 10:30:20] INFO: Database created successfully: /path/to/database.db
[2025-12-19 10:30:20] INFO: ✅ Creating BLAST Database completed successfully
[2025-12-19 10:30:20] INFO: ============================================================
[2025-12-19 10:30:20] INFO: STEP 2: Running BLAST Alignment
[2025-12-19 10:30:20] INFO: ============================================================
[2025-12-19 10:30:20] INFO: Processing sample 1/3: sample1
[2025-12-19 10:30:25] INFO: Sample sample1 alignment completed
[2025-12-19 10:30:30] INFO: Processing sample 2/3: sample2
[2025-12-19 10:30:35] INFO: Sample sample2 alignment completed
[2025-12-19 10:30:40] INFO: Processing sample 3/3: sample3
[2025-12-19 10:30:45] INFO: Sample sample3 alignment completed
[2025-12-19 10:30:45] INFO: BLAST alignment completed. Results: /path/to/results.tsv
[2025-12-19 10:30:45] INFO: ✅ Running BLAST Alignment completed successfully
[2025-12-19 10:30:45] INFO: ============================================================
[2025-12-19 10:30:45] INFO: STEP 3: Processing Results
[2025-12-19 10:30:45] INFO: ============================================================
[2025-12-19 10:30:45] INFO: Statistics:
[2025-12-19 10:30:45] INFO:   Total alignments: 1,234
[2025-12-19 10:30:45] INFO:   Samples count: 3
[2025-12-19 10:30:45] INFO:   Unique queries: 456
[2025-12-19 10:30:45] INFO:   Unique subjects: 789
[2025-12-19 10:30:45] INFO: ✅ Processing Results completed successfully
[2025-12-19 10:30:45] INFO: ============================================================
[2025-12-19 10:30:45] INFO: Pipeline Summary
[2025-12-19 10:30:45] INFO: ============================================================
[2025-12-19 10:30:45] INFO: Total runtime: 30.45 seconds
[2025-12-19 10:30:45] INFO: Pipeline completed successfully
[2025-12-19 10:30:45] INFO: Sample count: 3
[2025-12-19 10:30:45] INFO: BLAST type: blastn
[2025-12-19 10:30:45] INFO: Target database: /path/to/reference.fa
[2025-12-19 10:30:45] INFO: Output directory: /path/to/output
[2025-12-19 10:30:45] INFO: ✅ Results saved to: /path/to/output
```

### stderr输出内容 | stderr Content

```
[2025-12-19 10:30:25] WARNING: Low quality alignments detected in sample1
[2025-12-19 10:30:35] WARNING: High duplication rate detected in sample2
[2025-12-19 10:30:45] ERROR: Failed to process sample3: insufficient memory
```

## 💡 使用示例 | Usage Examples

### 示例1：基本单文件分析 | Example 1: Basic Single File Analysis

```bash
# 单个序列文件比对
biopytools blast -i query_sequence.fa -r target_database.fa -o single_analysis

# 输出：
# [INFO] Pipeline started
# [INFO] Processing 1 sample...
# [INFO] Results: single_analysis/blast_summary_results.tsv
```

### 示例2：批量目录分析 | Example 2: Batch Directory Analysis

```bash
# 批量处理多个FASTA文件
biopytools blast \
    -i /path/to/sequences/ \
    -r /path/to/database.fa \
    -o /path/to/results/ \
    -t 16 \
    -v

# 自动扫描目录中的*.fa文件
# 支持verbose输出显示详细信息
```

### 示例3：蛋白质比对分析 | Example 3: Protein Alignment Analysis

```bash
# 蛋白质序列比对
biopytools blast \
    -i proteins.fa \
    -r protein_database.fa \
    -o protein_analysis \
    --blast-type blastp \
    --target-db-type prot \
    --min-identity 85
```

### 示例4：跨物种同源搜索 | Example 4: Cross-Species Homology Search

```bash
# 跨物种同源基因搜索
biopytools blast \
    -i query_genes.fa \
    -r target_proteome.fa \
    -o homology_search \
    --blast-type blastx \
    --min-identity 60 \
    --evalue 1e-6
```

### 示例5：质量控制分析 | Example 5: Quality Control Analysis

```bash
# 高质量比对分析
biopytools blast \
    -i high_quality_sequences.fa \
    -r curated_database.fa \
    -o quality_analysis \
    --min-identity 95 \
    --min-coverage 90 \
    --high-quality-evalue 1e-15
```

### 示例6：调试和测试 | Example 6: Debugging and Testing

```bash
# 模拟运行（不实际执行）
biopytools blast \
    -i test.fa \
    -r test_db.fa \
    -o test_output \
    --dry-run

# 详细日志记录到文件
biopytools blast \
    -i test.fa \
    -r test_db.fa \
    -o test_output \
    -vv \
    --log-file debug.log
```

## 📁 输出文件说明 | Output Files Description

### 主要输出文件 | Main Output Files

```
output_directory/
├── blast_summary_results.tsv          # 汇总比对结果
│   Sample    qseqid    sseqid    pident    ...    evalue    bitscore
│   sample1   gene1    target1   95.2     ...    1e-10      156.8
│   sample2   gene2    target2   89.5     ...    5e-08      134.2
├── sample1_blastn_results.tsv        # 单个样品的详细结果
├── sample2_blastn_results.tsv        # 单个样品的详细结果
└── sample3_blastn_results.tsv        # 单个样品的详细结果
```

### 结果文件格式 | Result File Format

**汇总文件格式**：
```
Sample    qseqid    sseqid    pident    length    mismatch    gapopen    qstart    qend    sstart    send    evalue    bitscore
```

**字段说明**：
- `Sample`: 样品名称
- `qseqid`: 查询序列ID
- `sseqid`: 目标序列ID
- `pident`: 相似度百分比
- `length`: 比对长度
- `mismatch`: 错配数
- `gapopen`: Gap开放数
- `qstart/qend`: 查询序列起止位置
- `sstart/send`: 目标序列起止位置
- `evalue`: E-value
- `bitscore`: Bit score

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

**必需工具**：
- **BLAST+** (版本 2.12+)
  - `makeblastdb`: 数据库构建工具
  - `blastn`: 核酸-核酸比对
  - `blastp`: 蛋白质-蛋白质比对
  - `blastx`: 核酸-蛋白质比对
  - `tblastn`: 蛋白质-核酸比对
  - `tblastx`: 核酸-核酸翻译比对

**Python包**：
```bash
pip install click
```

### 安装BLAST+ | Installing BLAST+

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install ncbi-blast+

# CentOS/RHEL
sudo yum install ncbi-blast+

# Conda
conda install -c bioconda blast
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（大数据集推荐16GB以上）
- **存储**: 至少预留数据库大小2倍的磁盘空间
- **临时存储**: 额外10GB用于临时文件

## ⚠️ 注意事项 | Important Notes

1. **数据质量**: 确保输入FASTA文件格式正确
2. **参数优化**: 根据数据类型选择合适的BLAST算法和参数
3. **内存管理**: 大数据集分析时注意内存使用情况
4. **权限设置**: 确保对输入文件和输出目录有读写权限
5. **数据库准备**: 确保目标数据库已正确构建

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "BLAST command not found" 错误**
```bash
# 检查BLAST安装
which blastn
blastn -version

# 重新安装或添加到PATH
export PATH=$PATH:/path/to/blast/bin
```

**Q: "makeblastdb failed" 错误**
```bash
# 检查数据库文件
ls -la your_database.fa
head -5 your_database.fa

# 检查文件格式（必须是FASTA格式）
grep "^>" your_database.fa | head -5
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools blast -i input.fa -r db.fa -o results -t 2

# 或增加系统内存
```

**Q: 无比对结果**
```bash
# 降低E-value阈值
biopytools blast -i input.fa -r db.fa -o results -q 1e-3

# 降低相似度要求
biopytools blast -i input.fa -r db.fa -o results --min-identity 50
```

### 性能优化建议 | Performance Optimization

1. **线程配置**:
   - 小数据集：4-8线程
   - 中等数据集：16-32线程
   - 大数据集：32-64线程

2. **数据库优化**:
   - 使用压缩格式节省存储空间
   - 定期清理临时文件
   - 考虑数据库分割策略

3. **I/O优化**:
   - 使用SSD存储
   - 将临时文件放在高速存储
   - 避免网络存储瓶颈

## 📚 相关资源 | Related Resources

### 学术文献 | Academic Papers

- [Altschul SF et al. (1990) Basic local alignment search tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Altschul SF et al. (1997) Gapped BLAST and PSI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [NCBI BLAST+ User Manual](https://ftp.ncbi.nlm.nih.gov/blast/documents/blast+manual.pdf)

### 相关工具 | Related Tools

- [DIAMOND](https://github.com/bbuchfink/diamond) - 快速蛋白质序列比对
- [MMseqs2](https://mmseqs.org/) - 超快速序列搜索和聚类
- [HMMER](http://hmmer.org/) - 序列分析工具包

### 教程和文档 | Tutorials and Documentation

- [NCBI BLAST Tutorial](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/)
- [BLAST+ Command Line Applications User Manual](https://ftp.ncbi.nlm.nih.gov/blast/documents/blast+manual.pdf)
- [BioPython BLAST Module](https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec-blast)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关方法学文献：

```
Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. (1990)
Basic local alignment search tool.
J Mol Biol Biol 215:403-410.

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. (1997)
Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.
Nucleic Acids Res 25:3389-3402.
```