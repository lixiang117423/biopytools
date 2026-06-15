# Kmer丰度转VCF工具 | Kmer Abundance to VCF Converter

**高效的kmer丰度矩阵转VCF格式工具 | Efficient Kmer Abundance Matrix to VCF Format Converter**

## 功能概述 | Overview

Kmer丰度转VCF工具是一个专门用于将kmer丰度矩阵文件转换为标准VCF格式的命令行工具。该工具基于聚合频次过滤策略，能够有效去除低频噪音和多拷贝重复序列，适用于kmer GWAS分析、presence/absence变异（PAV）分析等群体遗传学研究。

## 主要特性 | Key Features

- **智能过滤**: 基于聚合频次的两级过滤策略，去除低频噪音和高频重复
- **流式处理**: 采用两遍扫描策略，内存友好，支持超大规模文件（亿级行数）
- **标准输出**: 生成符合VCF v4.2规范的标准格式文件
- **进度监控**: 实时显示处理进度，每100万行提示一次
- **灵活配置**: 支持自定义聚合频次阈值和输出路径

## 核心算法 | Core Algorithm

### 转换规则

1. **丰度转换**: 丰度值 > 0 转换为二进制1（存在），丰度值 = 0 转换为二进制0（缺失）
2. **聚合频次**: 对每个唯一kmer，统计其在所有样本中出现的总次数
3. **过滤范围**: 保留 `MIN_AGG_COUNT <= 聚合频次 <= 样本数` 的kmer
   - 过低：可能是测序错误或稀有kmer
   - 过高：可能是多拷贝重复序列
4. **VCF格式**:
   - CHROM: `Chr_kmer`（统一）
   - POS: 过滤后的行号（1, 2, 3...）
   - REF: `<kmer序列>`
   - ALT: `.`（表示缺失）
   - GT: `0/0`（有kmer）或 `1/1`（无kmer）

### 两遍扫描策略

- **Pass 1**: 统计每个唯一kmer的聚合频次（所有样本出现次数之和）
- **Pass 2**: 根据聚合频次过滤并转换为VCF格式输出

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 转换kmer丰度矩阵为VCF
biopytools kmer2vcf -i abundance_matrix.tsv -o output_prefix

# 指定最小聚合频次
biopytools kmer2vcf -i abundance_matrix.tsv -o output_prefix -m 5

# 指定输出目录
biopytools kmer2vcf -i abundance_matrix.tsv -o output_prefix -O ./results
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-matrix` | 输入kmer丰度矩阵文件（TSV格式）| `-i abundance_matrix.tsv` |
| `-o, --output-prefix` | 输出文件前缀 | `-o output_prefix` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --min-agg-count` | `3` | 最小聚合频次阈值（低于此值的kmer会被过滤）|
| `-O, --output-dir` | `./kmer2vcf_output` | 输出目录路径 |

## 输入文件格式 | Input File Format

### TSV丰度矩阵格式

输入文件为Tab分隔的文本文件（TSV格式）：

```tsv
kmer	test2	test3	test1	W681A	W680A
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGG	2	2	2	230	147
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGGG	2	2	2	255	183
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGGGC	2	2	2	255	205
```

**格式要求**:
- 第一行：表头，第一列为`kmer`，后续列为样本名
- 第一列：kmer序列（字符类型）
- 后续列：各样本的kmer丰度值（整数）

## 输出文件 | Output Files

### 1. VCF文件 (*.vcf)

标准VCF v4.2格式文件，包含变异信息：

```vcf
##fileformat=VCFv4.2
##source=kmer2vcf
##filter_agg_range=3_to_5
##INFO=<ID=KMER,Number=1,Type=String,Description="Kmer sequence">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=Chr_kmer,length=2147483647>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	test2	test3	test1	W681A	W680A
Chr_kmer	1	.	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGG	.	.	PASS	KMER=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGG	GT	0/0	0/0	0/0	0/0	0/0
Chr_kmer	2	.	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGGG	.	.	PASS	KMER=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGGG	GT	0/0	0/0	0/0	0/0	0/0
```

**基因型说明**:
- `0/0`: 纯合参考型，表示该样本有此kmer
- `1/1`: 纯合变异型，表示该样本缺失此kmer

### 2. TXT文件 (*.txt)

过滤后的kmer丰度矩阵（TSV格式），保留原始丰度值：

```tsv
KmerID	test2	test3	test1	W681A	W680A
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGG	2	2	2	230	147
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCAGGTTAGCTCGCCTGGG	2	2	2	255	183
```

### 3. 日志文件 (kmer2vcf.log)

包含详细的处理过程和统计信息：

```
2026-01-26 14:30:00.123 - INFO - ======================================
2026-01-26 14:30:00.124 - INFO -  Kmer丰度转VCF工具|Kmer Abundance to VCF Converter
2026-01-26 14:30:00.125 - INFO - ======================================
2026-01-26 14:30:00.126 - INFO - 输入文件|Input file: /path/to/abundance_matrix.tsv
2026-01-26 14:30:00.127 - INFO - 样本数|Number of samples: 5
...
```

## 使用示例 | Usage Examples

### 示例1：基本转换

```bash
# 转换kmer丰度矩阵，使用默认参数
biopytools kmer2vcf \
    -i abundance_matrix.tsv \
    -o my_analysis
```

**输出文件**:
- `./kmer2vcf_output/my_analysis.vcf`
- `./kmer2vcf_output/my_analysis.txt`
- `./kmer2vcf_output/kmer2vcf.log`

### 示例2：严格过滤

```bash
# 使用更高的聚合频次阈值，去除更多低频kmer
biopytools kmer2vcf \
    -i abundance_matrix.tsv \
    -o strict_output \
    -m 10 \
    -O ./strict_results
```

### 示例3：处理大规模数据

```bash
# 处理超大规模文件（9.58亿行）
biopytools kmer2vcf \
    -i huge_matrix.tsv \
    -o huge_output \
    -m 3
```

**处理时间参考**:
- 1亿行: ~5-10分钟
- 10亿行: ~50-100分钟
- 具体时间取决于CPU和磁盘IO性能

## 性能优化建议 | Performance Optimization Tips

### 1. 内存优化

工具采用流式处理，内存占用主要在Pass 1：
- 内存需求 ≈ 唯一kmer数 × 100字节
- 对于10万个唯一kmer，内存占用约10MB
- 对于1亿个唯一kmer，内存占用约10GB

### 2. 磁盘IO优化

- 使用SSD存储可显著提升速度（2-5倍）
- 确保输出目录有足够空间（约为输入文件的2-3倍）

### 3. 并行处理

当前版本为单线程处理，适合大多数场景。如需更高性能，建议：
- 将输入文件按染色体或样本拆分
- 并行运行多个实例
- 最后合并VCF文件

## 常见问题 | FAQ

### Q1: 如何选择最小聚合频次阈值？

**A**: 根据样本数量和数据质量选择：
- 小样本量（<10）: 使用默认值3
- 中等样本量（10-50）: 使用5-10
- 大样本量（>50）: 使用10-20
- 严格过滤: 设置更高值（如样本数的20%）

### Q2: 为什么有些kmer被过滤？

**A**: 两种情况：
1. **聚合频次过低** (< MIN_AGG_COUNT): 可能是测序错误或稀有kmer
2. **聚合频次过高** (> 样本数): 可能是多拷贝重复序列

### Q3: 输出的VCF文件如何用于下游分析？

**A**: 可用于：
- **kmer GWAS**: 使用PLINK、GEMMA等工具进行关联分析
- **PAV分析**: 研究基因presence/absence变异
- **群体结构**: 使用VCFtools、ADMIXTURE等分析群体遗传结构
- **系统发育**: 构建物种或品种进化树

### Q4: 如何处理超大文件（>100GB）？

**A**: 建议：
1. 检查可用内存（建议>32GB）
2. 使用SSD存储
3. 考虑先按kmer ID排序和去重
4. 分批处理，最后合并结果

### Q5: 基因型0/0和1/1是什么意思？

**A**:
- **0/0**: 纯合参考型，表示该样本**有**这个kmer
- **1/1**: 纯合变异型，表示该样本**无**这个kmer（缺失）
- 这是presence/absence变异的标准表示方法

## 与其他工具的对比 | Comparison with Other Tools

| 特性 | kmer2vcf | 自定义脚本 | kmtricks |
|------|----------|-----------|----------|
| 易用性 | 高（一键转换）| 低（需要编程）| 中（需要多个步骤）|
| 标准VCF输出 | 是 | 否 | 否 |
| 智能过滤 | 是 | 需自己实现 | 基础过滤 |
| 大规模支持 | 是（流式处理）| 取决于实现 | 是 |
| 文档完善 | 是 | 无 | 中等 |
| 可扩展性 | 高 | 低 | 中 |

## 技术细节 | Technical Details

### 聚合频次计算

对于每个唯一kmer，其聚合频次计算公式：

```
聚合频次(kmer) = Σ(样本i中该kmer的存在标记)
其中存在标记 = 1 if 丰度 > 0 else 0
```

### 过滤逻辑

```
保留条件:
MIN_AGG_COUNT <= 聚合频次(kmer) <= 样本数

过滤原因:
- 聚合频次 < MIN_AGG_COUNT: 低频噪音/测序错误
- 聚合频次 > 样本数: 多拷贝重复序列
```

### VCF编码规则

遵循VCF v4.2规范：
- REF: kmer序列（表示参考状态）
- ALT: `.`（表示deletion，即缺失该kmer）
- GT: `0/0`（有REF）或 `1/1`（有ALT/deletion）

## 引用信息 | Citation

如果在研究中使用此工具，建议引用：

```
Kmer2VCF: A tool for converting kmer abundance matrices to standard VCF format
for presence/absence variation (PAV) analysis.
```

## 许可证 | License

MIT License

## 更新日志 | Changelog

### Version 1.0.0 (2026-01-26)
- 初始版本发布
- 支持kmer丰度矩阵转VCF
- 实现聚合频次过滤策略
- 支持大规模文件流式处理
- 提供详细的日志和进度显示

## 相关资源 | Related Resources

- [VCF格式规范](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [kmer GWAS方法学](https://www.nature.com/articles/s41467-018-05198-5)
- [Presence/Absence变异分析](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3130188/)
- [kmtricks软件](https://github.com/tlebamer/kmtricks)
