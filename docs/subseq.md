# 🧬 序列子集提取模块

**专业的FASTA序列提取工具 | Professional FASTA Sequence Extraction Tool**

## 📖 功能概述 | Overview

序列子集提取模块是一个强大的FASTA序列提取工具，支持多种提取方式，包括基于ID列表、模式匹配和长度筛选。提供高效的序列处理能力、灵活的匹配选项和详细的统计报告，适用于各种生物序列处理和分析任务。

## ✨ 主要特性 | Key Features

- **📋 ID列表提取**: 根据提供的ID列表精确提取指定序列
- **🔍 模式匹配提取**: 支持包含、开头、结尾、正则表达式等多种模式匹配
- **📏 长度筛选提取**: 基于序列长度范围进行灵活筛选
- **🔄 顺序控制**: 支持保持ID列表顺序或原始FASTA顺序
- **⚙️ 灵活配置**: 可配置大小写敏感性、匹配模式等参数
- **📊 详细统计**: 提供完整的提取统计信息和成功率报告
- **🛡️ 错误处理**: 完善的错误处理和用户友好的提示信息
- **📝 双语日志**: 支持中英双语日志记录

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 根据ID列表提取序列
biopytools subseq -i input.fasta -l id_list.txt -o output.fasta

# 模式匹配提取（包含特定字符串）
biopytools subseq -i input.fasta -p "chr1" -o chr1_sequences.fasta

# 长度筛选提取（提取1000-5000bp的序列）
biopytools subseq -i input.fasta -o medium_sequences.fasta \
    --length-only --min-length 1000 --max-length 5000
```

### 高级用法 | Advanced Usage

```bash
# 正则表达式匹配
biopytools subseq -i input.fasta -p "^gene_\d+.*" \
    -o gene_sequences.fasta --pattern-type regex

# 忽略大小写匹配
biopytools subseq -i input.fasta -p "GENE" -o gene_case_insensitive.fasta \
    --ignore-case

# 以特定前缀开头的序列
biopytools subseq -i input.fasta -p "chr1_" -o chr1_sequences.fasta \
    --pattern-type startswith

# 以特定后缀结尾的序列
biopytools subseq -i input.fasta -p "_gene" -o gene_end_sequences.fasta \
    --pattern-type endswith
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTA文件路径 | `-i input.fasta` |
| `-o, --output` | 输出FASTA文件路径 | `-o output.fasta` |

### 提取方式选项（互斥）| Extraction Method Options (Mutually Exclusive)

| 参数 | 描述 | 示例 |
|------|------|------|
| `-l, --id-list` | ID列表文件路径 | `-l sequences.txt` |
| `-p, --pattern` | 模式匹配字符串 | `-p "chr1_"` |
| `--length-only` | 仅使用长度筛选模式 | `--length-only` |

### 模式匹配选项 | Pattern Matching Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--pattern-type` | `contains` | 模式类型：`contains`, `startswith`, `endswith`, `regex` |
| `--ignore-case` | `False` | 忽略大小写匹配 |

### 长度筛选选项 | Length Filtering Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-length` | `0` | 最小序列长度（bp） |
| `--max-length` | `无限制` | 最大序列长度（bp） |

### 其他选项 | Other Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-order` | `False` | 不保持ID列表顺序，使用FASTA原始顺序 |
| `--log-dir` | `当前目录` | 日志文件输出目录 |

## 📁 输入文件格式 | Input File Formats

### FASTA序列文件 | FASTA Sequence File

标准FASTA格式的序列文件：

```fasta
>seq1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr1_seq3
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>gene1
TTTTAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCC
>seq4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

### ID列表文件 | ID List File

每行一个序列ID的文本文件：

```txt
seq1
chr1_seq3
gene1
nonexistent_id
```

**文件要求**:
- 每行一个序列ID
- 空行会被自动忽略
- ID必须与FASTA文件中的序列ID完全匹配

## 💡 使用示例 | Usage Examples

### 示例1：基因组序列提取 | Example 1: Genome Sequence Extraction

```bash
# 提取特定染色体的序列
echo "chr1 chr2 chr3 chr4 chr5" > chromosome_ids.txt

biopytools subseq \
    -i human_genome.fasta \
    -l chromosome_ids.txt \
    -o selected_chromosomes.fasta
```

### 示例2：基因家族序列提取 | Example 2: Gene Family Sequence Extraction

```bash
# 提取所有基因家族成员
biopytools subseq \
    -i annotated_genome.fasta \
    -p "ABC_" \
    -o ABC_gene_family.fasta \
    --pattern-type startswith

# 忽略大小写匹配
biopytools subseq \
    -i all_genes.fasta \
    -p "ABC" \
    -o ABC_family_all.fasta \
    --ignore-case
```

### 示例3：正则表达式复杂匹配 | Example 3: Complex Regex Matching

```bash
# 提取特定命名规则的序列
biopytools subseq \
    -i transcriptome.fasta \
    -p "^gene_[A-Z]\d{3}.*$" \
    -o pattern_genes.fasta \
    --pattern-type regex

# 提取包含特定数字编号的序列
biopytools subseq \
    -i contigs.fasta \
    -p "contig_\d{3,}" \
    -o long_contigs.fasta \
    --pattern-type regex
```

### 示例4：长度筛选应用 | Example 4: Length Filtering Applications

```bash
# 提取完整基因序列（>500bp）
biopytools subseq \
    -i predicted_genes.fasta \
    -o complete_genes.fasta \
    --length-only \
    --min-length 500

# 提取短序列片段（<100bp，可能为片段）
biopytools subseq \
    -i assembly.fasta \
    -o short_fragments.fasta \
    --length-only \
    --max-length 100

# 提取中等长度序列（1kb-10kb）
biopytools subseq \
    -i all_sequences.fasta \
    -o medium_sequences.fasta \
    --length-only \
    --min-length 1000 \
    --max-length 10000
```

### 示例5：质控数据清理 | Example 5: Quality Control Data Cleaning

```bash
# 移除过短的低质量序列
biopytools subseq \
    -i raw_sequences.fasta \
    -o quality_sequences.fasta \
    --length-only \
    --min-length 200

# 提取特定标记基因
biopytools subseq \
    -i marker_genes.fasta \
    -o COX_genes.fasta \
    -p "COX1\|COX2\|COX3" \
    --pattern-type regex
```

### 示例6：批量处理工作流 | Example 6: Batch Processing Workflow

```bash
#!/bin/bash
# 批量处理多个FASTA文件

for fasta_file in *.fasta; do
    output_file="${fasta_file%.fasta}_extracted.fasta"

    # 提取所有基因序列
    biopytools subseq \
        -i "$fasta_file" \
        -p "gene" \
        -o "$output_file" \
        --ignore-case

    echo "Processed: $fasta_file -> $output_file"
done
```

## 📊 输出结果 | Output Results

### 输出文件格式 | Output File Format

输出为标准FASTA格式，包含提取的序列：

```fasta
>seq1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr1_seq3
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>gene1
TTTTAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCC
```

### 统计信息示例 | Statistics Example

```
📊 统计信息 | Statistics:
    📄 输入序列总数 | Total input sequences: 1000
    📋 请求提取ID数 | Requested IDs: 150
    ✅ 成功匹配ID数 | Successfully matched IDs: 148
    📤 成功提取序列数 | Successfully extracted sequences: 148
    📈 提取成功率 | Extraction success rate: 100.00%
```

### 日志文件 | Log Files

自动生成的日志文件包含：
- 📅 处理时间戳
- 📝 详细操作记录
- ⚠️ 警告和错误信息
- 📊 处理统计信息

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Python** (版本 3.7+)
- **BioPython** (版本 1.78+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 安装依赖 | Installing Dependencies

```bash
# 安装BioPython
pip install biopython

# 安装其他依赖
pip install click
```

### 硬件要求 | Hardware Requirements

- **CPU**: 单核即可，多核可提升处理速度
- **RAM**: 取决于输入FASTA文件大小
- **存储**: 需要足够空间存储输出文件和日志

## ⚠️ 注意事项 | Important Notes

1. **ID匹配**: 序列ID必须与FASTA文件中的ID完全匹配
2. **文件格式**: 确保FASTA文件格式正确
3. **内存使用**: 大文件处理时注意内存使用情况
4. **正则表达式**: 复杂的正则表达式可能影响处理速度
5. **编码格式**: 建议使用UTF-8编码的文本文件

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "找不到序列ID"错误**
```bash
# 检查ID是否完全匹配（包括大小写）
grep "seq1" input.fasta

# 使用忽略大小写选项
biopytools subseq -i input.fasta -l ids.txt -o output.fasta --ignore-case
```

**Q: "FASTA文件格式错误"**
```bash
# 检查FASTA文件格式
head -10 input.fasta
# 确保每条序列以 > 开头
```

**Q: 正则表达式不匹配**
```bash
# 测试正则表达式
echo "gene_001_test" | grep -E "^gene_\d{3}.*$"

# 使用简单的模式进行测试
biopytools subseq -i input.fasta -p "gene" -o test.fasta --pattern-type contains
```

**Q: 内存不足错误**
```bash
# 检查系统内存
free -h

# 分批处理大文件
split -l 1000 large_ids.txt small_ids_
for file in small_ids_*; do
    biopytools subseq -i large.fasta -l "$file" -o "output_$file.fasta"
done
```

**Q: 处理速度慢**
```bash
# 使用更简单的匹配模式
biopytools subseq -i input.fasta -p "prefix" -o output.fasta --pattern-type startswith

# 避免复杂正则表达式
biopytools subseq -i input.fasta -p "simple_pattern" -o output.fasta
```

## 📚 相关资源 | Related Resources

- [FASTA格式规范](https://en.wikipedia.org/wiki/FASTA_format)
- [正则表达式教程](https://www.regular-expressions.info/)
- [BioPython文档](https://biopython.org/wiki/Documentation)
- [生物序列处理最佳实践](https://www.ncbi.nlm.nih.gov/books/NBK279688/)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用BioPyTools：

```
Li, X. BioPyTools: A Python toolkit for bioinformatics analysis.
GitHub repository: https://github.com/lixiang117423/biopytools
```

以及相关的BioPython包：

```
Cock, P. J. A., et al. (2009). The Biopython project: an open-source toolkit for computational molecular biology.
Bioinformatics, 25(11), 1422-1423.
```