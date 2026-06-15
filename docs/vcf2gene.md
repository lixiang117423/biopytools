# VCF2Gene 变异基因注释工具

**专业的VCF变异基因区域注释工具 | Professional VCF Variant Gene Region Annotation Tool**

## 功能概述 | Overview

VCF2Gene变异基因注释工具是一个基于GFF基因组注释文件的VCF变异位点注释工具。它能够快速准确地判断每个变异位点位于外显子、内含子还是基因间区，并标注对应的基因ID，适用于各种基因组变异分析研究。

## 主要特性 | Key Features

- **精确定位**: 准确判断变异位点位于外显子/内含子/基因间区
- **高效处理**: 优化的算法设计，支持大规模VCF文件处理
- **灵活输入**: 支持标准VCF格式和GFF/GFF3格式文件
- **简洁输出**: 标准化的输出格式，便于下游分析
- **日志记录**: 详细的处理过程日志和进度追踪
- **压缩支持**: 自动识别并处理gzip压缩的VCF和GFF文件

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 变异基因注释
biopytools vcf2gene \
    -i variants.vcf \
    -g annotation.gff \
    -o annotated_variants.txt
```

### 使用压缩文件 | Using Compressed Files

```bash
# 使用gzip压缩的输入文件
biopytools vcf2gene \
    -i variants.vcf.gz \
    -g annotation.gff3.gz \
    -o annotated_variants.txt
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --vcf` | 输入VCF文件路径（支持.vcf.gz）| `-i variants.vcf` |
| `-g, --gff` | 输入GFF注释文件路径（支持.gff3.gz）| `-g annotation.gff` |
| `-o, --output` | 输出结果文件路径 | `-o annotated.txt` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `4` | 线程数（预留参数，暂未启用多线程） |

## 输入文件格式 | Input File Formats

### VCF变异文件 | VCF Variant File

标准VCF格式的变异调用结果：

```vcf
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
Chr1	1000	.	A	G	45.2	PASS	DP=30
Chr1	5000	.	T	C	62.8	PASS	DP=25
Chr2	10000	.	G	A	50.1	PASS	DP=20
```

### GFF注释文件 | GFF Annotation File

标准GFF3格式的基因组注释文件，需包含gene和exon特征：

```gff
##gff-version 3
Chr1	RefSeq	gene	1000	6000	.	+	.	ID=gene001;Name=Gene1
Chr1	RefSeq	exon	1000	1500	.	+	.	ID=exon1;Parent=gene001
Chr1	RefSeq	exon	2000	2500	.	+	.	ID=exon2;Parent=gene001
Chr1	RefSeq	exon	5000	5500	.	+	.	ID=exon3;Parent=gene001
Chr2	RefSeq	gene	8000	12000	.	-	.	ID=gene002;Name=Gene2
Chr2	RefSeq	exon	9000	9500	.	-	.	ID=exon4;Parent=gene002
```

**文件要求**:
- 标准GFF3格式
- 包含gene特征类型（用于定义基因范围）
- 包含exon特征类型（用于标注外显子区域）
- exon特征需要通过Parent属性或gene_id属性关联到基因

## 输出结果 | Output Results

### 输出格式 | Output Format

输出文件为制表符分隔的文本文件，包含以下列：

| 列名 | 描述 | 示例值 |
|------|------|--------|
| `Chr` | 染色体 | Chr1 |
| `Pos` | 变异位置 | 12345 |
| `Ref` | 参考碱基 | A |
| `Alt` | 变异碱基 | G |
| `Gene` | 基因ID或intergenic | gene001 |
| `Type` | 变异类型（exon/intron/NA）| exon |

### 输出示例 | Output Example

```
Chr	Pos	Ref	Alt	Gene	Type
Chr1	1200	A	G	gene001	exon
Chr1	3000	T	C	gene001	intron
Chr1	7000	G	A	intergenic	NA
Chr2	9200	C	T	gene002	exon
```

### 类型说明 | Type Description

| Type值 | 说明 | 优先级 |
|--------|------|--------|
| `exon` | 变异位于外显子区域 | 高 |
| `intron` | 变异位于基因范围内但不在外显子（即内含子）| 中 |
| `NA` | 变异位于基因间区 | 低 |

## 使用示例 | Usage Examples

### 示例1：基本注释 | Example 1: Basic Annotation

```bash
# 对植物基因组变异进行注释
biopytools vcf2gene \
    -i arabidopsis_variants.vcf \
    -g arabidopsis_annotation.gff3 \
    -o annotated_arabidopsis.txt
```

### 示例2：使用压缩文件 | Example 2: Using Compressed Files

```bash
# 处理压缩的VCF和GFF文件
biopytools vcf2gene \
    -i human_variants.vcf.gz \
    -g human_annotation.gff3.gz \
    -o annotated_human.txt
```

### 示例3：批量处理多个样本 | Example 3: Batch Processing

```bash
# 使用循环处理多个样本
for sample in sample1 sample2 sample3; do
    biopytools vcf2gene \
        -i ${sample}.vcf \
        -g annotation.gff3 \
        -o ${sample}_annotated.txt
done
```

## 算法说明 | Algorithm Description

### 注释逻辑 | Annotation Logic

工具采用以下优先级判断变异位置：

1. **外显子检查**: 首先检查变异是否位于任何外显子区域内
2. **基因范围检查**: 如果不在外显子，检查是否在基因起止位置范围内
3. **基因间区**: 如果既不在外显子也不在基因范围内，标记为基因间区

### 数据结构 | Data Structure

- **基因存储**: 使用字典按染色体组织基因信息
- **外显子存储**: 使用排序的列表按染色体存储外显子区间
- **区间查询**: 对已排序的外显子进行高效区间查询

## 性能说明 | Performance

### 处理速度 | Processing Speed

- 小规模VCF（< 1万个变异）: 秒级完成
- 中等规模VCF（1万-100万变异）: 分钟级完成
- 大规模VCF（> 100万变异）: 根据基因注释复杂度，可能需要更长时间

### 内存占用 | Memory Usage

内存使用主要取决于GFF文件的基因和外显子数量：
- 小型基因组（如拟南芥）: < 100MB
- 中型基因组（如水稻）: 100-500MB
- 大型基因组（如人类）: 500MB-2GB

## 注意事项 | Important Notes

1. **GFF格式要求**: 确保GFF文件包含gene和exon特征类型
2. **染色体命名一致性**: VCF和GFF文件的染色体命名必须完全一致
3. **坐标系统**: 使用1-based坐标系统（与GFF3标准一致）
4. **多ALT位点**: 对于有多个ALT的变异，会为每个ALT创建一行输出
5. **文件编码**: 输入文件应为UTF-8编码

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**问题1: "文件不存在"错误**
```bash
# 检查文件路径是否正确
ls -lh variants.vcf
ls -lh annotation.gff

# 使用绝对路径
biopytools vcf2gene \
    -i /full/path/to/variants.vcf \
    -g /full/path/to/annotation.gff \
    -o output.txt
```

**问题2: 染色体名称不匹配**
```bash
# 检查VCF中的染色体名称
bcftools view variants.vcf | grep -v "^#" | cut -f1 | sort -u

# 检查GFF中的染色体名称
grep "^>" annotation.gff | cut -f1 | sort -u

# 确保染色体名称一致（如Chr1 vs chr1 vs 1）
```

**问题3: 所有变异都标记为intergenic**
```bash
# 检查GFF文件格式是否正确
head -20 annotation.gff

# 确认包含gene和exon特征
grep -w "gene" annotation.gff | head
grep -w "exon" annotation.gff | head
```

**问题4: 内存占用过高**
```bash
# 对于非常大的基因组，考虑：
# 1. 按染色体分割VCF文件
# 2. 分批处理
# 3. 使用服务器或集群环境
```

## 技术细节 | Technical Details

### 实现语言 | Implementation Language

- Python 3.7+
- 标准库（无额外依赖）

### 依赖包 | Dependencies

- `gzip`: 处理压缩文件
- `argparse`: 命令行参数解析
- `logging`: 日志记录
- `pathlib`: 路径处理

### 代码架构 | Code Architecture

遵循BioPyTools标准模块化架构：
- `config.py`: 配置管理
- `utils.py`: 日志和工具函数
- `calculator.py`: 核心计算逻辑
- `main.py`: 命令行入口

## 版本历史 | Version History

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.0.0 | 2026-02-05 | 初始版本发布 |

## 许可证 | License

本项目采用MIT许可证 - 详见项目LICENSE文件

---

## 联系方式 | Contact

如有问题或建议，请联系：Xiang LI
