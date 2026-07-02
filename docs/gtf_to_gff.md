# GTF到GFF转换工具|GTF to GFF Converter

## 概述|Overview

GTF到GFF转换工具(`gtf2gff`)用于将GTF格式的基因注释文件转换为clean的GFF3格式。

## 功能|Features

- ✅ 将BRAKER3输出的GTF文件转换为GFF3格式
- ✅ 清理属性字段，只保留必要的ID、Name、Parent
- ✅ 自动转换transcript为mRNA格式
- ✅ 生成规范的exon和CDS ID
- ✅ 可选移除intron特征
- ✅ 支持start_codon和stop_codon

## 安装|Installation

该工具是biopytools的一部分，随biopytools一起安装：

```bash
cd /path/to/biopytools
pip install -e .
```

## 快速开始|Quick Start

### 基本用法|Basic Usage

```bash
# 转换GTF到GFF
biopytools gtf2gff -i braker.gtf -o clean.gff
```

### 常用选项|Common Options

```bash
# 移除intron特征
biopytools gtf2gff -i braker.gtf -o clean.gff --remove-introns

# 保留所有原始属性
biopytools gtf2gff -i braker.gtf -o clean.gff --keep-all-attributes

# 不清理属性（保留原始格式）
biopytools gtf2gff -i braker.gtf -o clean.gff --no-clean
```

## 参数说明|Parameters

### 必需参数|Required Parameters

| 参数| 说明| 示例|
|------|------|------|
| `-i, --input` | 输入GTF文件| `braker.gtf` |
| `-o, --output` | 输出GFF文件| `clean.gff` |

### 可选参数|Optional Parameters

| 参数 | 默认值 | 说明|
|------|--------|------|
| `--remove-introns` | False | 移除intron特征|
| `--keep-all-attributes` | False | 保留所有原始属性|
| `--no-clean` | False | 不清理属性字段|
| `-t, --threads` | 12 | 线程数（当前版本未使用）|

## 输入输出格式|Input/Output Format

### 输入格式（GTF）|Input Format (GTF)

```gtf
Chr12	AUGUSTUS	gene	1741	1953	.	-	.	g1
Chr12	AUGUSTUS	transcript	1741	1953	0.99	-	.	g1.t1
Chr12	AUGUSTUS	CDS	1741	1953	0.99	-	0	transcript_id "g1.t1"; gene_id "g1";
Chr12	AUGUSTUS	exon	1741	1953	.	-	.	transcript_id "g1.t1"; gene_id "g1";
```

### 输出格式（GFF3）|Output Format (GFF3)

```gff
##gff-version 3
#!gff-spec-version 1.21
#!source-braker BRAKER3
Chr12	AUGUSTUS	gene	1741	1953	.	-	.	ID=g1;Name=g1
Chr12	AUGUSTUS	mRNA	1741	1953	0.99	-	.	ID=g1.mRNA1;Name=g1.mRNA1;Parent=g1
Chr12	AUGUSTUS	CDS	1741	1953	0.99	-	0	ID=g1.mRNA1.cds1;Parent=g1.mRNA1
Chr12	AUGUSTUS	exon	1741	1953	.	-	.	ID=g1.mRNA1.exon1;Parent=g1.mRNA1
```

## ID转换规则|ID Conversion Rules

### BRAKER格式 → GFF格式|BRAKER Format → GFF Format

| BRAKER ID | GFF ID | 说明|
|-----------|---------|------|
| `g1` | `g1` | 基因ID保持不变|
| `g1.t1` | `g1.mRNA1` | transcript转换为mRNA，序号格式调整|
| - | `g1.mRNA1.exon1` | exon ID自动生成|
| - | `g1.mRNA1.cds1` | CDS ID自动生成|

## 属性清理规则|Attribute Cleaning Rules

### 默认清理模式（clean_attributes=True）|Default Clean Mode

**gene**:
- 只保留：`ID`, `Name`

**mRNA**:
- 只保留：`ID`, `Name`, `Parent`

**exon/CDS**:
- 只保留：`ID`, `Parent`

### 保留所有属性模式（keep_all_attributes=True）|Keep All Attributes Mode

保留GTF文件中的所有原始属性字段。

## 完整工作流示例|Complete Workflow Example

### BRAKER3 + GTF2GFF

```bash
# 1. 运行BRAKER3注释
biopytools braker \
  --genome genome.fa \
  --species Ov53 \
  --prot_seq proteins.fa \
  --isoseq long_reads.fa \
  --rnaseq_dirs rnaseq_dir \
  --threads 64

# 2. 转换GTF到GFF
biopytools gtf2gff \
  -i output/04_braker_annotation/braker.gtf \
  -o output/clean.gff \
  --remove-introns

# 3. （可选）使用gff_renamer重命名ID
biopytools gff-renamer \
  -i output/clean.gff \
  -o output/final.gff \
  -p CDRT \
  -s Ov
```

## 注意事项|Notes

1. **intron特征**：默认保留intron，使用`--remove-introns`可以移除
2. **属性清理**：默认清理属性，使用`--keep-all-attributes`保留原始属性
3. **ID格式**：transcript自动转换为mRNA格式（`g1.t1` → `g1.mRNA1`）
4. **GFF3版本**：默认输出GFF3格式，文件头包含`##gff-version 3`

## 故障排除|Troubleshooting

### 问题：转换后的GFF文件为空|Problem: Empty GFF file

**可能原因**：
- 输入GTF文件格式不正确
- GTF文件使用的是非标准特征名称

**解决方案**：
```bash
# 检查GTF文件格式
head -20 input.gtf

# 确认特征类型
cut -f3 input.gtf | sort -u
```

### 问题：缺少某些特征|Problem: Missing features

**可能原因**：
- BRAKER输出中确实缺少这些特征
- 使用了`--remove-introns`导致intron被移除

**解决方案**：
```bash
# 检查原始GTF文件中的特征类型
cut -f3 input.gtf | sort -u
```

## 版本历史|Version History

- **v1.0.0** (2025-02-04): 初始版本
  - 支持GTF到GFF3转换
  - 属性清理功能
  - transcript到mRNA ID转换
