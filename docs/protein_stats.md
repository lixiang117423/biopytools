# Protein Stats 蛋白质理化性质分析工具

## 功能概述|Overview

Protein Stats工具用于计算蛋白质序列的理化性质，包括序列长度、分子量、等电点等。

## 主要特性|Key Features

- 计算序列长度
- 计算分子量(Molecular Weight)
- 计算等电点(Isoelectric Point, pI)
- 可选：氨基酸组成
- 可选：不稳定指数
- 可选：脂肪指数
- 可选：芳香性
- 多种输出格式（TSV/CSV/Excel）

## 快速开始|Quick Start

### 基本用法

```bash
# 默认计算：长度、分子量、等电点
biopytools protein-stats -i proteins.fa

# 指定输出文件
biopytools protein-stats -i proteins.fa -o stats.tsv
```

### 高级用法

```bash
# 计算氨基酸组成
biopytools protein-stats -i proteins.fa --aa-composition

# 计算所有性质
biopytools protein-stats -i proteins.fa --aa-composition --instability-index --gravy --aromaticity

# 输出Excel格式
biopytools protein-stats -i proteins.fa --output-format excel -o stats.xlsx
```

## 参数说明|Parameters

### 必需参数|Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --protein-fasta` | 蛋白序列FASTA文件 | `-i proteins.fa` |

### 输出配置|Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-file` | `protein_stats.tsv` | 输出文件路径 |
| `--output-format` | `tsv` | 输出格式：tsv/csv/excel |

### 基本计算选项|Basic Calculation Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-length` | False | 不计算序列长度 |
| `--no-mw` | False | 不计算分子量 |
| `--no-pi` | False | 不计算等电点 |

### 高级计算选项|Advanced Calculation Options

| 参数 | 描述 |
|------|------|
| `--aa-composition` | 计算氨基酸组成（20种氨基酸的百分比） |
| `--instability-index` | 计算不稳定指数 |
| `--gravy` | 计算脂肪指数（疏水性） |
| `--aromaticity` | 计算芳香性 |

## 输出文件|Output File

### 基本输出列（默认）

| 列名 | 描述 |
|------|------|
| id | 蛋白ID |
| length | 序列长度（aa） |
| mw_da | 分子量 |
| pi | 等电点 |
| acid_base | 酸碱性（acidic/basic/neutral）|
| 酸碱性 | 酸碱性（酸性/碱性/中性）|

### 高级输出列（使用--aa-composition等）

| 列名 | 描述 |
|------|------|
| aa_A, aa_C, ... | 各氨基酸百分比（20种） |
| instability_index | 不稳定指数 |
| gravy | 脂肪指数（正值：疏水，负值：亲水） |
| aromaticity | 芳香性（芳香氨基酸百分比） |

## 使用示例|Usage Examples

### 示例1：基本分析

```bash
biopytools protein-stats -i proteins.fa -o nbarc_stats.tsv
```

输出：
```
#id	length	mw_da	pi	acid_base	酸碱性
Ccu05G0121-mRNA-1	1825	204532.45	6.52	acidic	酸性
Ccu01G2607-mRNA-1	2567	287123.89	8.34	basic	碱性
```

### 示例2：完整分析

```bash
biopytools protein-stats -i proteins.fa \
    --aa-composition \
    --instability-index \
    --gravy \
    --aromaticity \
    -o complete_stats.xlsx \
    --output-format excel
```

### 示例3：只计算等电点

```bash
biopytools protein-stats -i proteins.fa --no-length --no-mw -o pi_only.tsv
```

## 输出说明|Output Description

### 分子量
- 单位：道尔顿
- 基于平均原子量计算

### 等电点
- 蛋白质净电荷为零时的pH值
- 精度：小数点后2位
- 计算方法：Bjellqvist方法

### 酸碱性
- 基于等电点判断|Based on isoelectric point
- **acidic**: pI < 7（酸性蛋白）
- **basic**: pI > 7（碱性蛋白）
- **neutral**: pI = 7（中性蛋白）

### 不稳定指数
- 范围：通常<40为稳定
- >40：可能不稳定

### Gravy
- > 0：疏水蛋白
- < 0：亲水蛋白
- 典型范围：-2到+2

### 芳香性
- 芳香氨基酸（F, W, Y）的百分比
- 范围：0-1

## 注意事项|Important Notes

1. 输入必须是蛋白序列（FASTA格式）
2. 氨基酸组成选项会增加20列数据
3. Excel格式需要openpyxl包
4. 大文件（>10000条序列）可能需要较长时间

## 性能参考|Performance Reference

| 序列数量 | 预计时间 |
|---------|---------|
| 100 | < 1秒 |
| 1,000 | ~3秒 |
| 10,000 | ~30秒 |
| 100,000 | ~5分钟 |
