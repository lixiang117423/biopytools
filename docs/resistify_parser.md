# Resistify Parser Resistify结果解析工具

## 功能概述|Overview

Resistify Parser工具用于解析和处理Resistify软件的输出结果，将多个TSV文件合并为单一输出文件，支持序列提取和数据筛选。

## 主要特性|Key Features

- 解析results.tsv、domains.tsv、annotations.tsv文件
- 合并domain信息到主表格
- 按分类、长度、LRR长度筛选数据
- 提取NLR和NB-ARC序列
- 输出CSV和Excel格式

## 快速开始|Quick Start

### 基本用法

```bash
# 默认解析当前目录的Resistify输出
biopytools resistify-parser -i resistify_output/

# 指定输出前缀
biopytools resistify-parser -i resistify_output/ -o my_results
```

### 高级用法

```bash
# 提取NLR序列
biopytools resistify-parser -i resistify_output/ --extract-nlr

# 按分类筛选（如只保留TNL类）
biopytools resistify-parser -i resistify_output/ --filter-classification TN

# 按长度筛选
biopytools resistify-parser -i resistify_output/ --min-length 500 --max-length 3000

# 完整功能组合
biopytools resistify-parser -i resistify_output/ \
    --extract-nlr --extract-nbarc \
    --filter-classification CNL \
    --min-length 1000 \
    --min-lrr-length 200 \
    -o filtered_results
```

## 参数说明|Parameters

### 必需参数|Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-dir` | Resistify输出目录 | `-i resistify_output/` |

### 输出配置|Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-prefix` | `resistify_results` | 输出文件前缀 |
| `--output-dir` | `.` | 输出目录 |
| `--no-csv` | False | 不输出CSV文件 |
| `--no-excel` | False | 不输出Excel文件 |

### 序列提取选项|Sequence Extraction Options

| 参数 | 描述 |
|------|------|
| `--extract-nlr` | 提取NLR序列到FASTA文件 |
| `--extract-nbarc` | 提取NB-ARC序列到FASTA文件 |

### 筛选选项|Filtering Options

| 参数 | 描述 | 示例 |
|------|------|------|
| `--filter-classification` | 按分类筛选（如TN, CNL, NL, TNL等） | `--filter-classification CNL` |
| `--min-length` | 最小序列长度（aa） | `--min-length 500` |
| `--max-length` | 最大序列长度（aa） | `--max-length 3000` |
| `--min-lrr-length` | 最小LRR长度 | `--min-lrr-length 200` |

### 其他选项|Other Options

| 参数 | 描述 |
|------|------|
| `--include-motifs` | 包含motifs详情 |

## 输入文件要求|Input File Requirements

工具会在输入目录中查找以下文件：

- **results.tsv** (必需): Resistify主结果文件
- **domains.tsv** (必需): Domain信息文件
- **annotations.tsv** (可选): 注释信息文件
- **nlr.fasta** (可选): NLR序列文件（使用--extract-nlr时需要）
- **nbarc.fasta** (可选): NB-ARC序列文件（使用--extract-nbarc时需要）

## 输出文件|Output Files

### 主要输出文件

| 文件 | 描述 |
|------|------|
| `{prefix}.csv` | 主要结果CSV文件 |
| `{prefix}.xlsx` | 主要结果Excel文件 |
| `{prefix}_nlr.fasta` | 提取的NLR序列（使用--extract-nlr） |
| `{prefix}_nbarc.fasta` | 提取的NB-ARC序列（使用--extract-nbarc） |

### 主要输出列

| 列名 | 描述 |
|------|------|
| Sequence | 序列ID |
| Length | 序列长度（aa） |
| LRR_Length | LRR域长度 |
| Domains | 包含的domain列表 |
| Classification | NLR分类（TNL, CNL, NL等） |
| NBARC_motifs | NB-ARC motifs信息 |
| has_TIR | 是否包含TIR域 |
| has_NB_ARC | 是否包含NB-ARC域 |
| has_CC | 是否包含CC域 |
| has_LRR | 是否包含LRR域 |
| has_RPW8 | 是否包含RPW8域 |

## 使用示例|Usage Examples

### 示例1：基本解析

```bash
biopytools resistify-parser -i /path/to/resistify_output/
```

输出：
```
2026-02-28 10:00:00.000 - INFO - ============================================================
2026-02-28 10:00:00.001 - INFO - 开始Resistify结果解析|Starting Resistify results parsing
2026-02-28 10:00:00.002 - INFO - ============================================================
2026-02-28 10:00:00.003 - INFO - 解析results.tsv|Parsing results.tsv
2026-02-28 10:00:00.050 - INFO - 解析完成，共245条NLR基因|Parsing completed, 245 NLR genes
...
```

### 示例2：筛选特定分类

```bash
# 只保留TNL类NLR基因
biopytools resistify-parser -i resistify_output/ \
    --filter-classification TNL \
    -o tnl_only
```

### 示例3：完整流程（解析+筛选+序列提取）

```bash
biopytools resistify-parser -i resistify_output/ \
    --filter-classification CNL \
    --min-length 1000 \
    --max-length 4000 \
    --min-lrr-length 200 \
    --extract-nlr \
    --extract-nbarc \
    -o cnl_filtered
```

## 输出说明|Output Description

### NLR分类|NLR Classification

Resistify识别的主要NLR类型：

- **TNL**: TIR-NB-ARC-LRR
- **CNL**: CC-NB-ARC-LRR
- **NL**: NB-ARC-LRR
- **TN**: TIR-NB-ARC
- **CN**: CC-NB-ARC
- **N**: NB-ARC only
- **其他**: 其他组合

### Domain标记|Domain Markers

- **has_TIR**: Toll/Interleukin-1 receptor domain
- **has_NB_ARC**: Nucleotide-binding adapter shared by APAF-1, R proteins, and CED-4
- **has_CC**: Coiled-coil domain
- **has_LRR**: Leucine-rich repeat
- **has_RPW8**: RPW8 domain

## 注意事项|Important Notes

1. 输入目录必须包含results.tsv和domains.tsv
2. 提取序列功能需要相应的FASTA文件
3. Excel输出需要openpyxl包
4. 分类筛选支持部分匹配（如"TN"会匹配TN、TNL、TNX等）
5. 长度筛选为闭区间（包含边界值）

## 故障排除|Troubleshooting

### 错误：输入目录不存在

```
错误: Directory does not exist: /path/to/dir
```

**解决方法**：检查输入目录路径是否正确

### 错误：必需文件不存在

```
配置错误:
Required file not found: results.tsv
```

**解决方法**：确保输入目录包含所有必需文件

### 警告：Excel保存失败

```
WARNING: openpyxl not installed, skipping Excel output
```

**解决方法**：
```bash
pip install openpyxl
```

## 性能参考|Performance Reference

| NLR基因数量 | 预计时间 |
|---------|---------|
| 100 | < 1秒 |
| 1,000 | ~2秒 |
| 10,000 | ~15秒 |
| 100,000 | ~2分钟 |
