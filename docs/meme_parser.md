# MEME Parser MEME Motif发现和解析工具

## 功能概述|Overview

MEME Parser工具用于运行MEME软件进行motif发现，并解析MEME输出结果（XML/TXT格式）为结构化表格。

## 主要特性|Key Features

- 运行MEME软件进行motif发现
- 解析MEME XML格式输出
- 解析MEME TXT格式输出
- 输出TSV/CSV/Excel格式表格
- 仅解析模式（处理已有MEME输出）

## 快速开始|Quick Start

### 基本用法

```bash
# 运行MEME并解析结果（蛋白质序列）
biopytools meme-parser -i proteins.fa -protein

# 运行MEME并解析结果（DNA序列）
biopytools meme-parser -i sequences.fa -dna
```

### 仅解析模式

```bash
# 解析已有的MEME XML输出
biopytools meme-parser -i meme_out/meme.xml --parse-only

# 解析已有的MEME TXT输出
biopytools meme-parser -i meme_out/meme.txt --parse-only
```

### 高级用法

```bash
# 自定义MEME参数
biopytools meme-parser -i proteins.fa -protein \
    -nmotifs 15 \
    -minw 8 -maxw 60 \
    -mod zoops \
    -o my_motifs

# 指定MEME软件路径
biopytools meme-parser -i proteins.fa -protein \
    --meme-path /path/to/meme
```

## 参数说明|Parameters

### 必需参数|Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-file` | 输入FASTA文件或MEME输出文件(xml/txt) | `-i proteins.fa` |

### 运行模式|Run Mode

| 参数 | 描述 |
|------|------|
| `--parse-only` | 仅解析模式，不运行MEME |

### 输出配置|Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-prefix` | `meme_results` | 输出文件前缀 |
| `--output-dir` | `.` | 输出目录 |
| `--no-tsv` | - | 不输出TSV文件 |
| `--no-csv` | - | 不输出CSV文件 |
| `--no-excel` | - | 不输出Excel文件 |

### MEME软件路径|MEME Software Path

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--meme-path` | `meme` | MEME可执行文件路径 |

### MEME参数|MEME Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-protein` | False | 输入序列为蛋白质 |
| `-dna` | False | 输入序列为DNA |
| `-mod` | `zoops` | Motif分布模式（zoops/anr/oor） |
| `-nmotifs` | 10 | Motif数量 |
| `-minw` | 6 | 最小motif宽度 |
| `-maxw` | 50 | 最大motif宽度 |
| `-objfun` | `classic` | 目标函数（classic/de/ce/cd） |
| `-markov-order` | 0 | Markov链阶数 |

## 输出文件|Output Files

### 主要输出文件

| 文件 | 描述 |
|------|------|
| `{prefix}.tsv` | 主要结果TSV文件 |
| `{prefix}.csv` | 主要结果CSV文件 |
| `{prefix}.xlsx` | 主要结果Excel文件 |
| `{prefix}_meme_out/` | MEME软件输出目录 |

### 输出表格列|Output Table Columns

| 列名 | 描述 |
|------|------|
| sequence_id | 序列ID |
| length | 序列长度 |
| motif_id | Motif ID (如Motif-1) |
| start | Motif起始位置 |
| end | Motif结束位置 |
| width | Motif宽度 |
| strand | 链方向 (+/-) |
| pvalue | P值 |
| evalue | E值 |
| num_sites | Motif位点数 |

## 使用示例|Usage Examples

### 示例1：蛋白质motif分析

```bash
biopytools meme-parser -i nlr_proteins.fa -protein \
    -nmotifs 10 \
    -minw 6 -maxw 50 \
    -o nlr_motifs
```

输出：
```
2026-02-28 15:00:00.000 - INFO - ============================================================
2026-02-28 15:00:00.001 - INFO - MEME Parser流程开始|MEME Parser pipeline starting
2026-02-28 15:00:00.002 - INFO - 运行MEME|Running MEME
...
```

### 示例2：解析已有MEME输出

```bash
# 假设已有MEME输出在meme_out目录
biopytools meme-parser -i meme_out/meme.xml --parse-only -o parsed_motifs
```

### 示例3：DNA motif分析

```bash
biopytools meme-parser -i promoter_sequences.fa -dna \
    -nmotifs 5 \
    -minw 8 -maxw 20 \
    -mod anr \
    -o promoter_motifs
```

## 输出说明|Output Description

### Motif分布模式|Motif Distribution Mode

- **zoops**: Zero or one occurrence per sequence（每个序列零次或一次）
- **anr**: Any number of repetitions（任意重复次数）
- **oor**: One occurrence per sequence（每个序列一次）

### 目标函数|Objective Function

- **classic**: 经典目标函数
- **de**: Differential Enrichment
- **ce**: Central Enrichment
- **cd**: Central Depletion

### 输出表格示例|Output Table Example

```
sequence_id	length	motif_id	start	end	width	strand	pvalue	evalue	num_sites
seq1		500	Motif-1	50	75	25	+	1e-10	1e-5	2
seq1		500	Motif-2	200	220	20	+	1e-8	1e-3	1
seq2		450	Motif-1	30	55	25	+	1e-12	1e-6	3
```

## 注意事项|Important Notes

1. 输入文件可以是FASTA文件（运行模式）或MEME输出文件（解析模式）
2. 蛋白质序列使用`-protein`，DNA序列使用`-dna`
3. MEME运行时间取决于序列数量和长度
4. 解析模式不需要MEME软件
5. XML格式解析比TXT格式更准确

## 故障排除|Troubleshooting

### 错误：MEME软件未找到

```
Configuration error:
MEME software not found or not executable: meme
```

**解决方法**：
```bash
# 指定MEME完整路径
biopytools meme-parser -i proteins.fa -protein \
    --meme-path /share/org/YZWL/yzwl_lixg/miniforge3/envs/meme_v.5.5.9/bin/meme
```

### 错误：输入文件不存在

```
Parameter Error: File or directory does not exist: proteins.fa
```

**解决方法**：检查输入文件路径是否正确

### 警告：Excel保存失败

```
WARNING: openpyxl not installed, skipping Excel output
```

**解决方法**：
```bash
pip install openpyxl
```

## 性能参考|Performance Reference

| 序列数量 | 序列长度 | 预计时间 |
|---------|---------|---------|
| 10 | 500 aa | ~2分钟 |
| 50 | 500 aa | ~10分钟 |
| 100 | 500 aa | ~30分钟 |
| 500 | 500 aa | ~2小时 |

## MEME软件安装|MEME Software Installation

```bash
# 使用conda安装
conda install -c bioconda meme

# 或使用已安装的MEME
biopytools meme-parser -i proteins.fa -protein \
    --meme-path /path/to/your/meme
```

## 相关链接|Related Links

- MEME Suite: http://meme-suite.org/
- MEME文档: http://meme-suite.org/doc/meme.html
