# 基因组Gap统计模块

**统计基因组FASTA文件中Gap的位置和长度 | Statistics Gap Positions and Lengths in Genome FASTA Files**

## 功能概述 | Overview

基因组Gap统计模块是一个专门用于分析基因组组装中Gap（连续N）区域的工具。它能够精确识别每条序列中所有Gap的位置、长度，并提供详细的统计信息。支持单文件处理和批量处理模式，适用于基因组组装质量评估和比较分析。

## 主要特性 | Key Features

- **精确Gap识别**: 使用正则表达式精确匹配连续N，支持自定义最小长度
- **1-based坐标**: 输出使用1-based坐标系统，便于与基因组浏览器对照
- **单文件模式**: 处理单个FASTA文件，结果输出到文件或终端
- **批量处理模式**: 处理整个文件夹的FASTA文件，合并结果到一个输出文件
- **详细统计**: 提供总gap数量、总长度、各序列统计等详细信息
- **灵活输出**: 支持TSV格式输出，便于导入Excel或数据分析工具
- **自动样品识别**: 批量模式下自动从文件名提取样品名称

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 统计单个基因组文件的gap
biopytools gap-stat -i genome.fa -o gaps.txt

# 输出到终端
biopytools gap-stat -i genome.fa

# 只统计长度>=100的gap
biopytools gap-stat -i genome.fa --min-n 100 -o long_gaps.txt
```

### 高级用法 | Advanced Usage

```bash
# 批量处理文件夹中所有基因组文件
biopytools gap-stat -i genomes_folder/ -o all_gaps.txt

# 批量处理并只统计长gap
biopytools gap-stat -i genomes_folder/ -o long_gaps.txt --min-n 50
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTA文件或文件夹路径 | `-i genome.fa` |

### 可选参数 | Optional Parameters

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-o, --output` | 输出文件路径（不指定则输出到终端） | 输出到终端 | `-o gaps.txt` |
| `--min-n` | 最小连续N数量 | 1 | `--min-n 100` |

## 输出文件 | Output Files

### 单文件模式输出格式

```tsv
seq_name	start	end	gap_length
chr1	1000	1050	51
chr1	2500	2505	6
chr2	500	510	11
```

- **seq_name**: 序列名称
- **start**: Gap起始位置（1-based坐标）
- **end**: Gap终止位置（1-based坐标）
- **gap_length**: Gap长度（bp）

### 批量模式输出格式

```tsv
sample	seq_name	start	end	gap_length
sample1	chr1	1000	1050	51
sample1	chr2	500	510	11
sample2	chr1	2000	2010	11
```

增加了**sample**列标识样品来源。

### 统计信息

程序执行后会在stderr输出汇总信息：

```
============================================================
汇总|Summary
============================================================
总gap数量|Total gaps: 125
总gap长度|Total gap length: 15,320 bp
包含gap的序列数|Sequences with gaps: 12
============================================================
```

批量模式时还会显示各样品的独立统计。

## 使用场景 | Use Cases

### 场景1：评估基因组组装质量

```bash
# 查看组装中的gap情况
biopytools gap-stat -i assembly.fa -o assembly_gaps.txt

# 根据结果评估：
# - Gap总数越少越好
# - Gap长度分布应该集中在短gap
# - 长Gap可能需要进一步处理
```

### 场景2：比较不同组装版本

```bash
# 将不同版本的组装文件放在同一文件夹
mkdir assembly_comparison
cp v1_assembly.fa assembly_comparison/
cp v2_assembly.fa assembly_comparison/

# 批量统计
biopytools gap-stat -i assembly_comparison/ -o comparison_gaps.txt

# 在Excel中对比不同版本的gap数量和分布
```

### 场景3：只关注长Gap

```bash
# 只统计长度>=100bp的gap，这些可能影响基因注释
biopytools gap-stat -i genome.fa --min-n 100 -o long_gaps.txt
```

### 场景4：监控Gap填充效果

```bash
# Gap填充前后对比
biopytools gap-stat -i before_gapfill.fa -o before.txt
biopytools gap-stat -i after_gapfill.fa -o after.txt

# 比较两个文件：
# - Gap数量是否减少
# - Gap总长度是否缩短
# - 剩余gap的分布情况
```

## 技术细节 | Technical Details

### Gap识别算法

使用正则表达式 `N{min_n,}` 匹配连续N：
- 默认匹配所有连续N（min_n=1）
- 可设置min_n过滤短gap
- 大小写不敏感（N/n都匹配）

### 坐标系统

- 使用**1-based坐标系统**
- start: Gap第一个N的起始位置
- end: Gap最后一个N的终止位置
- gap_length = end - start + 1

### 批量处理逻辑

1. 扫描文件夹，识别所有FASTA文件（.fa, .fasta, .fna等）
2. 从文件名提取样品名称（去掉扩展名）
3. 逐个处理文件，收集所有gap信息
4. 在每条记录前添加sample列
5. 合并所有结果到单个输出文件
6. 生成各样品的独立统计信息

## 输出示例 | Output Examples

### 示例1：单文件输出

**输入文件 (test.fa)**:
```fasta
>chr1
ATCGATCGATNNNNNATCGATCGATCGATCGNNNNNATCG
>chr2
ATCGATCGATCGATCGATCGATCGATCGATCG
```

**输出**:
```tsv
seq_name	start	end	gap_length
chr1	11	15	5
chr1	32	36	5
```

### 示例2：批量处理

**输入文件夹结构**:
```
genomes/
  ├── sample1.fa
  ├── sample2.fa
  └── sample3.fa
```

**命令**:
```bash
biopytools gap-stat -i genomes/ -o all_gaps.txt
```

**输出 (all_gaps.txt)**:
```tsv
sample	seq_name	start	end	gap_length
sample1	chr1	100	150	51
sample2	chr1	200	210	11
sample3	chr2	500	520	21
```

## 常见问题 | FAQ

**Q: 1-based坐标是什么意思？**

A: 1-based坐标是指序列第一个位置从1开始计数，而不是0。这与基因组浏览器（如IGV、UCSC）和常用工具（如BED格式）一致，便于对照查看。

**Q: 如何判断gap统计结果是否合理？**

A:
- **高质量组装**: gap数量少（通常<100），gap长度短（大部分<100bp）
- **中等质量组装**: gap数量中等（100-1000），有一些长gap
- **低质量组装**: gap数量多（>1000），很多长gap（>1000bp）

**Q: 批量处理时如何识别不同样品的结果？**

A: 输出文件的第一列是**sample**列，显示从文件名提取的样品名称。可以根据这一列筛选、排序或统计分析。

**Q: 为什么统计信息输出到stderr而不是文件？**

A: 这样设计是为了让汇总信息不影响结果文件的导入。如果需要保存统计信息到文件，可以重定向stderr：
```bash
biopytools gap-stat -i genome.fa 2> statistics.txt
```

## 性能建议 | Performance Recommendations

| 基因组大小 | 预计运行时间 | 内存使用 |
|-----------|-------------|---------|
| <100 Mb | <1秒 | <100 MB |
| 100 Mb - 1 Gb | 1-5秒 | <500 MB |
| 1 Gb - 10 Gb | 5-30秒 | <1 GB |
| >10 Gb | 30秒-2分钟 | <2 GB |

*注：批量处理时间与文件数量成正比*

## 依赖要求 | Dependencies

### 必需依赖 | Required Dependencies

- Python 3.7+
- 标准库模块：re, sys, argparse, pathlib, collections

### 无外部依赖 | No External Dependencies

本模块仅使用Python标准库，无需安装额外的第三方包。

## 与其他工具的配合 | Integration with Other Tools

### 1. 配合Gap填充工具

```bash
# 填充前统计
biopytools gap-stat -i assembly.fa -o before.txt

# 使用TGS-GapCloser填充
biopytools tgsgapcloser -s assembly.fa -t ont -ir reads.fa -o filled

# 填充后统计
biopytools gap-stat -i assembly.gapcloser.fa -o after.txt

# 对比效果
diff before.txt after.txt
```

### 2. 导入Excel进行可视化

```bash
# 生成TSV文件
biopytools gap-stat -i genome.fa -o gaps.tsv

# 在Excel中：
# 1. 数据选项卡 -> 获取数据 -> 从文本/CSV
# 2. 选择制表符分隔
# 3. 创建数据透视表分析gap分布
```

### 3. 与agp2table配合

```bash
# AGP转表格查看scaffold结构
biopytools agp2table -i assembly.agp -o assembly_table.txt

# Gap统计查看gap详情
biopytools gap-stat -i assembly.fa -o gaps.txt

# 综合评估组装质量
```

## 引用 | Citation

如果使用本模块，请引用：

> BioPyTools Gap Statistics Module v1.0.0
> 仓库地址: https://github.com/your-org/biopytools

## 版本历史 | Version History

- **v1.0.0** (2026-03-20): 初始版本
  - 支持单文件和批量处理模式
  - 1-based坐标系统
  - 详细统计信息
  - 符合BioPyTools开发规范v2.11
