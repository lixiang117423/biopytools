# Kmer矩阵比较工具

**专业的kmer矩阵比较分析工具 | Professional Kmer Matrix Comparison Tool**

## 功能概述 | Overview

Kmer矩阵比较工具是一个高效的kmer矩阵文件比较分析工具，用于比较两个kmer矩阵文件，找出特有的kmer并计算窗口统计信息。使用流式处理技术，不会将整个文件载入内存，适用于大规模基因组数据的比较分析。

## 主要特性 | Key Features

- **流式处理**: 不会将整个文件载入内存，只载入kmer序列集合
- **特有kmer识别**: 通过比较kmer序列找出只在某个文件中出现的kmer
- **FoundAs标注**: 利用已有的FoundAs字段（forward/reverse/not_found）
- **窗口统计**: 按指定窗口大小统计unique kmer在各样品中的0/1比例
- **每个样品独立计算**: 每个样品的0/1比例单独计算
- **双文件输出**: 分别生成两个文件的统计结果

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 比较两个kmer矩阵文件
biopytools kmertools compare \
    -f1 hap1_group2_chr12_matrix.txt \
    -f2 hap2_group9_chr12_matrix.txt \
    -o chr12_comparison

# 指定窗口大小
biopytools kmertools compare \
    -f1 matrix1.txt \
    -f2 matrix2.txt \
    -o comparison \
    -w 50000
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-f1, --file1` | 第一个kmer矩阵文件路径 | `-f1 matrix1.txt` |
| `-f2, --file2` | 第二个kmer矩阵文件路径 | `-f2 matrix2.txt` |
| `-o, --output-prefix` | 输出文件前缀 | `-o comparison` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-w, --window-size` | `100000` | 窗口大小（行数）|Window size in lines |

## 输入文件格式 | Input File Format

kmer矩阵文件格式（TSV格式）：

```tsv
KmerID	Sequence	FoundAs	Sample1	Sample2	Sample3	...
hap1_1_51	AAAAAAA	forward	1	1	0	...
hap1_2_52	AAAAAAA	reverse	1	0	1	...
hap1_3_53	ATTTTTT	not_found	NA	NA	NA	...
```

**列说明**：
- **KmerID**: kmer唯一标识符
- **Sequence**: kmer序列（51bp）
- **FoundAs**: 查找方式（forward/reverse/not_found）
- **SampleN**: 各样品中的丰度（0或1）

## 输出结果 | Output Results

### 输出文件 | Output Files

生成两个文件：
1. `{output_prefix}_file1_stats.txt` - 文件1的统计结果
2. `{output_prefix}_file2_stats.txt` - 文件2的统计结果

### 输出格式 | Output Format

```tsv
File_ID	Start	End	Sample1_ratio	Sample2_ratio	Sample3_ratio	...
hap1_group2_chr12_matrix	1	100000	0.9876	0.9543	0.9234	...
hap1_group2_chr12_matrix	100001	200000	0.9123	0.8901	0.8765	...
```

**列说明**：
- **File_ID**: 文件标识
- **Start**: 窗口起始行（从1开始）
- **End**: 窗口结束行
- **SampleN_ratio**: 该窗口内unique kmer在该样品中值为1的比例

## 使用示例 | Usage Examples

### 示例1：基本比较 | Example 1: Basic Comparison

```bash
# 比较两个单倍型的kmer矩阵
biopytools kmertools compare \
    -f1 hap1_group2_chr12_matrix.txt \
    -f2 hap2_group9_chr12_matrix.txt \
    -o chr12_hap_comparison
```

### 示例2：自定义窗口大小 | Example 2: Custom Window Size

```bash
# 使用50000行的窗口进行统计
biopytools kmertools compare \
    -f1 matrix1.txt \
    -f2 matrix2.txt \
    -o comparison_w50000 \
    -w 50000
```

## 分析流程 | Analysis Workflow

1. **提取kmer序列集合**（流式处理）
   - 从两个文件中分别提取所有kmer序列
   - 存储在内存中的set中

2. **找出特有序列**
   - 计算两个集合的差集
   - unique1 = seqs1 - seqs2
   - unique2 = seqs2 - seqs1

3. **流式处理文件1**
   - 逐行读取文件1
   - 检查每个kmer是否在unique1中
   - 统计每个窗口中unique kmer在各个样品中的0/1比例

4. **流式处理文件2**
   - 同样的流程处理文件2

5. **输出统计结果**
   - 分别生成两个文件的统计报告

## 注意事项 | Important Notes

1. **内存使用**: 只载入kmer序列集合到内存，不是整个文件
2. **FoundAs字段**: 只统计FoundAs不是not_found的kmer
3. **Unique定义**: unique是指kmer序列只在一个文件中出现
4. **流式处理**: 适合处理大文件（GB级别）
5. **样品比例**: 每个样品的0/1比例独立计算

## 性能建议 | Performance Recommendations

- **大文件处理**: 建议在计算节点上运行，不要在登录节点
- **窗口大小**: 默认100000行适合大多数情况
- **磁盘空间**: 确保输出目录有足够的磁盘空间

## 常见问题 | FAQ

**Q: 为什么要用sequence而不是kmer_id来比较？**
A: 因为同一个kmer序列可能在不同位置出现，会有不同的kmer_id。

**Q: FoundAs字段有什么用？**
A: FoundAs标识kmer是如何被找到的（forward/reverse/not_found），我们只统计非not_found的kmer。

**Q: 窗口大小应该如何选择？**
A: 默认100000行适合大多数情况。如果想要更精细的统计，可以减小窗口大小；如果想要更粗粒度的统计，可以增大窗口大小。

**Q: 为什么每个样品的比例是独立计算的？**
A: 因为不同样品可能有不同的缺失率，独立计算可以反映每个样品的真实情况。

## 相关资源 | Related Resources

- [Kmer工具集文档](kmertools.md)
- [Intersect工具文档](kmertools_intersect.md)
- [开发规范文档](../scripts/develop_python_guides.md)
