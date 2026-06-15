# 基因组挂载率统计模块

**计算FASTA文件中序列的挂载率 | Calculate Genome Mount Rate from FASTA Files**

## 功能概述 | Overview

基因组挂载率统计模块用于计算FASTA文件中前N条（或最长N条）序列占总基因组长度的百分比。该工具常用于评估基因组组装质量，例如统计主要染色体（较长的序列）占总基因组的比例。

## 主要特性 | Key Features

- **灵活的统计方式**: 支持按文件顺序或按长度排序统计前N条序列
- **推荐排序模式**: 排序后可计算最长N条序列的占比，更准确反映基因组挂载情况
- **简洁高效**: 纯Python实现，无需外部依赖
- **详细输出**: 显示总序列数、基因组大小、目标序列长度和挂载率百分比
- **标准日志**: 符合BioPyTools开发规范的日志输出格式

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 计算前10条序列的挂载率（按文件顺序）
biopytools genome-mount-rate -i genome.fa -n 10

# 计算最长10条序列的挂载率（推荐）
biopytools genome-mount-rate -i genome.fa -n 10 --sort
```

### 高级用法 | Advanced Usage

```bash
# 统计主要染色体数量（假设有20条染色体）
biopytools genome-mount-rate -i assembled_genome.fa -n 20 --sort

# 检查 scaffold 挂载情况
biopytools genome-mount-rate -i scaffolds.fa -n 100 --sort
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入的FASTA文件路径 | `-i genome.fa` |
| `-n, --number` | 要计算的序列数量 | `-n 10` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sort` | `False` | 按长度从大到小排序后计算（计算最长N条的占比） |

## 使用示例 | Usage Examples

### 示例1：评估二倍体基因组组装质量 | Example 1: Assess Diploid Genome Assembly Quality

```bash
# 对于二倍体生物（如大多数动物），统计最长20条序列（10对染色体）
biopytools genome-mount-rate -i mammal_genome.fa -n 20 --sort
```

**输出解释**:
- 挂载率接近100%说明组装质量很高，主要染色体序列完整
- 挂载率较低说明存在大量未挂载的scaffold或contig

### 示例2：评估植物基因组组装质量 | Example 2: Assess Plant Genome Assembly Quality

```bash
# 对于植物（可能有多倍体），统计最长30-50条序列
biopytools genome-mount-rate -i plant_genome.fa -n 50 --sort
```

### 示例3：比较不同组装版本 | Example 3: Compare Different Assembly Versions

```bash
# 比较v1.0和v2.0版本的挂载率
biopytools genome-mount-rate -i genome_v1.0.fa -n 20 --sort
biopytools genome-mount-rate -i genome_v2.0.fa -n 20 --sort
```

### 示例4：按文件顺序统计 | Example 4: Statistic by File Order

```bash
# 如果FASTA文件已按重要性排序，可直接统计
biopytools genome-mount-rate -i sorted_genome.fa -n 10
```

## 输出结果 | Output Results

### 输出格式 | Output Format

```
----------------------------------------
总序列数|Total sequences: 1500
总基因组大小|Total genome size: 1,234,567,890 bp
统计目标|Target: 最长|longest 10 条序列|sequences
目标序列总长|Target sequences total length: 1,200,000,000 bp
----------------------------------------
占比|Mount rate: 97.20%
----------------------------------------
```

### 结果解读 | Result Interpretation

| 挂载率范围 | 含义 | 建议 |
|-----------|------|------|
| **95% - 100%** | 优秀 | 主要序列完整，组装质量很高 |
| **80% - 95%** | 良好 | 主要序列基本完整，有少量未挂载序列 |
| **50% - 80%** | 中等 | 存在较多未挂载序列，可能需要scaffold构建 |
| **< 50%** | 较低 | 组装碎片化严重，需要改进组装策略 |

## 应用场景 | Applications

### 1. 基因组组装质量评估

在基因组组装完成后，计算主要染色体序列的挂载率是评估组装质量的重要指标。

### 2. 组装版本比较

比较不同组装版本或不同组装软件的结果，选择挂载率更高的版本。

### 3. 染色体级别组装验证

验证是否成功将scaffold挂载到染色体级别。

### 4. 基因组注释前评估

在开始基因组注释前，评估组装的完整性。

## 技术细节 | Technical Details

### 算法说明

1. 读取FASTA文件，统计每条序列的长度
2. 如果启用`--sort`，按长度从大到小排序
3. 取前N条序列，计算其总长度
4. 计算前N条序列长度占总长度的百分比

### 性能说明

- 时间复杂度: O(n) 或 O(n log n)（启用排序时）
- 空间复杂度: O(n)
- n为序列数量

### 限制

- 仅支持标准FASTA格式
- 序列名称以`>`开头
- 空行会被自动跳过

## 常见问题 | FAQ

### Q1: 为什么要使用`--sort`参数？

**A**: 使用`--sort`参数可以计算最长N条序列的占比，这在基因组组装质量评估中更有意义。例如，如果你想统计主要染色体的占比，它们通常是最长的序列。

### Q2: 如何选择N值？

**A**: N值应该根据物种的染色体数量选择：
- 二倍体动物: N = 染色体对数 x 2（如人类：N=23x2=46）
- 多倍体植物: N = 染色体组数 x 单倍体染色体数
- 如果不确定，可以尝试多个N值观察挂载率变化

### Q3: 挂载率低应该怎么办？

**A**: 挂载率低通常意味着：
1. 存在大量未挂载的scaffold
2. 组装碎片化严重
3. 可能需要使用遗传图谱或Hi-Fi数据辅助挂载

### Q4: 支持压缩格式吗？

**A**: 当前版本仅支持未压缩的FASTA格式。如需处理压缩文件，请先解压。

### Q5: 可以统计所有序列吗？

**A**: 可以，但挂载率会始终为100%。建议设置合理的N值以获得有意义的结果。

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Python**: 3.7+
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 硬件建议 | Hardware Recommendations

- **RAM**: 根据FASTA文件大小，通常1GB足够
- **存储**: 仅需读取输入文件，无额外存储需求

## 相关资源 | Related Resources

- [FASTA格式规范](https://en.wikipedia.org/wiki/FASTA_format)
- [基因组组装质量评估](https://www.ncbi.nlm.nih.gov/pmc/articles/PMCPMC4545759/)
- [BioPyTools开发规范](/share/org/YZWL/yzwl_lixg/software/scripts/develop_python_guides.md)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用BioPyTools：

```
BioPyTools: A comprehensive bioinformatics toolkit
https://github.com/your-repo/biopytools
```
