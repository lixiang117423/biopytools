# 基因组装配统计工具|Genome Assembly Statistics Tool

## 版本|Version: 1.0.0

## 功能描述|Function Description

基因组装配统计工具用于报告FASTA和FASTQ文件的序列长度统计信息，包括N50、N60、N70、N80、N90、N100等组装质量指标。|Genome Assembly Statistics Tool reports sequence length statistics from fasta and/or fastq files, including assembly quality metrics such as N50, N60, N70, N80, N90, N100.

## 功能特性|Features

- 支持FASTA和FASTQ格式文件|Support FASTA and FASTQ format files
- 计算N50、N60、N70、N80、N90、N100统计值|Calculate N50, N60, N70, N80, N90, N100 statistics
- 统计Gap数量和N碱基数|Count gaps and N bases
- 多种输出格式（默认、grep友好、tab分隔）|Multiple output formats (default, grep-friendly, tab-delimited)
- 自动生成CSV和Excel报告|Automatically generate CSV and Excel reports
- 支持最小长度过滤|Support minimum length filtering

## 安装|Installation

```bash
# 安装biopytools|Install biopytools
pip install biopytools
```

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 分析单个FASTA文件|Analyze single FASTA file
biopytools assembly-stats genome.fa

# 分析多个文件|Analyze multiple files
biopytools assembly-stats genome1.fa genome2.fa genome3.fa

# 分析FASTQ文件|Analyze FASTQ file
biopytools assembly-stats reads.fastq
```

### 高级选项|Advanced Options

```bash
# 设置最小长度过滤|Set minimum length filtering
biopytools assembly-stats -l 1000 genome.fa

# Tab分隔输出（适合脚本处理）|Tab-delimited output (suitable for script processing)
biopytools assembly-stats -t genome.fa

# Grep友好格式|Grep-friendly format
biopytools assembly-stats -s genome.fa

# 无header的Tab输出|Tab-delimited output without header
biopytools assembly-stats -u genome.fa

# 指定输出目录|Specify output directory
biopytools assembly-stats -o results/ genome.fa
```

### Python代码调用|Python Code Usage

```python
from biopytools.assembly_stats import AssemblyStatsConfig, AssemblyStatsRunner
from biopytools.assembly_stats.utils import AssemblyStatsLogger

# 创建配置|Create configuration
config = AssemblyStatsConfig(
    input_files=['genome.fa'],
    min_length=1000,
    output_dir='./results'
)

# 创建日志|Create logger
logger_manager = AssemblyStatsLogger('./results/assembly_stats.log')
logger = logger_manager.get_logger()

# 运行分析|Run analysis
runner = AssemblyStatsRunner(config, logger)
runner.run()
```

## 输出格式|Output Format

### 默认输出格式|Default Output Format

```
stats for genome.fa
sum = 23328019, n = 16, ave = 1458001.19, largest = 3291936
N50 = 1687656, n = 5
N60 = 1472805, n = 7
N70 = 1445207, n = 8
N80 = 1343557, n = 10
N90 = 1067971, n = 12
N100 = 5967, n = 16
N_count = 0
Gaps = 0
```

### Tab分隔格式|Tab-delimited Format

```bash
biopytools assembly-stats -t genome.fa
```

输出|Output:
```
File	Sum	N	Ave	Largest	N50	N60	N70	N80	N90	N100	N_count	Gaps
genome.fa	23328019	16	1458001.19	3291936	1687656	1472805	1445207	1343557	1067971	5967	0	0
```

## 输出字段说明|Output Field Description

| 字段|Field | 说明|Description |
|------|------|------|----------|
| sum|总和|所有序列的总长度|Total length of all sequences |
| n|序列数|序列的数量|Number of sequences |
| ave|平均长度|序列的平均长度|Average sequence length |
| largest|最大长度|最长的序列长度|Length of the longest sequence |
| N50|N50值|50%的组装长度包含在N50长度的序列中|50% of assembly is contained in sequences of length N50 or longer |
| N60|N60值|60%的组装长度包含在N60长度的序列中|60% of assembly is contained in sequences of length N60 or longer |
| N70|N70值|70%的组装长度包含在N70长度的序列中|70% of assembly is contained in sequences of length N70 or longer |
| N80|N80值|80%的组装长度包含在N80长度的序列中|80% of assembly is contained in sequences of length N80 or longer |
| N90|N90值|90%的组装长度包含在N90长度的序列中|90% of assembly is contained in sequences of length N90 or longer |
| N100|N100值|100%的组装长度包含在N100长度的序列中|100% of assembly is contained in sequences of length N100 or longer |
| N_count|N碱基数|所有未确定的N碱基总数|Total number of undetermined N bases |
| Gaps|Gap数|连续N碱基的Gap总数|Total number of gaps (consecutive Ns) |

## 报告文件|Report Files

工具会自动在输出目录生成以下报告文件|The tool automatically generates the following report files in the output directory:

- `assembly_stats.csv`: CSV格式的统计报告|Statistics report in CSV format
- `assembly_stats.xlsx`: Excel格式的统计报告|Statistics report in Excel format
- `assembly_stats.log`: 运行日志文件|Run log file

## 示例|Example

```bash
# 分析疟原虫参考基因组|Analyze Plasmodium falciparum reference genome
biopytools assembly-stats Pf3D7_v3.fasta
```

输出|Output:
```
stats for Pf3D7_v3.fasta
sum = 23328019, n = 16, ave = 1458001.19, largest = 3291936
N50 = 1687656, n = 5
N60 = 1472805, n = 7
N70 = 1445207, n = 8
N80 = 1343557, n = 10
N90 = 1067971, n = 12
N100 = 5967, n = 16
N_count = 0
Gaps = 0
```

解释|Interpretation:
- 总长度为23,328,019 bp，共16条序列|Total length is 23,328,019 bp across 16 sequences
- 平均序列长度为1,458,001 bp|Average sequence length is 1,458,001 bp
- 最长序列为3,291,936 bp|Longest sequence is 3,291,936 bp
- N50为1,687,656 bp，意味着50%的组装（约11.6 Mb）包含在5条序列中|N50 is 1,687,656 bp, meaning 50% of the assembly (~11.6 Mb) is contained in 5 sequences
- 没有未确定的N碱基|No undetermined N bases

## 参数说明|Parameter Description

| 参数|Parameter | 简写|Short | 类型|Type | 默认值|Default | 说明|Description |
|--------|--------|------|------|----------|------|
| files|--files|-|str|必需|Required|输入FASTA/FASTQ文件列表|List of input FASTA/FASTQ files |
| min-length|--min-length|-l|int|1|最小序列长度过滤|Minimum sequence length cutoff |
| grep-friendly|--grep-friendly|-s|flag|False|Grep友好输出格式|Print grep-friendly output |
| tab-delimited|--tab-delimited|-t|flag|False|Tab分隔输出|Print tab-delimited output |
| no-header|--no-header|-u|flag|False|Tab分隔输出且无header|Print tab-delimited output without header |
| output-dir|--output-dir|-o|string|./assembly_stats_output|输出目录|Output directory |

## 故障排除|Troubleshooting

### 问题1: 文件格式错误|Issue 1: File Format Error

**错误信息|Error Message:**
```
无法解析文件格式|Cannot parse file format: xxx.xxx
```

**解决方法|Solution:**
确保输入文件是标准FASTA或FASTQ格式|Ensure input file is in standard FASTA or FASTQ format

### 问题2: 没有有效序列|Issue 2: No Valid Sequences

**错误信息|Error Message:**
```
未找到有效序列|No valid sequences found in: xxx.fa
```

**解决方法|Solution:**
- 检查最小长度设置是否过高|Check if minimum length setting is too high
- 使用`-l`参数降低最小长度阈值|Use `-l` parameter to lower the minimum length threshold

## 版本历史|Version History

- **v1.0.0** (2025-01-04): 初始版本|Initial version
  - 支持FASTA/FASTQ文件分析|Support FASTA/FASTQ file analysis
  - 计算Nx统计值|Calculate Nx statistics
  - 生成CSV和Excel报告|Generate CSV and Excel reports

## 许可证|License

MIT License

## 作者|Author

biopytools development team

## 联系方式|Contact

如有问题或建议，请提交issue|For questions or suggestions, please submit an issue at:
https://github.com/yourusername/biopytools/issues
