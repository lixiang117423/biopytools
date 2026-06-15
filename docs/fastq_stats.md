# FASTQ文件统计工具|FASTQ File Statistics Tool

基于seqkit的高性能FASTQ文件统计工具，支持自动配对双末端文件，输出CSV/Excel格式
High-performance FASTQ file statistics tool based on seqkit, supports automatic paired-end file matching, outputs CSV/Excel format

## 功能特性|Features

- 支持单端和双末端FASTQ文件统计|Supports single-end and paired-end FASTQ file statistics
- 自动识别配对文件|Automatic paired file recognition
- 多线程处理|Multi-threaded processing
- 输出CSV或Excel格式|Outputs CSV or Excel format
- 提供详细的统计信息|Provides detailed statistics including:
  - reads数量|read count
  - 序列长度统计|sequence length statistics (min/max/mean)
  - 总碱基数|total bases
  - 文件大小（压缩和解压后）|file size (compressed and uncompressed)
  - R1/R2分别统计|separate R1/R2 statistics

## 依赖要求|Dependencies

- seqkit: `conda install -c bioconda seqkit`
- pandas和openpyxl (仅Excel输出时需要)|pandas and openpyxl (required only for Excel output): `pip install pandas openpyxl`

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 统计目录下所有FASTQ文件|Statistics for all FASTQ files in directory
biopytools fq-stats -i /data/fastq/ -o results.csv

# 使用模式匹配双末端文件|Use pattern matching for paired-end files
biopytools fq-stats -i /data/fastq/ -o results.csv -p "*_1.clean.fq.gz"

# 输出Excel格式|Output in Excel format
biopytools fq-stats -i /data/fastq/ -o stats.xlsx

# 指定线程数|Specify thread count
biopytools fq-stats -i /data/fastq/ -o results.csv -t 16
```

### 参数说明|Parameters

| 参数|Parameter | 说明|Description |
|-----------|------|-----------|
| `-i, --input` | 输入FASTQ文件或目录|Input FASTQ file or directory |
| `-o, --output` | 输出文件路径(.csv或.xlsx)|Output file path (.csv or .xlsx) |
| `-p, --pattern` | FASTQ文件匹配模式|FASTQ file matching pattern |
| `-t, --threads` | 线程数(默认:12)|Number of threads (default: 12) |

## 模式匹配说明|Pattern Matching

模式参数用于自动识别配对的R1和R2文件：
The pattern parameter is used to automatically identify paired R1 and R2 files:

- `*` 通配符代表样品名称|`*` wildcard represents sample name
- 常见的R1模式示例|Common R1 pattern examples:
  - `*_1.clean.fq.gz` → R2: `*_2.clean.fq.gz`
  - `*_R1.fastq.gz` → R2: `*_R2.fastq.gz`
  - `*.1.fq.gz` → R2: `*.2.fq.gz`

## 输出文件格式|Output Format

### CSV格式示例|CSV Format Example

```csv
sample_name,R1_file,R2_file,R1_reads,R1_min_length,R1_max_length,R1_mean_length,R1_total_bases,R1_file_size_formatted,R1_uncompressed_size_formatted,R2_reads,...
sample1,sample1_1.fq.gz,sample1_2.fq.gz,1000000,50,150,125,125000000,1.2GB,5.8GB,1000000,...
```

### Excel格式|Excel Format

Excel输出包含相同的列，格式化更美观：
Excel output contains the same columns with better formatting:

- **sample_name**: 样品名称|Sample name
- **R1_file**: R1文件名|R1 file name
- **R2_file**: R2文件名|R2 file name
- **R1_reads**: R1 reads数量|R1 read count
- **R1_min_length**: R1最小序列长度|R1 minimum sequence length
- **R1_max_length**: R1最大序列长度|R1 maximum sequence length
- **R1_mean_length**: R1平均序列长度|R1 mean sequence length
- **R1_total_bases**: R1总碱基数|R1 total bases
- **R1_file_size_formatted**: R1文件大小（压缩）|R1 file size (compressed)
- **R1_uncompressed_size_formatted**: R1文件大小（解压后）|R1 file size (uncompressed)
- **R2_***: R2相关统计|R2 related statistics
- **total_reads**: 总reads数|Total read count
- **total_bases**: 总碱基数|Total bases

## 示例|Examples

### 示例1: 统计单个文件|Example 1: Single File

```bash
biopytools fq-stats -i sample_R1.fastq.gz -o sample_stats.csv
```

### 示例2: 统计目录下所有文件|Example 2: All Files in Directory

```bash
biopytools fq-stats -i /data/fastq/ -o all_samples.xlsx
```

### 示例3: 使用模式匹配双末端文件|Example 3: Pattern Matching for Paired-end Files

```bash
biopytools fq-stats -i /data/fastq/ -o paired_results.csv -p "*_R1.fastq.gz"
```

### 示例4: 高性能多线程处理|Example 4: High-performance Multi-threaded Processing

```bash
biopytools fq-stats -i /data/fastq/ -o results.csv -p "*_1.clean.fq.gz" -t 24
```

## 注意事项|Notes

1. **seqkit安装**: 确保seqkit已安装并在PATH中|Ensure seqkit is installed and in PATH
2. **文件格式**: 支持.fastq/.fq/.fastq.gz/.fq.gz格式|Supports .fastq/.fq/.fastq.gz/.fq.gz formats
3. **内存使用**: 大量文件或大文件可能需要较多内存|Large number of files or large files may require more memory
4. **Excel输出**: 需要安装pandas和openpyxl库|Requires pandas and openpyxl libraries for Excel output

## 故障排除|Troubleshooting

### 错误: 未找到seqkit|Error: seqkit not found

```bash
# 使用conda安装|Install with conda
conda install -c bioconda seqkit

# 或从官网下载|Or download from official website
# https://bioinf.shenwei.me/seqkit/download/
```

### 错误: 需要pandas|Error: pandas required

```bash
pip install pandas openpyxl
```

### 未找到匹配的文件|No matching files found

- 检查输入路径是否正确|Check if input path is correct
- 检查文件匹配模式是否正确|Check if file matching pattern is correct
- 确保文件扩展名被支持|Ensure file extensions are supported

## 版本历史|Version History

- **v1.0.0** (2026-01-20): 初始版本|Initial version
  - 支持单端和双末端FASTQ文件统计|Supports single-end and paired-end FASTQ file statistics
  - 自动识别配对文件|Automatic paired file recognition
  - 多线程处理|Multi-threaded processing
  - CSV/Excel输出|CSV/Excel output
