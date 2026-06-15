# BAM to FASTQ批量转换工具|BAM to FASTQ Batch Conversion Tool

## 功能描述|Function Description

使用bam2fastq将BAM文件批量转换为FASTQ格式的工具|Batch convert BAM files to FASTQ format using bam2fastq

## 功能特点|Features

- 支持批量转换多个BAM文件|Support batch conversion of multiple BAM files
- 支持并行处理提高转换速度|Support parallel processing for faster conversion
- 自动生成R1、R2 FASTQ文件|Automatically generate R1, R2 FASTQ files
- 输出gzip压缩的FASTQ文件|Output gzip compressed FASTQ files
- 详细的转换日志和统计信息|Detailed conversion logs and statistics
- 正确处理软链接文件|Properly handle symlinked files

## 依赖要求|Requirements

- Python 3.7+
- bam2fastq (需要安装|Installation required)

### 安装bam2fastq|Install bam2fastq

```bash
conda install -c bioconda pbmm2 pbbam
```

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 转换单个目录中的所有BAM文件|Convert all BAM files in a directory
biopytools bam2fastq -i ./bam_dir -o ./fastq_dir
```

### 高级用法|Advanced Usage

```bash
# 指定线程数|Specify thread count
biopytools bam2fastq -i ./bam_dir -o ./fastq_dir -t 32

# 并行处理多个BAM文件|Process multiple BAM files in parallel
biopytools bam2fastq -i ./bam_dir -o ./fastq_dir -j 4 -t 16

# 使用自定义bam2fastq路径|Use custom bam2fastq path
biopytools bam2fastq -i ./bam_dir -o ./fastq_dir --bam2fastq-path /path/to/bam2fastq
```

## 参数说明|Parameters

| 短参数|长参数|类型|默认值|说明|Description|
|-------|------|------|--------|------|-----------|
|-i|--input|str|必需|输入文件夹路径(包含BAM文件)|Input directory path (containing BAM files)|
|-o|--output|str|必需|输出文件夹路径|Output directory path|
|-t|--threads|int|64|每个BAM文件转换使用的线程数|Threads per BAM file conversion|
|-j|--jobs|int|1|并行处理的BAM文件数量|Number of parallel BAM file processing|
|--bam2fastq-path|str|'bam2fastq'|bam2fastq可执行文件路径|bam2fastq executable path|

## 输出文件|Output Files

对于每个BAM文件`sample.bam`，bam2fastq会自动生成以下FASTQ文件：

For each BAM file `sample.bam`, bam2fastq will automatically generate the following FASTQ files:

- `sample.fastq.gz` - 主FASTQ文件|Main FASTQ file
- 可能会生成其他后缀的文件（如`_1.fastq.gz`, `_2.fastq.gz`）取决于BAM文件内容|Other suffix files may be generated depending on BAM content

## 使用示例|Examples

### 示例1: 基本转换|Example 1: Basic Conversion

```bash
biopytools bam2fastq -i ./bam_files -o ./fastq_output
```

**输出|Output:**
```
找到|Found 10 个BAM文件|BAM files
每个文件使用|Each file uses 64 个线程|threads
并行处理|Parallel processing 1 个文件|files
正在处理|Processing: sample1.bam
完成|Completed: sample1.bam
...
==================================================
转换完成|Conversion completed!
总数|Total: 10
成功|Success: 10
失败|Failed: 0
成功率|Success rate: 100.00%
==================================================
```

### 示例2: 高性能转换|Example 2: High Performance Conversion

```bash
# 使用4个并行任务，每个任务32个线程
biopytools bam2fastq -i ./bam_files -o ./fastq_output -j 4 -t 32
```

### 示例3: 使用自定义bam2fastq路径|Example 3: Using Custom bam2fastq Path

```bash
biopytools bam2fastq -i ./bam_files -o ./fastq_output \
  --bam2fastq-path /share/org/YZWL/yzwl_lixg/miniforge3/envs/pbbam_v.2.4.0/bin/bam2fastq
```

## 注意事项|Notes

1. **输入目录|Input Directory**: 输入目录应包含`.bam`扩展名的文件|The input directory should contain files with `.bam` extension

2. **线程设置|Thread Settings**:
   - `-t` 控制每个BAM文件转换使用的线程数|`-t` controls threads per BAM file conversion
   - `-j` 控制同时转换的BAM文件数量|`-j` controls number of BAM files converted simultaneously
   - 总线程数 = `-t` × `-j`|Total threads = `-t` × `-j`

3. **输出文件|Output Files**: 所有输出文件都是gzip压缩格式|All output files are in gzip compressed format

4. **日志文件|Log Files**: 转换日志保存在输出目录的`bam2fastq_conversion.log`|Conversion logs are saved in `bam2fastq_conversion.log` in the output directory

5. **软链接处理|Symlink Handling**: 程序会自动解析软链接文件，支持相对路径和绝对路径的软链接|The program automatically resolves symlink files, supporting both relative and absolute path symlinks

## 故障排除|Troubleshooting

### 问题1: 未找到bam2fastq|Issue 1: bam2fastq not found

**错误信息|Error Message:**
```
未找到bam2fastq|bam2fastq not found
```

**解决方法|Solution:**
```bash
conda install -c bioconda pbmm2 pbbam
```

### 问题2: 输入目录为空|Issue 2: Empty input directory

**错误信息|Error Message:**
```
未找到BAM文件|No BAM files found
```

**解决方法|Solution:** 检查输入目录是否包含`.bam`文件|Check if the input directory contains `.bam` files

### 问题3: 软链接解析失败|Issue 3: Symlink resolution failed

**错误信息|Error Message:**
```
软链接解析失败|Symlink resolution failed
```

**解决方法|Solution:** 确保软链接的目标文件存在且可访问|Ensure the symlink target file exists and is accessible

## 技术细节|Technical Details

### bam2fastq命令参数|bam2fastq command parameters

```bash
bam2fastq \
  -o output_prefix \
  -j threads \
  -c 6 \
  input.bam
```

- `-o`: 输出文件前缀|Output file prefix
- `-j`: 线程数|Thread count
- `-c`: gzip压缩级别|Gzip compression level (1-9)

## 版本历史|Version History

| 版本|Version | 日期|Date | 说明|Description |
|------|---------|------|------|
|1.0.0|2026-01-06|初始版本|Initial version |
|2.0.0|2026-01-06|从samtools切换到bam2fastq|Switched from samtools to bam2fastq |
