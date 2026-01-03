# RNA-seq分析流程模块

**完整的RNA-seq数据分析流程|Complete RNA-seq Data Analysis Pipeline**

## 功能概述|Overview

RNA-seq分析流程模块是一个完整的转录组数据分析工具，基于HISAT2和StringTie实现从序列比对到基因表达定量的全流程自动化处理，计算FPKM和TPM值，支持批量样本处理和灵活的参数配置。

## 主要特性|Key Features

- **多种比对算法**: 基于HISAT2的快速准确序列比对
- **灵活输入**: 支持目录自动扫描和文件模式匹配
- **超时保护**: 单样本处理超时自动跳过，避免卡住整个流程
- **表达定量**: StringTie转录本定量和FPKM/TPM计算
- **批量处理**: 自动样本检测和批量并行处理
- **结果整合**: 自动生成表达矩阵和统计报告
- **空间优化**: 可选删除BAM中间文件节省存储空间

## 快速开始|Quick Start

### 基本用法|Basic Usage

```bash
# 基本RNA-seq分析
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i /data/fastq/ -o rnaseq_results

# 指定FASTQ文件命名模式
biopytools rnaseq -g genome.fa -f genes.gtf \\
    --input ./samples/ --output ./analysis/ \\
    --pattern "*.R1.fastq.gz"

# 处理后删除BAM文件节省空间
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i /data/rna_samples/ -o results/ \\
    -t 32 --remove yes
```

### 高级用法|Advanced Usage

```bash
# 自定义超时时间（3小时）
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i /data/rna_samples/ -o results/ \\
    --sample-timeout 10800

# 高性能配置
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i samples/ -o results/ \\
    --threads 64 --pattern "*_1.clean.fq.gz" \\
    --remove yes -v
```

## 参数说明|Parameters

### 必需参数|Required Parameters

| 参数|描述|示例|
|---|---|---|
| `-g, --genome`|基因组FASTA文件路径|`-g genome.fa`|
| `-f, --gtf`|基因注释GTF文件路径|`-f genes.gtf`|
| `-i, --input`|输入FASTQ文件目录或样本信息文件|`-i /data/fastq/`|
| `-o, --output`|输出目录路径|`-o results/`|

### 输入选项|Input Options

| 参数|默认值|描述|
|---|---|---|
| `-p, --pattern`|`None`|FASTQ文件命名模式（如：`*.R1.fastq.gz`或`*_1.fq.gz`），`*`代表样本名|

### 处理参数|Processing Parameters

| 参数|默认值|描述|
|---|---|---|
| `-t, --threads`|`8`|线程数|
| `-r, --remove`|`no`|处理后删除BAM文件（yes/no）|
| `--sample-timeout`|`21600`|单个样本处理超时时间（秒），默认6小时|

### 日志选项|Logging Options

| 参数|默认值|描述|
|---|---|---|
| `-v, --verbose`|`0`|增加输出详细程度（-v: INFO, -vv: DEBUG）|
| `--quiet`|`False`|静默模式，仅输出错误信息|
| `--log-file`|`None`|日志文件路径|
| `--log-level`|`INFO`|日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL）|

### 执行控制|Execution Control

| 参数|默认值|描述|
|---|---|---|
| `--dry-run`|`False`|试运行模式，不实际执行|

## 输出文件|Output Files

### 表达矩阵|Expression Matrix

- **`all.fpkm.tpm.txt`**: 所有样本的FPKM和TPM矩阵
  - 格式：基因ID × 样本的FPKM和TPM值
  - 用途：差异表达分析、基因表达量比较

### HISAT2索引|HISAT2 Index

- **`hisat2_index/`**: HISAT2索引文件目录
  - 位置：基因组文件所在目录
  - 命名：`{genome_name}.hisat2.index.*`

### StringTie输出|StringTie Output

- **`stringtie_output/`**: 每个样本的GTF文件
  - 文件：`sample1.gtf`, `sample2.gtf`, ...
  - 内容：转录本组装和定量结果

### FPKM输出|FPKM Output

- **`fpkm_output/`**: 每个样本的FPKM文件
  - 文件：`sample1.fpkm.txt`, `sample2.fpkm.txt`, ...
  - 格式：基因ID、基因名称、FPKM值

### BAM文件|BAM Files (optional)

- **`*.sorted.bam`**: 比对结果BAM文件
  - 可通过`--remove`参数控制删除
  - 删除后可节省大量存储空间

### 日志文件|Log File

- **`rnaseq_processing_YYYYMMDD_HHMMSS.log`**: 运行日志
  - 包含完整的处理流程和错误信息

## 文件命名模式|File Naming Pattern

### 支持的配对文件格式|Supported Paired-end Formats

| Read1模式|Read2对应|示例|
|---|---|---|
| `*.R1.fastq.gz`|`*.R2.fastq.gz`|`sample1.R1.fastq.gz`, `sample1.R2.fastq.gz`|
| `*_1.fq.gz`|`*_2.fq.gz`|`sample_1.fq.gz`, `sample_2.fq.gz`|
| `*_R1.fq.gz`|`*_R2.fq.gz`|`sample_R1.fq.gz`, `sample_R2.fq.gz`|
| `*_1.fastq.gz`|`*_2.fastq.gz`|`sample_1.fastq.gz`, `sample_2.fastq.gz`|
| `*.1.fq.gz`|`*.2.fq.gz`|`sample.1.fq.gz`, `sample.2.fq.gz`|
| `*_f1.fq.gz`|`*_r2.fq.gz`|`sample_f1.fq.gz`, `sample_r2.fq.gz`|

### 自动模式匹配|Auto Pattern Matching

当不指定`--pattern`参数时，程序会自动尝试以上常见格式。

## 分析流程|Analysis Pipeline

### 1. 构建HISAT2索引|Build HISAT2 Index

从基因组FASTA文件构建HISAT2索引，支持剪接位点和外显子注释。

### 2. 解析输入样本|Parse Input Samples

- 支持目录扫描和文件模式匹配
- 自动配对R1/R2文件
- 验证文件完整性

### 3. 处理所有样本|Process All Samples

对每个样本执行：
- **HISAT2比对**: 序列比对到参考基因组
- **StringTie定量**: 转录本组装和表达定量
- **FPKM值提取**: 提取FPKM表达值
- **BAM文件处理**: 根据参数保留或删除

### 4. 合并表达矩阵|Merge Expression Matrix

- 合并所有样本的FPKM值
- 计算TPM值
- 生成统一的表达矩阵

### 5. 生成总结报告|Generate Summary Report

- 处理样本统计
- 成功/失败样本列表
- 运行时间统计

## 表达量指标|Expression Metrics

### FPKM|Fragments Per Kilobase per Million

- 定义：每百万个片段中，每个千个碱基长度的转录本的片段数
- 用途：校正基因长度和测序深度
- 适用：基因水平比较

### TPM|Transcripts Per Million

- 定义：每百万个转录本中，某个转录本的拷贝数
- 特点：所有样本TPM之和相同
- 适用：跨样本表达量比较

## 使用示例|Examples

### 示例1：标准分析流程|Standard Pipeline

```bash
biopytools rnaseq \\
    -g /data/genome.fa \\
    -f /data/genes.gtf \\
    -i /data/rna_seq/fastq/ \\
    -o /data/rna_seq/results/ \\
    --pattern "*_1.clean.fq.gz" \\
    --threads 32
```

### 示例2：删除中间文件节省空间|Save Disk Space

```bash
biopytools rnaseq \\
    -g genome.fa \\
    -f annotation.gtf \\
    -i fastq_files/ \\
    -o results/ \\
    --remove yes \\
    --threads 64
```

### 示例3：设置样本超时|Set Sample Timeout

```bash
# 3小时超时
biopytools rnaseq \\
    -g genome.fa -f genes.gtf \\
    -i samples/ -o results/ \\
    --sample-timeout 10800

# 12小时超时
biopytools rnaseq \\
    -g genome.fa -f genes.gtf \\
    -i samples/ -o results/ \\
    --sample-timeout 43200
```

### 示例4：详细日志模式|Verbose Logging

```bash
biopytools rnaseq \\
    -g genome.fa -f genes.gtf \\
    -i samples/ -o results/ \\
    -vv --log-file rnaseq.log
```

## 性能建议|Performance Tips

### 线程数优化|Thread Optimization

- HISAT2比对: 建议 16-64 线程
- StringTie定量: 建议 8-16 线程
- 总线程数建议不超过CPU核心数

### 存储空间|Storage Space

- BAM文件占用大量空间（通常为FASTQ的2-3倍）
- 使用`--remove yes`可节省50-70%的存储空间
- 建议使用SSD存储I/O密集文件

### 内存要求|Memory Requirements

- HISAT2索引: ~8-16 GB
- StringTie定量: ~4-8 GB/样本
- 建议至少32 GB可用内存

## 故障排除|Troubleshooting

### 样本处理超时|Sample Processing Timeout

**问题**: 某个样本处理时间过长卡住整个流程

**解决方案**:
```bash
# 使用默认6小时超时
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i samples/ -o results/

# 自定义超时时间
biopytools rnaseq -g genome.fa -f genes.gtf \\
    -i samples/ -o results/ --sample-timeout 10800
```

### 样本配对失败|Sample Pairing Failed

**问题**: R1/R2文件配对不成功

**解决方案**:
- 检查文件命名是否符合支持的模式
- 使用`--pattern`参数指定正确的模式
- 确保R1和R2文件在同一目录

### 内存不足|Insufficient Memory

**问题**: 处理大样本时内存溢出

**解决方案**:
- 减少并行样本数
- 减少线程数
- 增加系统交换空间

## 注意事项|Notes

### 输入要求|Input Requirements

- 需要HISAT2和StringTie已安装并配置在PATH中
- 输入FASTQ文件应为配对末端测序数据
- GTF文件应包含完整的基因注释信息
- 基因组FASTA文件应与GTF注释版本匹配

### 输出说明|Output Description

- 表达矩阵包含所有成功处理的样本
- 失败的样本会被记录在日志中但不会出现在矩阵中
- FPKM和TPM值保留4位小数
- 基因ID来源于GTF文件

### 流程特点|Pipeline Features

- 支持断点续传：已处理的样本不会重复处理
- 自动错误处理：单个样本失败不影响其他样本
- 超时保护：避免单个样本卡住整个流程
- 详细日志：记录每个步骤的详细信息

## 版本历史|Version History

### v0.2.0 (2026-01-03)

- 添加样本处理超时功能
- 优化日志输出格式（中英文对照）
- 改进错误处理和异常捕获
- 添加执行时间统计
- 优化样本解析逻辑

### v0.1.0 (2025-12-20)

- 初始版本发布
- 支持HISAT2 + StringTie标准流程
- 实现FPKM/TPM计算
- 支持批量样本处理
