# 三代转录组比对工具 (Long RNA-seq Alignment)

**版本|Version**: 1.3.0
**作者|Author**: BioPyTools Team
**日期|Date**: 2026-01-06

## 模块简介|Module Introduction

三代转录组比对工具用于将PacBio/Nanopore等长读长转录组数据比对到参考基因组。该模块基于NCBI优化的minimap2流程，专门针对三代转录组数据的特点进行优化，支持长intron检测和可变剪接分析。

Long RNA-seq alignment tool for aligning long-read transcriptome data (PacBio/Nanopore) to reference genome. Based on NCBI-optimized minimap2 pipeline, specifically optimized for third-generation transcriptome data, supporting long intron detection and alternative splicing analysis.

## 功能特点|Features

- **文件夹批量处理|Batch Directory Processing**: 支持输入文件夹，自动批量处理内所有BAM/FASTQ文件
- **多格式输入支持|Multiple Input Format Support**: 支持BAM和FASTQ格式输入，自动检测文件类型
- **智能配对检测|Smart Paired-end Detection**: 自动识别单端/双端测序，自动查找配对的R2文件
- **样本名自动提取|Auto Sample Name Extraction**: 自动从输入文件名提取样本名称
- **优化比对参数|Optimized Alignment Parameters**: 采用NCBI推荐的minimap2参数，专门针对三代转录组数据优化
- **BAM转FASTQ|BAM to FASTQ**: 自动将BAM格式转换为FASTQ格式
- **剪接比对|Splice Alignment**: 使用`-ax splice`模式，支持长intron检测（默认100kb）
- **质量控制|Quality Control**: 自动进行比对质量检查和统计
- **并行处理|Parallel Processing**: 支持多线程，提高处理速度
- **详细统计|Detailed Statistics**: 生成flagstat、stats和覆盖度统计

## 安装依赖|Installation Dependencies

```bash
# conda安装
conda install -c bioconda minimap2 samtools

# 或使用bioconda
conda install -c bioconda/minimap2 minimap2
conda install -c bioconda samtools
```

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 单个文件模式
biopytools longrnaseq -i input.bam -r genome.fa -o output_dir

# 文件夹批量模式
biopytools longrnaseq -i /path/to/bam_folder -r genome.fa -o output_dir

# 或使用Python模块
python -m biopytools.longrnaseq -i input.bam -r genome.fa -o output_dir
```

### 参数说明|Parameters

| 参数|Parameter | 简写|Short | 必需|Required | 说明|Description | 默认值|Default |
|---|---|---|---|---|---|---|
| `--input-file` | 输入文件或文件夹 | `-i` | ✓ | Input file or directory (BAM/FASTQ) | - |
| `--ref-genome` | 参考基因组文件 | `-r` | ✓ | Reference genome FASTA file | - |
| `--output-dir` | 输出目录 | `-o` | ✓ | Output directory | - |
| `--sample-name` | 样本名称 | `-s` | ✗ | Sample name (auto-extracted, folder name for directory mode) | - |
| `--threads` | 线程数 | `-t` | ✗ | Number of threads | 64 |
| `--max-intron` | 最大intron长度(bp) | - | ✗ | Maximum intron length | 100000 |
| `--min-mapq` | 最小mapping quality | - | ✗ | Minimum mapping quality (0-60) | 20 |
| `--no-secondary` | 不输出次优比对 | - | ✗ | Do not output secondary alignments | False |
| `--minimap2-path` | minimap2路径 | - | ✗ | minimap2 executable path | "minimap2" |
| `--samtools-path` | samtools路径 | - | ✗ | samtools executable path | "samtools" |

### 高级用法|Advanced Usage

```bash
# 文件夹批量处理
biopytools longrnaseq \
  -i /path/to/bam_folder \
  -r genome.fa \
  -o output_dir \
  -t 64

# 自定义比对参数
biopytools longrnaseq \
  -i input.bam \
  -r genome.fa \
  -o output_dir \
  -t 64 \
  --max-intron 50000 \
  --min-mapq 30

# 指定样本名称
biopytools longrnaseq \
  -i input.bam \
  -r genome.fa \
  -o output_dir \
  -s custom_name

# 仅保留最优比对（不输出次优比对）
biopytools longrnaseq \
  -i input.bam \
  -r genome.fa \
  -o output_dir \
  --no-secondary

# 指定工具路径
biopytools longrnaseq \
  -i input.bam \
  -r genome.fa \
  -o output_dir \
  --minimap2-path /path/to/minimap2 \
  --samtools-path /path/to/samtools
```

## 输出文件|Output Files

### 目录结构|Directory Structure

```
output_dir/
├── alignment_results/
│   ├── logs/
│   │   └── alignment_sample1.log          # 日志文件|Log file
│   ├── stats/
│   │   ├── sample1_stats.txt              # 基础统计|Basic statistics
│   │   ├── sample1_detail_stats.txt       # 详细统计|Detailed statistics
│   │   └── sample1_summary.txt            # 样本汇总|Sample summary
│   ├── tmp/
│   │   └── sample1.fastq                  # 临时FASTQ|Temporary FASTQ
│   ├── sample1.sorted.bam                # 最终BAM文件|Final BAM file
│   ├── sample1.sorted.bam.bai            # BAM索引|BAM index
│   └── alignment_summary.txt             # 汇总报告|Summary report
```

### 主要输出文件说明|Main Output Files

#### 1. BAM文件|BAM File
- **文件名|Filename**: `{sample_name}.sorted.bam`
- **说明|Description**: 排序后的比对结果BAM文件，包含所有比对上的reads
- **用途|Usage**: 可直接用于后续分析，如转录本组装、表达定量等

#### 2. 统计文件|Statistics Files

**sample1_stats.txt** - flagstat统计
```
在总|in total: 1000000
比对上|mapped (99.85%): 998500
```

**sample1_detail_stats.txt** - 详细统计
包含samtools stats的所有统计信息

**alignment_summary.txt** - 汇总报告
包含分析参数、运行时间、样本信息等汇总

## 流程说明|Workflow Description

### 分析步骤|Analysis Steps

1. **准备FASTQ文件|Prepare FASTQ File**
   - 如果输入是BAM：使用samtools fastq将BAM转换为FASTQ格式
   - 如果输入是FASTQ：直接使用
   - 自动检测单端或双端测序
   - 统计reads数量

2. **minimap2比对|Align with minimap2**
   - 使用`-ax splice`模式进行剪接比对
   - 参数优化：
     - `-uf`: 全长转录本模式
     - `--eqx`: 使用=和X表示匹配/错配
     - `-G {max_intron}`: 设置最大intron长度
     - `-p 0.9`: 最小次优/最优比对分数比
     - `-N 25`: 最多报告25个次优比对
   - 管道处理：minimap2 → samtools view → samtools sort

3. **建立索引|Build Index**
   - 使用samtools index建立BAM索引

4. **生成统计|Generate Statistics**
   - flagstat: 基础比对统计
   - stats: 详细比对统计
   - depth: 覆盖度统计和平均深度计算

5. **质量控制|Quality Control**
   - 检查BAM文件完整性
   - 检查比对率（<70%会警告）

### minimap2参数说明|minimap2 Parameters

```bash
-ax splice           # 剪接比对模式|Splice alignment mode
-uf                  # 全长转录本模式|Full-length transcript mode
--eqx                # 使用=和X表示匹配/错配|Use = and X for match/mismatch
--MD                 # 生成MD标签|Generate MD tag
--cs                 # 生成CIGAR字符串|Generate CIGAR string
-Y                   # 使用软剪切|Use soft clipping
--sam-hit-only       # 只输出有比对的reads|Only output aligned reads
-G {max_intron}      # 最大intron长度|Maximum intron length
-p 0.9               # 最小次优/最优比对分数比|Min secondary/optimal score ratio
-N 25                # 最多报告25个次优比对|Report up to 25 secondary alignments
-t {threads}         # 线程数|Number of threads
```

## 示例|Examples

### 示例1：单个文件分析|Example 1: Single File Analysis

```bash
biopytools longrnaseq \
  -i jingyeA.sreads.bam \
  -r OV53_1.hifihic.chr.fa \
  -o alignment_results
```

### 示例2：文件夹批量处理|Example 2: Batch Directory Processing

```bash
biopytools longrnaseq \
  -i /path/to/bam_folder \
  -r genome.fa \
  -o alignment_results \
  -t 64
```

**说明|Notes**:
- 自动识别文件夹内所有 `.bam`、`.fastq`、`.fq.gz` 等文件
- 自动排除R2文件（只处理R1）
- 每个文件自动提取样本名进行比对
- 输出汇总报告显示成功/失败文件列表

### 示例3：自定义参数|Example 3: Custom Parameters

```bash
biopytools longrnaseq \
  -i pacbio.bam \
  -r genome.fa \
  -o results \
  -t 64 \
  --max-intron 50000 \
  --min-mapq 30
```

### 示例4：指定样本名称|Example 4: Specify Sample Name

```bash
biopytools longrnaseq \
  -i input.bam \
  -r genome.fa \
  -o results \
  -s custom_name
```

## 后续分析建议|Downstream Analysis Recommendations

1. **转录本组装|Transcript Assembly**
   ```bash
   stringtie alignments.bam -G genes.gtf -o transcripts.gtf
   ```

2. **基因表达定量|Gene Expression Quantification**
   ```bash
   featureCounts -a genes.gtf -o counts.txt alignments.bam
   ```

3. **可变剪接分析|Alternative Splicing Analysis**
   - rMATS: 检测可变剪接事件
   - SUPPA2: 快速可变剪接分析

4. **融合基因检测|Fusion Gene Detection**
   - STAR-Fusion: 融合基因检测
   - FusionCatcher: 融合转录本检测

5. **长reads特异分析|Long-read Specific Analysis**
   - TALON: 长reads转录本注释
   - SQANTI3: 长reads转录本质量控制和分类

## 常见问题|FAQ

### Q: 如何设置最大intron长度？
**A**: 使用`--max-intron`参数。对于大多数植物，可以使用100000（100kb）；对于动物，通常使用200000（200kb）。

### Q: 比对率较低怎么办？
**A**:
1. 检查参考基因组是否正确
2. 确认测序数据质量
3. 尝试调整`--max-intron`参数
4. 降低`--min-mapq`阈值

### Q: 如何只保留最优比对？
**A**: 使用`--no-secondary`参数。

### Q: 支持哪些输入格式？
**A**: 支持BAM和FASTQ格式（.fq, .fastq, .fq.gz, .fastq.gz）。FASTQ会自动检测单端或双端测序，双端测序会自动查找配对的R2文件。也支持输入文件夹进行批量处理。

### Q: 如何批量处理文件夹中的多个BAM文件？
**A**: 直接指定文件夹路径即可：
```bash
biopytools longrnaseq -i /path/to/bam_folder -r genome.fa -o output_dir
```
工具会自动识别文件夹内所有BAM和FASTQ文件并依次处理。

### Q: 文件夹批量处理时如何命名样本？
**A**: 每个文件会自动从其文件名提取样本名称。自动去除以下标记：
- **测序标记**: `_R1`, `_R2`, `_1`, `_2`, `-R1`, `-R2` 等
- **中间标记**: `.clean`, `.trimmed`, `.filtered` 等
- **扩展名**: `.bam`, `.fastq`, `.fq`, `.gz` 等

示例：
- `sample_1.fq.gz` → `sample`
- `sample_1.clean.fq.gz` → `sample`
- `sample_1.fastq` → `sample`
- `sample.fq` → `sample`
- `sample.fq.gz` → `sample`
- `sample_R1.fastq.gz` → `sample`

## 性能优化建议|Performance Optimization

1. **线程数|Threads**: 建议设置为CPU核心数的80-90%
2. **内存|Memory**: samtools sort使用`-m 4G`，可根据可用内存调整
3. **存储|Storage**: 确保有足够的临时空间（约为输入文件的2-3倍）

## 更新日志|Changelog

### Version 1.3.0 (2026-01-06)
- 支持文件夹批量处理|Support batch directory processing
- `-i`参数可接受文件或文件夹路径|`-i` parameter accepts file or directory path
- 自动识别并处理文件夹内所有BAM/FASTQ文件|Auto-detect and process all BAM/FASTQ files in directory
- 自动排除R2文件（只处理R1）|Auto-exclude R2 files (process only R1)
- 批量处理汇总报告|Batch processing summary report

### Version 1.2.0 (2026-01-06)
- 样本名称自动提取|Auto-extract sample name from input filename
- `-s`参数改为可选|`-s` parameter changed to optional

### Version 1.1.0 (2026-01-06)
- 新增FASTQ格式输入支持|Add FASTQ format input support
- 自动检测单端/双端测序|Auto-detect single-end/paired-end sequencing
- 自动查找配对R2文件|Auto-find paired R2 file
- 简化输出文件名（去除_aligned后缀）|Simplify output filename (remove _aligned suffix)
- 参数顺序优化（-i参数移至第一位）|Optimize parameter order (-i moved to first position)
- 默认线程数调整为64|Default threads changed to 64

### Version 1.0.0 (2026-01-06)
- 初始版本发布|Initial release
- 支持三代转录组比对|Support long RNA-seq alignment
- 优化的minimap2参数|Optimized minimap2 parameters
- 完整的统计和质量控制|Complete statistics and quality control

## 参考资源|References

- [minimap2文档](https://lh3.github.io/minimap2/minimap2.html)
- [samtools文档](http://www.htslib.org/doc/samtools.html)
- [NCBI转录组比对流程](https://github.com/ncbi/graft)
