# BWA 全基因组比对工具

**Burrows-Wheeler Aligner 全基因组比对分析工具 | BWA Whole Genome Alignment Tool**

## 📖 功能概述 | Overview

BWA (Burrows-Wheeler Aligner) 是广泛应用的基因组序列比对工具，特别适合将大规模测序数据（short reads）比对到参考基因组。本工具基于BWA-MEM算法封装，提供从基因组索引构建、序列比对、BAM文件处理到覆盖度分析的完整流程，支持批量样本处理和结果汇总。

## ✨ 主要特性 | Key Features

- **🎯 BWA-MEM算法**: 基于BWA-MEM算法，适合70bp-1Mbp的reads比对
- **📊 批量处理**: 支持目录批量分析，自动识别PE reads配对
- **🧬 完整流程**: 索引构建 → 比对 → 排序 → 去重 → 覆盖度分析
- **📈 统计汇总**: 自动生成比对统计和覆盖度报告
- **🔄 断点续传**: 支持跳过已完成样本的续传功能
- **⚙️ 灵活配置**: 支持BWA算法参数、打分参数、覆盖度参数的全面配置
- **📏 滑窗分析**: 支持基于滑动窗口的覆盖度分布分析
- **🚀 高性能**: 多线程并行处理，充分利用计算资源

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 批量比对分析（PE reads）
biopytools bwa \
    -g genome.fa \
    -i ./fastq_dir/ \
    -p "_1.clean.fq.gz" \
    -t 24 \
    -o ./bwa_results
```

### 高级用法 | Advanced Usage

```bash
# 自定义BWA参数 + 去重 + 覆盖度分析
biopytools bwa \
    -g genome.fa \
    -i ./fastq_dir/ \
    -p "_R1.fq.gz" \
    --markdup \
    --window-size 1000000 \
    --step-size 100000 \
    -t 32 \
    -o ./results
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 参考基因组FASTA文件 | `-g genome.fa` |
| `-i, --input` | 输入FASTQ文件目录 | `-i ./fastq/` |
| `-p, --pattern` | FASTQ文件匹配模式（reads 1） | `-p "_1.clean.fq.gz"` |

**文件命名规则**：
- PE reads文件名必须成对，如：
  - `sampleA_1.clean.fq.gz` 和 `sampleA_2.clean.fq.gz`
  - `sampleB_R1.fq.gz` 和 `sampleB_R2.fq.gz`
- `-p`参数指定reads 1的识别模式，reads 2通过自动替换识别
  - `_1` → `_2`
  - `_R1` → `_R2`

### 输出参数 | Output Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./bwa_output` | 📁 输出目录路径 |

### 性能参数 | Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `88` | ⚙️ 线程数（建议≤CPU核心数） |

### BWA算法参数 | BWA Algorithm Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--bwa-k` | `19` | 🎯 最小种子长度 |
| `--bwa-w` | `100` | 📏 带宽 |
| `--bwa-d` | `100` | 📐 对角线X-dropoff值 |
| `--bwa-r` | `1.5` | 🔢 内部种子因子 |
| `--bwa-c` | `500` | 🔍 种子出现次数阈值 |
| `--bwa-D` | `0.50` | ⏬ 短链丢弃比例 |
| `--bwa-W` | `0` | 🔗 最小链长 |
| `--bwa-m` | `50` | 🔄 配对拯救轮数 |
| `--bwa-S` | `False` | ⏭️ 跳过配对拯救 |
| `--bwa-P` | `False` | ⏭️ 跳过配对 |

### BWA打分参数 | BWA Scoring Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--bwa-A` | `1` | ✅ 匹配得分 |
| `--bwa-B` | `4` | ❌ 错配罚分 |
| `--bwa-O` | `"6,6"` | 🔓 Gap开放罚分（read1,read2） |
| `--bwa-E` | `"1,1"` | ➡️ Gap延伸罚分（read1,read2） |
| `--bwa-L` | `"5,5"` | ✂️ 末端剪切罚分 |
| `--bwa-U` | `17` | ⚠️ 未配对罚分 |

### BWA输出参数 | BWA Output Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--bwa-M` | `False` | 🏷️ 标记次要比对（适合Picard） |
| `--bwa-T` | `30` | 🎯 最小输出得分 |
| `--bwa-a` | `False` | 📋 输出所有比对 |
| `--bwa-C` | `False` | 📝 附加FASTQ注释 |
| `--bwa-V` | `False` | 🧬 输出参考序列头 |
| `--bwa-Y` | `False` | 🔄 软剪切补充比对 |

### 后处理参数 | Post-processing Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--markdup` | `False` | 🏷️ 标记重复序列（使用samtools markdup） |
| `--remove-dup` | `False` | 🗑️ 移除重复序列 |

### 覆盖度参数 | Coverage Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-base-quality` | `0` | 🎨 最小碱基质量（samtools depth -q） |
| `--min-mapping-quality` | `0` | 🧭 最小比对质量（samtools depth -Q） |
| `--max-depth` | `0` | 📊 最大深度限制（0=无限制） |

### 滑窗参数 | Window Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--window-size` | `1000000` | 📏 窗口大小（bp） |
| `--step-size` | `100000` | 👣 步长（bp） |

### 其他参数 | Other Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--resume` | `False` | 🔄 断点续传（跳过已完成样本） |
| `--keep-sam` | `False` | 💾 保留SAM文件 |

## 💡 使用示例 | Usage Examples

### 示例1：基本比对流程 | Example 1: Basic Alignment Pipeline

```bash
# 使用默认参数进行批量比对
biopytools bwa \
    -g reference.fa \
    -i ./clean_data/ \
    -p "_1.clean.fq.gz" \
    -t 24 \
    -o ./alignment_output
```

### 示例2：标记重复序列 | Example 2: Mark Duplicates

```bash
# 比对后标记重复序列
biopytools bwa \
    -g reference.fa \
    -i ./fastq/ \
    -p "_R1.fq.gz" \
    --markdup \
    -t 32 \
    -o ./bwa_markdup
```

### 示例3：自定义BWA参数 | Example 3: Custom BWA Parameters

```bash
# 调整种子长度和打分系统，适合高变异物种
biopytools bwa \
    -g reference.fa \
    -i ./reads/ \
    -p "_1.fq.gz" \
    --bwa-k 15 \
    --bwa-A 2 \
    --bwa-B 3 \
    --bwa-O "5,5" \
    -t 24 \
    -o ./custom_params
```

### 示例4：精细覆盖度分析 | Example 4: Fine-grained Coverage Analysis

```bash
# 使用小窗口进行覆盖度分布分析
biopytools bwa \
    -g reference.fa \
    -i ./data/ \
    -p "_1.fq.gz" \
    --window-size 500000 \
    --step-size 50000 \
    --min-base-quality 20 \
    --min-mapping-quality 20 \
    -t 24 \
    -o ./coverage_analysis
```

### 示例5：断点续传 | Example 5: Resume Processing

```bash
# 批量处理中断后继续
biopytools bwa \
    -g reference.fa \
    -i ./large_dataset/ \
    -p "_R1.fq.gz" \
    --resume \
    -t 48 \
    -o ./batch_alignment
```

### 示例6：最小参数模式 | Example 6: Minimal Output Mode

```bash
# 使用较宽松的输出阈值
biopytools bwa \
    -g reference.fa \
    -i ./reads/ \
    -p "_1.fq.gz" \
    --bwa-T 20 \
    --bwa-a \
    -t 16 \
    -o ./sensitive_align
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
bwa_output/
├── bam/                       # BAM文件目录|BAM files directory
│   ├── sample1.bam           # 最终BAM文件（已排序、索引）|Final BAM (sorted, indexed)
│   ├── sample1.bam.bai       # BAM索引文件|BAM index file
│   └── sample2.bam
├── coverage/                   # 覆盖度文件目录|Coverage files directory
│   ├── sample1.depth.txt     # 覆盖度深度文件|Coverage depth file
│   └── sample2.depth.txt
├── windows/                    # 滑窗统计目录|Window statistics directory
│   ├── sample1.windows.txt   # 窗口覆盖度统计|Window coverage stats
│   └── sample2.windows.txt
├── stats/                      # 统计文件目录|Statistics directory
│   ├── sample1.flagstat.txt  # 比对统计|Alignment statistics
│   ├── sample1.idxstats.txt  # 索引统计|Index statistics
│   └── alignment_summary.txt # 汇总报告|Summary report
└── logs/                       # 日志目录|Log directory
    └── bwa_alignment.log      # 运行日志|Run log
```

### 输出文件说明 | Output File Description

#### 1. BAM文件
- **`*.bam`**: 排序并索引的最终BAM文件，可直接用于下游分析
- **`*.bam.bai`**: BAM索引文件，用于快速随机访问

#### 2. 覆盖度文件
- **`*.depth.txt`**: 碱基级别的覆盖度深度文件
  ```
  chr1    1       25
  chr1    2       26
  chr1    3       24
  ...
  ```
  - 列1: 染色体
  - 列2: 位置
  - 列3: 覆盖深度

#### 3. 滑窗统计文件
- **`*.windows.txt`**: 滑动窗口覆盖度统计
  ```
  chr1    1       1000000     25.3
  chr1    100001  2000000     26.1
  ...
  ```
  - 列1: 染色体
  - 列2: 窗口起始
  - 列3: 窗口结束
  - 列4: 平均覆盖度

#### 4. 统计文件
- **`*.flagstat.txt`**: SAMtools flagstat统计结果
  ```
  1245032 + 0 in total ( QC-passed reads + QC-failed reads )
  1200456 + 0 mapped ( 96.39% : N/A )
  ...
  ```

- **`*.idxstats.txt`**: BAM索引统计
  ```
  chr1    230208000       11523456        23
  chr2    190000000       9876543         21
  ...
  ```
  - 列1: 参考序列名
  - 列2: 序列长度
  - 列3: 比对reads数
  - 列4: 未比对reads数

- **`alignment_summary.txt`**: 所有样本的汇总统计表
  ```
  Sample      Total_Reads    Mapped_Reads    Map_Rate    Cov_1X    Cov_5X    Cov_10X
  sample1     1245032        1200456         96.39%      95.2%     92.1%     88.5%
  ...
  ```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **BWA** (版本 0.7.x 或更新)
  - 安装: `conda install -c bioconda bwa`
  - 或: `sudo apt-get install bwa`

- **Samtools** (版本 1.10 或更新)
  - 安装: `conda install -c bioconda samtools`
  - 或: `sudo apt-get install samtools`

- **Python依赖**:
  - Python >= 3.7
  - pandas (统计汇总)

### 硬件要求 | Hardware Requirements

| 数据规模 | 内存 | 磁盘空间 | CPU |
|----------|------|----------|-----|
| 小规模 (<10 Gb) | 8-16 GB | 50 GB | 8-16 核 |
| 中规模 (10-50 Gb) | 32-64 GB | 200 GB | 16-32 核 |
| 大规模 (>50 Gb) | 128+ GB | 500+ GB | 32+ 核 |

### 环境配置 | Environment Setup

```bash
# 创建Conda环境
conda create -n bwa_env python=3.9
conda activate bwa_env

# 安装BWA和Samtools
conda install -c bioconda bwa samtools

# 验证安装
bwa
samtools --version
```

## ⚠️ 注意事项 | Important Notes

1. **文件命名**: PE reads文件必须严格配对命名，工具通过替换识别reads 2
2. **内存需求**: BWA索引需要加载到内存，大基因组建议至少64GB RAM
3. **磁盘空间**: 中间文件（SAM）很大，确保有足够临时空间
4. **线程设置**: 建议设置为CPU核心数的50-80%，避免系统过载
5. **基因组索引**: 首次运行会自动构建索引，耗时较长
6. **覆盖度分析**: 大基因组的深度覆盖度分析会消耗大量时间和磁盘
7. **样本续传**: 使用`--resume`时需确保之前运行的输出目录完整

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "bwa: command not found" 错误**

```bash
# 安装BWA
conda install -c bioconda bwa

# 验证安装
which bwa
bwa
```

**Q: 找不到配对的reads文件**

```bash
# 检查文件命名
ls -1 ./fastq/

# 确保reads文件成对:
# sampleA_1.fq.gz  <-- pattern匹配
# sampleA_2.fq.gz  <-- 自动识别

# 或:
# sampleA_R1.fq.gz  <-- pattern匹配
# sampleA_R2.fq.gz  <-- 自动识别
```

**Q: 内存不足错误**

```bash
# 减少线程数以降低内存使用
biopytools bwa -g genome.fa -i ./fastq/ -p "_1.fq.gz" -t 8 -o ./output
```

**Q: SAM文件过大导致磁盘空间不足**

```bash
# 管道式处理避免中间SAM文件（需修改配置或手动处理）
# 或使用更大的磁盘空间
```

**Q: 覆盖度统计耗时过长**

```bash
# 使用窗口统计代替碱基级别统计
# 或限制最大深度
biopytools bwa ... --window-size 1000000 --max-depth 1000
```

**Q: 某些样本比对失败**

```bash
# 检查日志文件确定失败原因
cat bwa_output/logs/bwa_alignment.log | grep -i error

# 单独重新运行失败的样本
# （需要手动处理或使用resume功能）
```

## 📚 相关资源 | Related Resources

- [BWA官方文档](http://bio-bwa.sourceforge.net/)
- [BWA-MEM论文](https://arxiv.org/abs/1303.3997)
- [Samtools文档](http://www.htslib.org/)
- [SAM格式规范](https://samtools.github.io/hts-specs/)
- [BWA使用教程](https://bio-bwa.sourceforge.net/bwa.shtml#3)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

BWA和Samtools软件遵循其原始许可证。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用BWA工具，请引用原始文献：

```
Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.
arXiv preprint arXiv:1303.3997. 2013.

Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform.
Bioinformatics, 2009, 25(14): 1754-1760.
doi: 10.1093/bioinformatics/btp324
```

如果在学术研究中使用Samtools，请引用：

```
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al.
The Sequence Alignment/Map format and SAMtools.
Bioinformatics, 2009, 25(16): 2078-2079.
doi: 10.1093/bioinformatics/btp352
```
