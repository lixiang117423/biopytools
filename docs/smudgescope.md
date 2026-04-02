# SmudgeScope基因组分析工具 | SmudgeScope Genome Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2026-04-01

## 概述 | Overview

SmudgeScope工具整合了**Jellyfish**、**GenomeScope 2.0**和**Smudgeplot**，提供完整的基因组k-mer分析流程，用于：
- 基因组大小估计
- 杂合度和倍性分析
- 重复序列含量评估
- 古老多倍体识别

SmudgeScope integrates **Jellyfish**, **GenomeScope 2.0**, and **Smudgeplot** to provide a complete genomic k-mer analysis pipeline for:
- Genome size estimation
- Heterozygosity and ploidy analysis
- Repeat content assessment
- Ancient polyploid identification

## 功能特点 | Features

- 🧬 **完整流程**: Jellyfish计数 → GenomeScope分析 → Smudgeplot倍性推断一站式完成
- 📊 **多维度分析**: k-mer分布、基因组大小、杂合度、倍性全面评估
- 🎯 **智能推断**: 自动推断倍性水平（1-6倍体）
- 💾 **断点续传**: 支持中断恢复，自动跳过已完成步骤
- ⚡ **高性能**: 支持多线程并行计算，优化内存使用
- 📈 **结果汇总**: 自动汇总所有样品的分析结果到all/目录
- 🔧 **灵活配置**: 丰富的参数配置选项，适应不同数据类型

## 分析流程 | Analysis Pipeline

### 步骤1: Jellyfish Count | K-mer计数

使用Jellyfish进行k-mer计数，生成k-mer频率分布。

Performs k-mer counting using Jellyfish, generates k-mer frequency distribution.

**输出文件 | Output Files:**
- `*.jf` - Jellyfish二进制格式k-mer计数文件
- `*.jellyfish.histo` - k-mer频率直方图

### 步骤2: GenomeScope 2.0 | 基因组特征分析

基于k-mer分布分析基因组特征，估计基因组大小、杂合度、重复率等。

Analyzes genomic characteristics based on k-mer distribution, estimates genome size, heterozygosity, repeat rate, etc.

**输出文件 | Output Files:**
- `model.txt` - GenomeScope拟合模型参数
- `linear_plot.png` - 线性化k-mer分布图
- `transformed_linear_plot.png` - 变换后的线性图
- `log_plot.png` - 对数尺度k-mer分布图
- `transformed_log_plot.png` - 变换后的对数图
- `summary.txt` - 分析结果摘要

**关键参数 | Key Parameters:**
- `kcov` - k-mer coverage，用于后续Smudgeplot分析
- `genome_size` - 估计的基因组大小
- `het` - 杂合度比例
- `repeat` - 重复序列比例

### 步骤3: FastK | 生成FastK表（Smudgeplot输入）

运行FastK生成Smudgeplot所需的k-mer对表格。

Runs FastK to generate k-mer pair table required for Smudgeplot.

**输出文件 | Output Files:**
- `fastk_table.ktab` - FastK二进制表格
- `fastk_table.hist` - k-mer频率直方图

### 步骤4: Smudgeplot Hetmers | 异质k-mer对分析

识别异质k-mer对，用于推断倍性和基因组结构。

Identifies heterozygous k-mer pairs for ploidy and genomic structure inference.

**输出文件 | Output Files:**
- `*.kmerpairs.smu` - k-mer对数据
- `*.sma` - 辅助文件

### 步骤5: Smudgeplot Plot | 倍性分析可视化

生成Smudgeplot图形，直观展示基因组倍性特征。

Generates Smudgeplot plots to visualize genomic ploidy characteristics.

**输出文件 | Output Files:**
- `*_smudgeplot.png` - 线性尺度Smudgeplot图
- `*_smudgeplot_log10.png` - 对数尺度Smudgeplot图
- `*.smudge_report.tsv` - 分析报告

## 安装和使用 | Installation and Usage

### 前置要求 | Prerequisites

#### 软件依赖 | Software Dependencies

```bash
# Jellyfish (k-mer计数)
conda install -c bioconda jellyfish

# GenomeScope 2.0 (基因组特征分析)
conda create -n genomescope_v.2.0.1 -c conda-forge genomescope2 r-base r-data.table

# Smudgeplot (倍性分析)
conda create -n smudgeplot -c conda-forge smudgeplot r-smudgeplot r-ggplot2 r-gridextra r-cowplot

# FastK (Smudgeplot的C后端)
# 随smudgeplot自动安装
```

#### Python环境要求 | Python Environment Requirements

```bash
# biopytools会自动处理Python依赖
# 确保biopytools已安装
pip install biopytools
```

### 基本用法 | Basic Usage

#### 单个样品分析 | Single Sample Analysis

```bash
# 基本用法
biopytools smudgescope -i clean/R0040_1_1.clean.fq.gz -o output/R0040_1

# 指定读长和k-mer大小
biopytools smudgescope -i clean/R0040_1_1.clean.fq.gz -o output/R0040_1 -l 150 -k 21

# 跳过Smudgeplot分析（仅运行GenomeScope）
biopytools smudgescope -i clean/R0040_1_1.clean.fq.gz -o output/R0040_1 --skip-smudgeplot
```

#### 批量分析 | Batch Analysis

```bash
# 指定输入目录和read1文件后缀模式
biopytools smudgescope -i clean/ -o output/ --read1-suffix "*_1.clean.fq.gz"

# 指定GenomeScope conda环境
biopytools smudgescope -i clean/ -o output/ --genomescope-env genomescope_v.2.0.1

# 指定倍性（不使用Smudgeplot自动推断）
biopytools smudgescope -i clean/ -o output/ --ploidy 2
```

### 参数说明 | Parameter Description

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 输入FASTQ文件或目录 | 必需 |
| `-o, --output-dir` | 输出目录 | 必需 |
| `-l, --read-length` | 测序读长 | 150 |
| `-k, --kmer-size` | K-mer大小 | 21 |
| `-t, --threads` | 线程数 | 12 |
| `-s, --hash-size` | Jellyfish哈希表大小 | 10G |
| `-c, --max-kmer-cov` | 最大k-mer覆盖度 | 1000 |
| `--skip-smudgeplot` | 跳过Smudgeplot倍性分析 | False |
| `--ploidy` | 基因组倍性(1-6) | 2（或由Smudgeplot推断） |
| `--genomescope-env` | GenomeScope conda环境名 | genomescope_v.2.0.1 |
| `--read1-suffix` | Read1文件后缀模式 | *_1.clean.fq.gz |

## 输出文件结构 | Output File Structure

```
output/
├── R0040_1/                          # 单个样品目录
│   ├── 00_pipeline_info/             # 流程信息
│   │   └── software_versions.yml
│   ├── 01_jellyfish/                # Jellyfish输出
│   │   ├── R0040_1.jellyfish.jf
│   │   └── R0040_1.jellyfish.histo
│   ├── 02_genomescope/              # GenomeScope输出
│   │   ├── model.txt
│   │   ├── linear_plot.png
│   │   ├── log_plot.png
│   │   └── summary.txt
│   ├── fastk/                       # FastK输出
│   │   ├── fastk_table.ktab
│   │   └── fastk_table.hist
│   └── 03_smudgeplot/               # Smudgeplot输出
│       ├── R0040_1.kmerpairs.smu
│       ├── R0040_1_smudgeplot.png
│       └── R0040_1.smudge_report.tsv
└── all/                              # 汇总目录
    ├── R0040_1_genomescope.png      # GenomeScope图形汇总
    ├── R0040_1_smudgeplot.png       # Smudgeplot图形汇总
    └── ...
```

## 结果解读 | Result Interpretation

### GenomeScope结果 | GenomeScope Results

**model.txt文件包含关键参数：**
- `len` - 估计的基因组大小
- `het` - 杂合度比例（0-1）
- `rep` - 重复序列比例
- `kcov` - k-mer coverage

**线性图解读：**
- 横轴：k-mer深度（coverage）
- 纵轴：k-mer数量（对数尺度）
- 不同颜色的线代表不同的基因组组分

### Smudgeplot结果 | Smudgeplot Results

**倍性判断：**
- **二倍体(2x)**: 左侧有明显的1:1异质k-mer对分布
- **三倍体(3x)**: 左侧有明显的AA:AB:BB=1:2:1分布
- **四倍体(4x)**: 左侧有明显的AA:AB:BB:CCC=1:2:1:...分布
- **古多倍体**: 多个分布峰，显示基因组杂交历史

**图形解读：**
- 左侧区域：低频异质k-mer对
- 中间区域：同质k-mer对
- 右侧区域：高拷贝数k-mer对

## 常见问题 | FAQ

### 1. Smudgeplot失败怎么办？

**问题**：Smudgeplot步骤失败，提示"可能是纯合二倍体"。

**原因**：
- 纯合二倍体基因组杂合度极低
- k-mer coverage过低
- 阈值设置过高

**解决方案**：
- 这是正常现象，纯合基因组无法进行Smudgeplot分析
- GenomeScope结果仍然有效，可以用于基因组大小估计
- 检查k-mer coverage是否合理（建议>30）

### 2. 内存不足错误

**问题**：Jellyfish或FastK报错内存不足。

**解决方案**：
```bash
# 减小哈希表大小
biopytools smudgescope -i input/ -o output/ -s 5G

# 减少线程数（降低内存使用）
biopytools smudgescope -i input/ -o output/ -t 32
```

### 3. 如何选择k-mer大小？

**推荐设置：**
- **短读长(150bp)**: k=21
- **长读长(>10kb)**: k=31
- **超长读长(>50kb)**: k=51

**原则**：
- k值越大，特异性越强，但k-mer数量越少
- k值应小于读长的1/2
- 常用值：21、31、51

### 4. 如何解读倍性结果？

**Smudgeplot倍性推断：**
- **自动推断**：基于k-mer对分布模式自动判断
- **手动指定**：如果已知倍性，可以使用`--ploidy`参数

**注意事项：**
- 杂合度<1%的基因组可能无法准确推断倍性
- 古老多倍体可能显示复杂的k-mer对分布模式
- 结合其他证据（如染色体数目）综合判断

### 5. 为什么GenomeScope返回非0退出码？

**现象**：GenomeScope成功完成但returncode=1，错误信息显示"Loading custom .Rprofile..."。

**原因**：
- conda run加载了用户home目录的.Rprofile
- R启动信息被误认为错误

**解决方案**：
- 这是已修复的问题（v1.0.0+）
- 如果仍然遇到，删除或重命名~/.Rprofile文件
- 或者设置环境变量`R_PROFILE_USER=""`

### 6. 断点续传如何工作？

**自动跳过的步骤：**
- Jellyfish Count: 检查`.jf`文件
- Jellyfish Histo: 检查`.jellyfish.histo`文件
- GenomeScope: 检查`model.txt`文件
- FastK: 检查`.ktab`文件
- Smudgeplot Hetmers: 检查`.smu`文件
- Smudgeplot Plot: 检查`.png`或`.pdf`文件

**重新运行方法：**
- 删除对应步骤的输出文件
- 或删除整个样品目录

## 性能优化建议 | Performance Optimization Tips

### 1. 内存优化 | Memory Optimization

```bash
# 减小Jellyfish哈希表
-s 5G  # 默认10G

# 减少FastK内存
--fastk-memory 50G  # 默认100G
```

### 2. 速度优化 | Speed Optimization

```bash
# 增加线程数
-t 64  # 默认12

# 跳过不需要的步骤
--skip-smudgeplot  # 仅运行GenomeScope
```

### 3. 批量分析优化 | Batch Analysis Optimization

```bash
# 使用相对路径减少路径长度
cd /path/to/project
biopytools smudgescope -i clean/ -o output/

# 确保输出目录在快速存储设备上（SSD）
-o /fast/storage/output/
```

## 引用 | Citation

如果使用本工具，请引用以下软件：

If you use this tool, please cite the following software:

- **Jellyfish**: Guillaume & Marçais  (2019)
- **GenomeScope 2.0**: Ranquet et al. (2022)
- **Smudgeplot**: R. Jaron & M. (2019)

## 更新日志 | Changelog

### v1.0.0 (2026-04-01)
- 初始版本发布
- 支持Jellyfish + GenomeScope + Smudgeplot完整流程
- 支持断点续传和批量分析
- 修复PATH污染问题
- 修复断点续传文件检查问题
- 添加结果文件汇总功能
