# rMVP GWAS分析工具 | rMVP GWAS Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2026-04-01

## 概述 | Overview

rMVP工具是基于**rMVP R包**的全基因组关联分析（GWAS）Python包装器，支持多种统计模型（GLM、MLM、FarmCPU），提供完整的GWAS分析流程和可视化功能。

The rMVP tool is a Python wrapper for the **rMVP R package** for Genome-Wide Association Studies (GWAS), supporting multiple statistical models (GLM, MLM, FarmCPU), providing a complete GWAS analysis pipeline and visualization functionality.

## 功能特点 | Features

- 🧬 **多模型支持**: 支持GLM、MLM、FarmCPU三种GWAS模型
- 📊 **完整流程**: VCF数据转换、GWAS分析、结果可视化一站式完成
- 🗺️ **多种可视化**: 曼哈顿图、QQ图、环形图等多种图表类型
- 🎯 **批量分析**: 支持多表型批量分析
- 💾 **断点续传**: 自动跳过已完成的步骤，支持中断恢复
- ⚡ **高性能**: 支持多线程并行计算
- 📈 **结果解析**: 自动提取显著信号位点和统计信息
- 🔧 **灵活配置**: 丰富的参数配置选项

## 分析模型 | Analysis Models

### GLM (General Linear Model | 一般线性模型)

**特点 | Characteristics:**
- 计算速度快
- 适用于初步探索性分析
- 不考虑群体结构和亲缘关系

**适用场景 | Use Cases:**
- 样本量较小（< 1000）
- 群体结构简单
- 快速筛查候选基因

### MLM (Mixed Linear Model | 混合线性模型)

**特点 | Characteristics:**
- 控制群体结构和亲缘关系
- 假阳性率低
- 计算速度较慢

**适用场景 | Use Cases:**
- 样本量较大
- 存在明显的群体结构
- 需要精确的关联分析

### FarmCPU (Fixed and random model Circulating Probability Unification)

**特点 | Characteristics:**
- 固定效应和随机效应循环使用
- 在速度和准确性之间取得平衡
- 参数优化较为复杂

**适用场景 | Use Cases:**
- 大规模数据集（> 10K样本）
- 需要兼顾速度和准确性
- 推荐作为主要分析方法

## 安装和使用 | Installation and Usage

### 前置要求 | Prerequisites

#### R环境要求 | R Environment Requirements

```bash
# 创建conda环境
conda create -n rMVP r-base
conda activate rMVP

# 安装rMVP包
R -e "install.packages('rMVP', repos='https://cloud.r-project.org')"
```

#### Python环境要求 | Python Environment Requirements

```bash
# biopytools会自动处理Python依赖
# 确保biopytools已安装
pip install biopytools
```

### 作为biopytools模块使用 | Using as biopytools module

```bash
# 基本用法
biopytools rmvp -i input.vcf -p phenotype.txt -o output/

# 使用未压缩VCF（推荐）
biopytools rmvp -i input.vcf -p phenotype.txt -o output/

# 指定模型和线程数
biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -m GLM FarmCPU -t 24

# 使用特定conda环境
biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -r my_r_env
```

### 作为Python模块使用 | Using as Python module

```python
from biopytools.rmvp import RMVPConfig, RMVPAnalyzer

# 创建配置
config = RMVPConfig(
    vcf_file="input.vcf",
    pheno_file="phenotype.txt",
    output_dir="output",
    models=["GLM", "MLM", "FarmCPU"],
    ncpus=24
)

# 运行分析
analyzer = RMVPAnalyzer(config)
success = analyzer.run_analysis()
```

## 命令行参数 | Command Line Arguments

### 必需参数 | Required Arguments

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --vcf` | VCF格式基因型文件（支持压缩但推荐未压缩） | `-i input.vcf` |
| `-p, --pheno` | 表型文件（TSV格式） | `-p phenotype.txt` |
| `-o, --output-dir` | 输出目录 | `-o output/` |

### 可选参数 | Optional Arguments

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-m, --models` | `GLM,MLM,FarmCPU` | 分析模型（逗号分隔） | `-m GLM FarmCPU` |
| `-t, --threads` | `12` | CPU核心数 | `-t 24` |
| `-r, --r-env` | `rMVP` | conda环境名称 | `-r my_r_env` |
| `--log-level` | `INFO` | 日志级别 | `--log-level DEBUG` |

## 表型文件格式 | Phenotype File Format

表型文件应为TSV格式（制表符分隔）：

```
SampleID	Trait1
Sample1	0.85
Sample2	0.72
Sample3	0.91
Sample4	NA
```

**要求 | Requirements:**
- 第一列：样本ID（必须与VCF文件中的样本名一致）
- 第二列起：表型值（数值型，缺失值用NA表示）
- 分隔符：制表符（`\t`）

## 输出文件 | Output Files

### 数据转换输出 | Data Conversion Output

- `{prefix}.geno.bin` - 基因型二进制文件
- `{prefix}.geno.desc` - 基因型描述文件
- `{prefix}.geno.map` - SNP位置信息
- `{prefix}.geno.ind` - 样本信息
- `{prefix}.phe` - 表型数据
- `{prefix}.kin.bin` / `{prefix}.kin.desc` - Kinship矩阵
- `{prefix}.pc.bin` / `{prefix}.pc.desc` - PCA结果

### GWAS分析输出 | GWAS Analysis Output

#### 数据文件 | Data Files

- `{prefix}_{trait}.{model}.csv` - 完整SNP结果
  - 列：SNP, CHROM, POS, A1, A2, MAF, Effect, SE, p值
  - 示例：`RMVP_Result_DI.glm.csv`

- `{prefix}_{trait}_signals.{model}.csv` - 显著信号位点
  - 示例：`RMVP_Result_DI_signals.glm.csv`

- `{prefix}_{trait}.RData` - 单表型分析结果对象
- `{prefix}_all_results.RData` - 所有表型分析结果

#### 图片文件 | Figure Files

- `{prefix}_{trait}.{model}.Rectangular-Manhattan.{trait}.jpg` - 矩形曼哈顿图
- `{prefix}_{trait}.{model}.QQplot.{trait}.jpg` - QQ图
- `{prefix}_{trait}.GLM.{trait}.MLM.{trait}.FarmCPU.Circular-Manhattan.{trait}.jpg` - 环形曼哈顿图
- `{prefix}_{trait}.GLM.{trait}.MLM.{trait}.FarmCPU.Multracks-Manhattan.{trait}.jpg` - 多轨道曼哈顿图
- `{prefix}_{trait}.Phe_Dist.{trait}.jpg` - 表型分布图
- `{prefix}_{trait}.{trait}.PCA_2D.jpg` - PCA二维图

## 断点续传功能 | Checkpoint Resume Feature

本工具支持断点续传，自动跳过已完成的步骤：

### 检查机制 | Check Mechanism

- **数据转换步骤**: 检查 `{prefix}.geno.desc` 文件
- **GWAS分析步骤**: 检查所有表型-模型组合的 `.pmap` 文件

### 使用示例 | Usage Example

```bash
# 第一次运行 - 执行完整流程
biopytools rmvp -i input.vcf -p phenotype.txt -o output/
# 输出：[3/5] 数据转换|Data conversion

# 中断后重新运行 - 自动跳过已完成步骤
biopytools rmvp -i input.vcf -p phenotype.txt -o output/
# 输出：跳过已完成步骤|Skipping completed step: 数据转换|Data conversion
```

## 常见问题 | FAQ

### 1. VCF文件支持压缩格式吗？

**问题 | Question**: VCF文件支持压缩格式吗？

**答 | Answer**: 代码层面支持`.vcf.gz`，但rMVP包内部不支持压缩格式。

**解决方案 | Solution**:
```bash
# 推荐使用未压缩的VCF文件
gunzip -c input.vcf.gz > input.vcf
biopytools rmvp -i input.vcf -p phenotype.txt -o output/
```

### 2. 如何解读GWAS结果？

**问题 | Question**: 如何解读GWAS结果？

**答 | Answer**: 主要查看以下文件：

1. **曼哈顿图** - 识别显著关联的染色体区域
   - X轴：染色体位置
   - Y轴：-log10(p值)
   - 显著阈值线：通常为 -log10(1e-5) 到 -log10(1e-8)

2. **QQ图** - 评估模型拟合度
   - X轴：期望p值
   - Y轴：观察p值
   - 越接近对角线说明模型拟合越好

3. **signals CSV文件** - 提取显著信号位点
   - 通常筛选条件：p < 1e-5
   - 包含SNP位置、效应值、MAF等信息

### 3. 如何自定义绘图？

**问题 | Question**: 如何自定义绘图？

**答 | Answer**: 使用RData文件重新绘图：

```r
# 加载结果
load("output/RMVP_Result_DI.RData")

# 使用rMVP内置绘图功能
library(rMVP)

# 重新绘制曼哈顿图
MVP.Report(
    phe = imvp$phe,
    geno = imvp$geno,
    map = imvp$map,
    imvp_results = imvp,
    file.type = "pdf",  # 输出PDF格式
    dpi = 600
)
```

### 4. 内存不足怎么办？

**问题 | Question**: 内存不足怎么办？

**答 | Answer**:

1. **减少线程数**
   ```bash
   biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -t 8
   ```

2. **只运行单个模型**
   ```bash
   biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -m FarmCPU
   ```

3. **使用计算节点**
   - 大型数据集建议在专用计算节点运行
   - 避免在登录节点运行

### 5. 分析需要多长时间？

**问题 | Question**: 分析需要多长时间？

**答 | Answer**:

| 数据规模 | CPU核心数 | 数据转换 | GWAS分析 | 总计 |
|---------|----------|----------|----------|------|
| 1K样本，100K SNP | 12 | 5分钟 | 5分钟 | 10分钟 |
| 10K样本，1M SNP | 24 | 30分钟 | 30分钟 | 1小时 |
| 50K样本，10M SNP | 48 | 2小时 | 3小时 | 5小时 |

**注意**: 使用断点续传功能后，重复运行可跳过数据转换步骤。

## 结果解读指南 | Result Interpretation Guide

### 曼哈顿图解读 | Manhattan Plot Interpretation

1. **X轴**: 染色体位置
2. **Y轴**: -log10(p值)，值越大表示关联性越强
3. **阈值线**:
   - 蓝色：建议阈值（p = 1e-5）
   - 红色：严格阈值（p = 1e-8）
4. **显著峰**: 超过阈值线的峰值区域

### QQ图解读 | QQ Plot Interpretation

1. **X轴**: 期望p值（Expected p-value）
2. **Y轴**: 观察p值（Observed p-value）
3. **对角线**: 理论分布
4. **偏离程度**:
   - 轻微偏离：正常
   - 严重偏离：可能存在群体结构或假阳性

### 信号位点筛选 | Signal Site Filtering

推荐筛选标准：
```r
# 加载signals文件
signals <- read.csv("RMVP_Result_DI_signals.glm.csv")

# 筛选显著信号（p < 1e-5）
significant <- signals[signals$DI.GLM < 1e-5, ]

# 进一步筛选MAF > 0.05
filtered <- significant[significant$MAF > 0.05, ]

# 按p值排序
sorted <- filtered[order(filtered$DI.GLM), ]
```

## 性能优化建议 | Performance Optimization Recommendations

### 1. 数据准备优化 | Data Preparation Optimization

- 使用未压缩的VCF文件
- 确保样本ID在VCF和表型文件中一致
- 预先过滤低质量SNP（可选）

### 2. 计算资源配置 | Computing Resource Configuration

| 数据规模 | 推荐CPU核心数 | 推荐内存 |
|---------|--------------|---------|
| < 1K样本 | 12 | 16 GB |
| 1K-10K样本 | 24 | 32 GB |
| > 10K样本 | 48+ | 64 GB+ |

### 3. 模型选择策略 | Model Selection Strategy

**初筛阶段**:
```bash
# 只运行GLM快速筛查
biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -m GLM
```

**精细分析**:
```bash
# 使用FarmCPU精确分析
biopytools rmvp -i input.vcf -p phenotype.txt -o output/ -m FarmCPU -t 48
```

## 依赖项 | Dependencies

### R包依赖 | R Package Dependencies

- rMVP >= 1.0.0
- bigmemory
- 其他依赖由rMVP自动安装

### Python依赖 | Python Dependencies

- Python >= 3.8
- pathlib
- dataclasses
- logging
- subprocess

### 系统依赖 | System Dependencies

- conda >= 4.0
- R >= 4.0.0

## 参考资源 | References

### rMVP相关

- **rMVP GitHub**: https://github.com/zhangwenyu931120/rMVP
- **rMVP论文**: https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.15451
- **rMVP文档**: https://github.com/zhangwenyu931120/rMVP/wiki

### GWAS分析方法

- **GWAS方法综述**: https://www.nature.com/articles/nrg3021
- **MLM方法**: https://www.nature.com/articles/ng.698
- **FarmCPU方法**: https://www.genetics.org/content/206/3/1393

## 更新日志 | Changelog

### v1.0.0 (2026-04-01)

- 初始版本发布
- 支持GLM、MLM、FarmCPU三种模型
- 实现断点续传功能
- 完整的日志分离（stdout/stderr）
- 自动R环境检测
- 结果文件自动解析和整合

## 许可证 | License

本工具遵循biopytools项目的许可证。

## 联系方式 | Contact

作者 | Author: Xiang LI <lixiang117423@gmail.com>

项目地址 | Project: https://github.com/your-org/biopytools
