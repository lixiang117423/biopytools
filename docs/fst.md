# Fst遗传分化计算模块

**专业的群体遗传分化分析工具 | Professional Population Genetic Differentiation Analysis Tool**

## 功能概述 | Overview

Fst遗传分化计算模块是一个专业的群体遗传学分析工具，基于PLINK软件构建，提供从VCF文件到Fst统计量计算的完整流程。支持自动化的数据转换、灵活的质量控制和多种结果输出格式，适用于各种群体遗传学研究。

## 主要特性 | Key Features

- **完整分析流程**: VCF转换→质量控制→Fst计算三步骤自动化
- **灵活质控参数**: 可配置的MAF、缺失率、HWE等质控阈值
- **自动群体检测**: 智能解析群体文件，支持制表符和逗号分隔
- **多群体支持**: 支持两群体及多群体的Fst计算
- **双格式输出**: 长格式表格和矩阵格式结果
- **中间文件保留**: 可选择保留或删除中间处理文件
- **详细日志记录**: 完整的处理过程日志和错误追踪
- **高效处理**: 优化的处理流程，支持大规模基因组数据

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 计算两个群体间的Fst
biopytools fst \
    -i variants.vcf \
    -p population.txt \
    -o fst_output

# 使用自定义质控参数
biopytools fst \
    -i variants.vcf \
    -p population.txt \
    -o fst_output \
    --maf 0.01 \
    --geno 0.2
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --vcf-file` | VCF文件路径 | `-i variants.vcf` |
| `-p, --pop-file` | 群体文件路径（样本ID + 群体标签） | `-p population.txt` |
| `-o, --output-dir` | 输出目录 | `-o fst_output` |

### 群体文件格式 | Population File Format

群体文件应为两列格式，第一列为样本ID，第二列为群体标签：

```text
sample1    POP1
sample2    POP1
sample3    POP2
sample4    POP2
```

支持以下分隔符：
- 制表符（`\t`）
- 逗号（`,`）
- 空格（` `）

### 质量控制参数 | Quality Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--maf` | `0.05` | 最小等位基因频率阈值 |
| `--geno` | `0.1` | 位点缺失率阈值（call rate < 1-geno） |
| `--mind` | `0.1` | 样本缺失率阈值（call rate < 1-mind） |
| `--hwe` | `1e-6` | Hardy-Weinberg平衡p值阈值 |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--plink-path` | `auto-detect` | PLINK软件路径（自动检测conda环境） |

### 输出控制 | Output Control

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-keep-intermediate` | `False` | 不保留中间文件 |

## 输入文件格式 | Input File Formats

### VCF文件 | VCF File

标准VCF格式的变异数据：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000000>
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2  sample3  sample4
chr1    1000 .   A   G   45.2  PASS    DP=30  GT:DP   0/0:30  0/1:25  1/1:28  0/0:32
chr1    2000 .   T   C   62.8  PASS    DP=25  GT:DP   0/1:25  1/1:22  0/0:28  0/1:30
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

#### 直接计算模式 | Direct Calculation Mode

```
fst_output/
├── 00_intermediate/          # 中间文件（如果保留）
│   ├── mydata.bed             # PLINK binary文件
│   ├── mydata.bim
│   ├── mydata.fam
│   ├── mydata_qc.bed          # 质控后文件（如果启用质控）
│   ├── mydata_qc.bim
│   └── mydata_qc.fam
├── fst_result.fst             # PLINK Fst原始输出
├── population_filtered.plink  # PLINK格式群体文件
├── fst_long_format.txt        # 长格式结果
├── fst_matrix.txt             # 矩阵格式结果
└── 99_logs/
    └── fst_calculation.log    # 运行日志
```

#### Bootstrap模式 | Bootstrap Mode

```
fst_output/
├── 00_intermediate/          # 中间文件（如果保留）
│   ├── mydata.bed
│   ├── mydata.bim
│   ├── mydata.fam
│   ├── mydata_qc.bed          # 质控后文件（如果启用质控）
│   ├── mydata_qc.bim
│   └── mydata_qc.fam
├── bootstrap_iterations/      # Bootstrap迭代文件
│   ├── fst_iter1.fst          # 第1次迭代的Fst结果
│   ├── fst_iter1.log
│   ├── fst_iter1.nosex
│   ├── population_iter1.plink
│   ├── fst_iter2.fst          # 第2次迭代的Fst结果
│   ├── ...
│   └── fst_iter100.fst        # 第100次迭代的Fst结果
├── all_fst_results.txt        # Bootstrap统计汇总（主要结果文件）
└── 99_logs/
    └── fst_calculation.log    # 运行日志
```

### 输出文件说明 | Output File Description

#### 1. all_fst_results.txt（Bootstrap模式主要结果文件）

Bootstrap统计汇总文件，包含所有群体对的Fst统计量：

```text
Population1    Population2    Mean_Fst    Std      Min      Max      Median
China_North    Korea          0.023456    0.001234 0.020123 0.026789 0.023456
China_North    Japan          0.019876    0.001567 0.017234 0.023456 0.019876
```

**列说明**：
- `Population1`, `Population2`: 群体对名称
- `Mean_Fst`: 100次迭代的平均Fst值
- `Std`: 标准差
- `Min`: 最小Fst值
- `Max`: 最大Fst值
- `Median`: 中位数

#### 2. fst_long_format.txt（直接计算模式）

长格式表格，每行表示一对群体的Fst值：

```text
Population1    Population2    Fst
POP1           POP2           0.123456
```

#### 2. fst_matrix.txt（矩阵格式）

矩阵格式，展示所有群体对之间的Fst值：

```text
        POP1        POP2
POP1    0.000000    0.123456
POP2    0.123456    0.000000
```

#### 3. fst_result.fst（PLINK原始输出）

PLINK生成的原始Fst文件，包含每个SNP位点的Fst值。

## 使用示例 | Usage Examples

### 示例1：基本Fst计算

```bash
# 计算两个野生群体的Fst
biopytools fst \
    -i wild_population.vcf \
    -p wild_pop_groups.txt \
    -o wild_fst_results
```

### 示例2：严格质控参数

```bash
# 使用严格的质控参数
biopytools fst \
    -i variants.vcf \
    -p population.txt \
    -o strict_fst \
    --maf 0.01 \
    --geno 0.05 \
    --mind 0.05 \
    --hwe 1e-10
```

### 示例3：多群体比较

```bash
# 计算多个地理群体间的Fst
biopytools fst \
    -i geo_populations.vcf \
    -p geo_groups.txt \
    -o geo_fst_results \
    --maf 0.02
```

### 示例4：不保留中间文件

```bash
# 节省存储空间，不保留中间文件
biopytools fst \
    -i variants.vcf \
    -p population.txt \
    -o fst_results \
    --no-keep-intermediate
```

## 分析流程 | Analysis Pipeline

### 步骤1：VCF转PLINK格式

将VCF文件转换为PLINK binary格式（.bed/.bim/.fam）

### 步骤2：质量控制

根据指定参数进行质量控制：
- 过滤低频变异（MAF阈值）
- 过滤高缺失率位点
- 过滤高缺失率样本
- 过滤偏离HWE的位点

### 步骤3：Fst计算

使用Weir & Cockerham (1984)方法计算Fst值：
- 支持两群体比较
- 支持多群体两两比较
- 生成逐位点Fst值和全基因组平均值

### 步骤4：结果格式化

- 生成长格式表格（便于统计分析）
- 生成矩阵格式（便于可视化）

## 注意事项 | Important Notes

1. **群体文件格式**: 确保群体文件的样本ID与VCF文件中的样本名完全一致
2. **质控参数**: 根据数据质量调整质控参数，过严可能导致可用位点过少
3. **样本量**: 每个群体至少需要2个样本
4. **内存需求**: 大规模VCF文件需要足够的内存（建议16GB以上）
5. **计算时间**: Fst计算时间与样本量和位点数成正比

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "样本不在VCF文件中"错误**
```bash
# 检查群体文件的样本ID是否与VCF文件一致
bcftools query -l variants.vcf | sort > vcf_samples.txt
cut -f1 population.txt | sort > pop_samples.txt
diff vcf_samples.txt pop_samples.txt
```

**Q: Fst值为NA或缺失**
```bash
# 可能原因：
# 1. 质控后可用位点过少 -> 降低--maf阈值
# 2. 样本量不足 -> 检查每个群体的样本数
# 3. 群体标签错误 -> 检查群体文件格式
```

**Q: 计算时间过长**
```bash
# 解决方案：
# 1. 使用更严格的质控参数减少位点数
# 2. 考虑使用--no-keep-intermediate节省I/O时间
# 3. 在高性能计算集群上运行
```

## 方法说明 | Method Description

### Fst计算方法

本工具使用Weir & Cockerham (1984)方法计算Fst值：

- **方法**: Weir & Cockerham's Fst
- **引用**: Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. Evolution, 38(6), 1358-1370.
- **特点**: 考虑了样本量和群体大小的影响，适用于各种群体遗传学研究

### Fst值解释

- **Fst < 0.05**: 遗传分化很小
- **0.05 ≤ Fst < 0.15**: 遗传分化中等
- **0.15 ≤ Fst < 0.25**: 遗传分化较大
- **Fst ≥ 0.25**: 遗传分化很大

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **PLINK** (版本 1.9 或更新)
  - 安装: conda install -c bioconda plink
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上）
- **RAM**: 最少4GB（推荐16GB以上）
- **存储**: 预留VCF文件大小3倍的磁盘空间

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关文献：

```
Weir, B.S. & Cockerham, C.C. (1984).
Estimating F-statistics for the analysis of population structure.
Evolution, 38(6), 1358-1370.
```

同时建议引用PLINK软件：

```
Purcell, S. et al. (2007).
PLINK: a tool set for whole-genome association and population-based linkage analyses.
American Journal of Human Genetics, 81(3), 559-575.
```
