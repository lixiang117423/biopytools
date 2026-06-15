# SV多样性分析工具 (sv_diversity)

## 🧬 概述

SV多样性分析工具用于计算结构变异(Structural Variant, SV)的群体遗传学多样性指标，包括Kappa_S(丰富度)和Kappa_Pi(分歧度)等核心指标。

The SV Diversity Analysis Tool calculates population genetics diversity metrics for Structural Variants (SVs), including core metrics like Kappa_S (richness) and Kappa_Pi (divergence).

## 📋 功能特点

- **群体多样性分析**: 计算群体水平的多样性指标
- **个体差异评估**: 提供每个样本的详细统计信息
- **标准化支持**: 支持基于基因组大小的标准化指标
- **统计检验友好**: 输出格式便于进行T检验、ANOVA等统计分析
- **详细报告**: 自动生成指标说明和使用建议

## 🚀 使用方法

### 基本用法

```bash
# 基本分析
biopytools sv-diversity -i sv_matrix.csv -m metadata.txt

# 指定输出前缀
biopytools sv-diversity -i sv_matrix.csv -m metadata.txt -p my_analysis

# 指定基因组大小进行标准化
biopytools sv-diversity -i sv_matrix.csv -m metadata.txt -g 2500000000
```

### 高级用法

```bash
# 完整参数示例
biopytools sv-diversity \
    -i population_svs.csv \
    -m sample_groups.txt \
    -p human_diversity \
    -g 3000000000 \
    -o ./results
```

## 📁 输入文件格式

### SV矩阵文件 (CSV/TSV格式)

```
,Sample1,Sample2,Sample3,Sample4
SV1,1,0,1,1
SV2,0,1,0,0
SV3,1,1,1,0
SV4,0,0,1,1
```

- 第一行为样本ID
- 第一列为SV ID
- 矩阵值: 1=存在该SV, 0=不存在该SV

### 元数据文件 (CSV/TSV格式)

```
SampleID,Group
Sample1,Population_A
Sample2,Population_A
Sample3,Population_B
Sample4,Population_B
```

- 必须包含SampleID和Group两列
- SampleID必须与SV矩阵中的列名一致

## 📊 输出文件说明

### 1. {prefix}_group_stats.tsv (群体水平)

| 列名 | 说明 |
|------|------|
| Group | 群体名称 |
| N | 样本数量 |
| Kappa_S_Raw | 丰富度(原始值) |
| Kappa_Pi_Raw | 分歧度(原始值) |
| Kappa_S_PerBP | 标准化丰富度(SVs/bp) |
| Kappa_Pi_PerBP | 标准化分歧度(/bp) |

### 2. {prefix}_sample_stats.tsv (个体水平)

| 列名 | 说明 |
|------|------|
| SampleID | 样本ID |
| Group | 所属群体 |
| SV_Count | SV携带数量(负荷) |
| SV_Density_PerBP | SV密度 |
| Mean_Pairwise_Diff | 平均成对差异 |
| Mean_Pairwise_Diff_PerBP | 标准化平均成对差异 |

### 3. {prefix}_README.txt

详细的指标说明、计算公式和使用建议。

## 🧮 指标解释

### Kappa_S (丰富度指标)
- **定义**: 经样本大小校正的分离SV总数
- **公式**: `|SVs| / Harmonic_Number(n)`
- **解读**: 群体的变异"总库存"，反映群体遗传多样性总量

### Kappa_Pi (分歧度指标)
- **定义**: 两个随机个体间平均差异数
- **公式**: `Σ(2 * k * m) / (n * (n-1))`，其中k为携带该SV的样本数，m为未携带的样本数
- **解读**: 群体内个体间的平均遗传差异程度

### 标准化指标
- 当提供基因组大小时，会计算PerBP(每碱基对)标准化指标
- 便于不同基因组大小的物种间比较

## 📈 统计学应用

### 群体比较分析
```R
# R语言示例
group_data <- read.table("group_stats.tsv", header=TRUE, sep="\t")

# 群体丰富度比较
wilcox.test(Kappa_S_PerBP ~ Group, data=group_data)

# 群体分歧度比较
wilcox.test(Kappa_Pi_PerBP ~ Group, data=group_data)
```

### 个体差异分析
```R
# R语言示例
sample_data <- read.table("sample_stats.tsv", header=TRUE, sep="\t")

# 正态性检验
shapiro.test(sample_data$SV_Density_PerBP)

# T检验比较两组个体
t.test(SV_Density_PerBP ~ Group, data=sample_data)
```

### 可视化建议
- **箱线图**: 展示群体间多样性分布差异
- **散点图**: 显示个体多样性指标分布
- **热图**: 展示样本间的遗传距离

## ⚙️ 参数详解

| 参数 | 短参数 | 长参数 | 类型 | 默认值 | 说明 |
|------|--------|--------|------|--------|------|
| SV矩阵 | `-i` | `--input` | str | 必需 | SV矩阵文件路径(CSV/TSV) |
| 元数据 | `-m` | `--meta` | str | 必需 | 元数据文件路径 |
| 前缀 | `-p` | `--prefix` | str | sv_diversity | 输出文件前缀 |
| 基因组大小 | `-g` | `--genome-size` | int | None | 基因组大小(bp)，用于标准化 |
| 输出目录 | `-o` | `--output-dir` | str | ./sv_diversity_output | 输出目录 |

## 🔬 适用场景

- **群体遗传学研究**: 评估不同群体的遗传多样性
- **进化分析**: 研究种群历史和进化关系
- **比较基因组学**: 比较不同物种或品系的遗传差异
- **医学遗传学**: 评估患者群体的结构变异负荷
- **育种研究**: 分析育种群体的遗传多样性水平

## ⚠️ 注意事项

1. **数据质量**: 确保SV检测结果准确可靠
2. **样本数量**: 每个群体建议至少有2个样本
3. **基因组大小**: 提供准确的基因组大小以获得标准化指标
4. **样本对应**: 确保SV矩阵和元数据中的样本ID一致
5. **统计学考虑**: 小样本时建议使用非参数检验方法

## 🐛 故障排除

### 常见错误

1. **"没有匹配的样本"**
   - 检查SV矩阵的列名是否与元数据中的SampleID一致
   - 确认文件格式正确(CSV/TSV)

2. **"群体X样本数少于2"**
   - 确保每个群体至少有2个样本
   - 检查Group列的拼写是否一致

3. **内存不足**
   - 对于大型数据集，考虑减少SV数量或样本数量
   - 确保系统有足够的可用内存

## 📚 参考文献

本工具基于群体遗传学经典多样性指标开发，相关理论请参考：
- Population Genetics: A Concise Guide
- Principles of Population Genetics
- Molecular Evolution and Phylogenetics

## 📞 技术支持

如遇到问题，请检查：
1. 输入文件格式是否正确
2. 样本ID是否完全匹配
3. 系统资源是否充足
4. 参数设置是否合理

技术问题请联系生物信息学分析团队。