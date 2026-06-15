# Dsuite Dtrios基因渗入分析工具 (dsuit)

## 🧬 概述

Dsuite Dtrios分析工具用于进行ABBA-BABA检验(D-statistics)，检测群体间的基因渗入事件。该工具基于Dsuite软件，计算D统计量、f4-ratio等群体遗传学指标，用于评估不同群体间的基因流模式。

The Dsuite Dtrios analysis tool performs ABBA-BABA tests (D-statistics) to detect introgression events between populations. This tool is based on Dsuite software and calculates population genetics metrics like D-statistics and f4-ratio to evaluate gene flow patterns between different populations.

## 📋 功能特点

- **ABBA-BABA检验**: 计算D统计量检测基因渗入
- **f4-ratio估计**: 量化基因流比例
- **多群体分析**: 自动计算所有可能的三元组组合
- **灵活过滤**: 支持基于等位基因数和变异类型的VCF过滤
- **系统发育树整合**: 可选的基于树结构的分析
- **Jackknife重采样**: 提供统计显著性评估
- **ABBA聚类分析**: 检测高度分化物种间的基因流
- **详细报告**: 生成完整的分析报告和可视化建议

## 🚀 使用方法

### 基本用法

```bash
# 基本分析
biopytools dsuit -i variants.vcf.gz -s sets.txt -p analysis1

# 指定过滤参数
biopytools dsuit -i variants.vcf.gz -s sets.txt -m 2 -M 2 -v snps

# 使用系统发育树
biopytools dsuit -i variants.vcf.gz -s sets.txt -t tree.nwk
```

### 高级用法

```bash
# 完整参数示例
biopytools dsuit \
    -i population_variants.vcf.gz \
    -s population_sets.txt \
    -p introgression_analysis \
    -t phylogeny.nwk \
    -k 30 \
    --ABBAclustering \
    -d /path/to/Dsuite \
    -o ./results

# ABBA聚类分析 (适用于高度分化物种)
biopytools dsuit \
    -i ancient_genomes.vcf.gz \
    -s species_sets.txt \
    -p deep_time_analysis \
    --ABBAclustering \
    -k 25
```

## 📁 输入文件格式

### 1. VCF文件

- 支持压缩和未压缩格式 (.vcf, .vcf.gz)
- 必须包含双等位SNP位点
- 样本名称必须与分组文件一致

### 2. 分组文件 (SETS.txt)

```
Sample1    Species1
Sample2    Species1
Sample3    Species2
Sample4    Species2
Sample5    Outgroup
Sample6    Outgroup
```

- **格式**: 两列，用制表符分隔
- **第一列**: 样本ID
- **第二列**: 群体/物种名称
- **外群**: 必须包含名为"Outgroup"的群体

## 📊 输出文件说明

### 1. 主要结果文件

| 文件名 | 说明 |
|--------|------|
| `{prefix}_BBAA.txt` | D统计量、f4-ratio、显著性检验结果 |
| `{prefix}_Dmin.txt` | 最小D统计量结果 |
| `{prefix}_tree.txt` | 基于系统发育树的结果 (如果提供树文件) |
| `{prefix}_combine.txt` | 组合分析结果 (用于DtriosCombine) |
| `{prefix}_analysis_summary.txt` | 分析总结报告 |

### 2. BBAA文件列解释

| 列名 | 说明 |
|------|------|
| P1 | 群体1 |
| P2 | 群体2 |
| P3 | 群体3 |
| BBAA | BBAA模式计数 |
| BABA | BABA模式计数 |
| ABBA | ABBA模式计数 |
| Dstat | D统计量 |
| Zscore | Z分数 (统计显著性) |
| pvalue | p值 |
| f4ratio | f4-ratio (基因流比例) |

## 🧮 指标解释

### D统计量 (D-statistic)
- **范围**: -1 到 +1
- **解读**:
  - D > 0: P1-P2间存在基因流
  - D < 0: P3与P1/P2间存在基因流
  - D ≈ 0: 不符合基因流模式或符合树状进化

### f4-ratio
- **范围**: 0 到 1
- **解读**: 基因流比例估计
  - 0: 无基因流
  - 0.1: 10%的基因组来源于基因流
  - 1: 完全的基因流

### 统计显著性
- **Z-score**: |Z| > 3 通常认为统计显著
- **p-value**: p < 0.05 通常认为显著

### ABBA聚类分析
- **目的**: 区分真实的基因流和同塑性产生的假阳性
- **适用**: 高度分化物种间、远古基因流事件
- **p-values**:
  - `clustering_sensitive`: 更高统计功效，可能假阳性
  - `clustering_robust`: 较低功效，更稳健

## ⚙️ 参数详解

| 参数 | 短参数 | 长参数 | 类型 | 默认值 | 说明 |
|------|--------|--------|------|--------|------|
| VCF文件 | `-i` | `--vcf` | str | 必需 | VCF文件路径 |
| 分组文件 | `-s` | `--sets` | str | 必需 | 群体分组文件路径 |
| 前缀 | `-p` | `--prefix` | str | dsuite_analysis | 输出文件前缀 |
| Dsuite路径 | `-d` | `--dsuite-bin` | str | /share/org/.../Dsuite | Dsuite可执行文件 |
| bcftools | `-b` | `--bcftools-bin` | str | bcftools | bcftools命令路径 |
| 最小等位基因 | `-m` | `--min-alleles` | int | 2 | 最小等位基因数 |
| 最大等位基因 | `-M` | `--max-alleles` | int | 2 | 最大等位基因数 |
| 变异类型 | `-v` | `--variant-type` | choice | snps | snps/indels/both |
| Jackknife块数 | `-k` | `--jk-num` | int | 20 | Jackknife重采样块数 |
| Jackknife窗口 | | `--jk-window` | int | None | Jackknife窗口大小 |
| 系统发育树 | `-t` | `--tree` | str | None | Newick格式树文件 |
| 运行名称 | `-n` | `--run-name` | str | None | 分析运行名称 |
| 禁用f4-ratio | | `--no-f4-ratio` | flag | False | 不计算f4-ratio |
| ABBA聚类 | | `--ABBAclustering` | flag | False | 启用ABBA聚类分析 |
| 输出目录 | `-o` | `--output-dir` | str | ./dsuite_output | 输出目录 |

## 🔬 适用场景

- **群体遗传学研究**: 检测物种间的基因渗入事件
- **进化生物学**: 重建物种演化历史和基因流模式
- **保护遗传学**: 评估杂交对遗传多样性的影响
- **作物育种**: 分析栽培品种与野生种的基因交流
- **人类进化**: 研究现代人与古人类间的基因流
- **病原体研究**: 追踪病原菌的传播和重组事件

## 📈 结果可视化

### R脚本示例
```R
# 读取D统计结果
library(ggplot2)

results <- read.table("dsuite_analysis_BBAA.txt", header=TRUE, sep="\t")

# D统计量分布图
ggplot(results, aes(x=Dstat)) +
  geom_histogram(bins=50, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  labs(title="D-statistic Distribution",
       x="D-statistic", y="Count") +
  theme_minimal()

# 显著性检验
significant <- results[abs(results$Zscore) > 3, ]
print(paste("Significant trios:", nrow(significant)))

# f4-ratio与D统计量关系
ggplot(results, aes(x=Dstat, y=f4ratio)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se=TRUE) +
  labs(title="f4-ratio vs D-statistic",
       x="D-statistic", y="f4-ratio") +
  theme_minimal()
```

### 系统发育树可视化
```bash
# 使用Dsuite自带的绘图工具
/path/to/Dsuite/utils/dtools.py fbranch_results.txt phylogeny.nwk
```

## ⚠️ 注意事项

### 数据质量要求
1. **VCF格式**: 确保VCF文件格式正确，样本名称一致
2. **群体分组**: 必须包含外群(Outgroup)
3. **样本数量**: 每个群体建议至少有2-3个样本
4. **变异数量**: 建议至少有数千个高质量SNP位点

### 统计考量
1. **Jackknife块数**: 对于全基因组分析，建议至少20个块
2. **多重检验**: 进行多个三元组分析时考虑多重检验校正
3. **连锁不平衡**: 考虑位点间的连锁关系
4. **样本量**: 小样本时谨慎解释结果

### 生物学解释
1. **D统计量解读**: 考虑系统发育关系的影响
2. **基因流方向**: D统计量的符号指示基因流方向
3. **时间尺度**: 不同时间尺度的基因流需要不同解释方法
4. **地理因素**: 结合地理分布解释基因流模式

## 🐛 故障排除

### 常见错误

1. **"外群不存在"**
   - 检查SETS.txt文件中是否包含"Outgroup"群体
   - 确认外群样本的命名准确

2. **"过滤后无变异位点"**
   - 检查VCF文件的变异质量
   - 调整过滤参数 (m, M, v)
   - 确认VCF文件包含双等位SNP

3. **"Dsuite不可执行"**
   - 检查Dsuite安装路径
   - 确认文件具有执行权限
   - 使用-d参数指定正确路径

4. **"内存不足"**
   - 减少样本数量或位点数量
   - 增加Jackknife块数减少内存使用
   - 使用分区域分析后合并结果

### 性能优化

1. **并行化**: 对于大型数据集，考虑使用DtriosParallel
2. **分区域分析**: 按染色体或基因组区域分别分析后合并
3. **内存管理**: 适当调整过滤条件减少数据量
4. **存储空间**: 确保有足够磁盘空间存储中间文件

## 📚 参考文献

1. **Dsuite原始论文**:
   Malinsky, M., Matschiner, M. & Svardal, H. (2021) Dsuite ‐ fast D‐statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21, 584–595.

2. **ABBA-BABA检验理论**:
   Patterson, N., et al. (2012) Ancient admixture in human history. Genetics 192, 1065-1093.

3. **ABBA聚类分析**:
   Koppetsch, T., Malinsky, M. & Matschiner, M. (2024) Towards Reliable Detection of Introgression in the Presence of Among-Species Rate Variation. Systematic Biology, syae028.

4. **f-branch方法**:
   Malinsky, M., et al. (2018) Evolution of genomic diversity in the African cichlid fish adaptive radiation. Nature Ecology & Evolution 2, 1217-1228.

## 📞 技术支持

如遇到问题，请检查：
1. 输入文件格式和路径是否正确
2. Dsuite和bcftools是否正确安装
3. 系统资源是否充足 (内存、磁盘空间)
4. 参数设置是否合理

技术问题请联系生物信息学分析团队或参考Dsuite官方文档。