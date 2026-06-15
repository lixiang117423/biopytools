# GWAS GEC分析工具 - 使用指南|GWAS GEC Analysis Tool - User Guide

## 🎯 功能简介

基于GEC（Genome-wide Error Correction）算法，计算GWAS研究的**有效独立检验数**和**校正后的显著性阈值**。

**使用VCF格式的参考文件进行LD结构分析和阈值计算！**

Calculate effective independent test numbers and adjusted significance thresholds for GWAS studies using GEC algorithm.

---

## ✨ 核心优势

- ✅ **简单易用**：提供VCF文件和GWAS P值即可运行
- ✅ **基于LD结构**：更准确的显著性阈值校正
- ✅ **自动计算**：自动计算有效检验数和阈值
- ✅ **详细报告**：包含染色体级别统计信息

---

## 🚀 快速开始

```bash
biopytools gwas-gec \
  -i gwas_domain_p_file.txt \
  -r input.vcf.gz \
  -o gec_output \
  -t 64 \
  -m 100g
```

**工具会自动完成：**
1. ✅ 读取VCF格式的参考基因组数据
2. ✅ 运行GEC算法计算有效检验数
3. ✅ 计算校正后的显著性阈值
4. ✅ 生成详细的汇总报告

---

## 📋 准备工作

### 需要两个文件：

#### 1. GWAS P值汇总统计文件

必需列（列名可自定义）：

---

## 📝 输入文件

### 1. GWAS P值文件

必需列（列名可自定义）：

| 列名|Column | 说明|Description | 示例|Example |
|---------|---------|----------------|-------------------|
| CHR | 染色体|Chromosome | Chr01, Chr02 ... |
| BP | 物理位置|Physical position | 9572, 752566 |
| P | P值|P-value | 0.9799, 1.2e-5 |

**自定义列名示例：**
```bash
biopytools gwas-gec -i results.txt -r input.vcf.gz \
  --chrom-col Chr --pos-col Pos --p-col Pval
```

#### 2. VCF格式的参考基因组文件

- `input.vcf.gz` 或 `input.vcf`
- 包含研究群体的基因型数据
- 用于计算LD结构

---

## 🔧 参数说明

### 必需参数

| 参数|Parameter | 说明|Description |
|------------|-------------|----------------------------------------|
| `-i, --pfile` | GWAS P值文件|GWAS P-value file |
| `-r, --reference` | 参考VCF文件|Reference VCF file |

### 可选参数

| 参数|默认值|说明|
|------------|-------------|-------------|
| `-o, --output-dir` | `./gec_output` | 输出目录|
| `-t, --threads` | `16` | 线程数|
| `-m, --memory` | `8g` | Java内存分配（根据数据量调整）|
| `--maf-filter` | `0.05` | MAF过滤阈值|
| `--alpha` | `0.05` | 显著性水平(FWER)|
| `--chrom-col` | `CHR` | 染色体列名|
| `--pos-col` | `BP` | 位置列名|
| `--p-col` | `P` | P值列名|

---

## 📊 输出文件

### 1. GEC原始输出

- `gwas_gec.effective.size.txt.gz` - 每个LD块的有效检验数

**列说明：**
- `Chrom`: 染色体
- `StartPos`: LD块起始位置
- `EndPos`: LD块终止位置
- `Num`: LD块内SNP数量
- `EffectiveNum`: 有效独立检验数

### 2. 汇总报告

- `gwas_gec_summary.txt` - 人类可读的汇总报告

**包含内容：**
- 总LD块数
- 总有效检验数
- 校正后显著性阈值
- 各染色体统计

### 3. 日志文件

- `gwas_gec.log` - 详细运行日志

---

## 📈 结果解读

### 阈值计算公式

```
校正后阈值 = Alpha / 总有效检验数
Adjusted Threshold = Alpha / Σ EffectiveNum
```

### 示例

假设分析结果：
- 总有效检验数：**456,789.12**
- Alpha水平：**0.05**

**校正后阈值：**
```
阈值 = 0.05 / 456,789.12 = 1.09e-07
```

### 与传统Bonferroni对比

| 方法|Method | 阈值|Threshold | 说明|Note |
|---------|------------------|-------------|
| 传统Bonferroni | 5×10⁻⁸ | 过于保守|Too conservative |
| **GEC校正** | **1.09×10⁻⁷** | **更准确|More accurate** |

**结论**：GEC校正后阈值放宽约5倍，能发现更多真实关联！

---

## 💡 实用建议

### 1. 内存分配

根据VCF文件大小调整内存：

- 小数据集（<50万SNP）：4-8g
- 大数据集（50-100万SNP）：16-32g
- 超大数据集（>100万SNP）：64g+

```bash
# 示例：大数据集使用更多内存
biopytools gwas-gec \
  -i gwas.txt \
  -r input.vcf.gz \
  -t 64 \
  -m 100g
```

### 2. 染色体格式一致性

确保P值文件和VCF文件的染色体命名格式一致：

**常见问题：**
- P值文件使用 `Chr01`，VCF使用 `1` → 不匹配
- P值文件使用 `1`，VCF使用 `Chr01` → 不匹配

**解决方案：**
保持两者格式一致，或使用 `--chrom-col` 和 `--pos-col` 参数指定列名

### 3. 输出文件说明

- 小数据集（<50万SNP）：4-8g
- 大数据集（>50万SNP）：16-32g
- 超大数据集（>100万SNP）：64g+

---

## 🚀 完整工作流程示例

```bash
# 进入工作目录
cd /path/to/your/project

# 运行GEC分析
biopytools gwas-gec \
  -i gwas_domain_p_file.txt \
  -r input.vcf.gz \
  -o gec_output \
  -t 64 \
  -m 100g

# 查看结果
cat gec_output/gwas_gec_summary.txt
```

**输出示例：**
```
总有效检验数: 456,789.12
校正后显著性阈值: 1.09e-07
```

---

## ❓ 常见问题

### Q1: KGGSee只支持VCF格式吗？

**A:** 是的。KGGSee（GEC算法实现）只支持VCF格式的参考文件。需要提供VCF格式的基因型数据。

### Q2: 内存不足怎么办？

**A:**
- 减少 `-t` 线程数
- 增加交换空间
- 或者按染色体分别分析

### Q3: 计算需要多长时间？

**A:**
- 小数据集（<50万SNP）：30分钟-1小时
- 大数据集（>50万SNP）：1-3小时
- 具体取决于数据量和硬件配置

### Q4: 可以使用压缩的VCF文件吗？

**A:** 可以。KGGSee支持 `.vcf.gz` 格式的压缩文件。

### Q5: P值文件格式有什么要求？

**A:**
- 必须是制表符分隔的文本文件
- 至少包含三列：染色体、位置、P值
- 表头必须包含列名（可自定义）
- 建议先排序：按染色体和位置排序

---

## 📚 更新日志

| 版本|Version | 日期|Date | 更新内容|Updates |
|---------|------------|----------|----------|
| 2.2.0 | 2026-01-14 | **简化设计|Simplified Design** |
| | | - ✅ 移除不必要的PLINK转换|Removed unnecessary PLINK conversion |
| | | - ✅ KGGSee仅支持VCF格式|KGGSee only supports VCF format |
| | | - ✅ 简化参数和配置|Simplified parameters and configuration |
| | | - ✅ 更清晰的错误提示|Clearer error messages |
| 2.1.0 | 2026-01-14 | **自动转换功能|Auto-conversion Feature** |
| | | - VCF自动转PLINK binary|Auto-convert VCF to PLINK binary |
| | | - 智能格式检测|Smart format detection |
| 2.0.0 | 2026-01-14 | **重大更新|Major Update** |
| | | - 重构代码，简化流程|Refactored code, simplified pipeline |
| 1.0.0 | 2026-01-14 | 初始版本|Initial release |

---

## 🎓 参考文献

Li, M., et al. (2021). GEC: a fast and efficient tool for genome-wide correction for multiple testing. *Bioinformatics*, 37(18), 2763-2769.

---

## 📧 技术支持

遇到问题？检查：
1. P值文件和参考文件的染色体格式是否一致
2. 内存分配是否足够
3. 文件路径是否正确

需要帮助？提供以下信息：
- 错误日志
- 输入文件格式
- 命令行参数
