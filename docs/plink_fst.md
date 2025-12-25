# 🧊 PLINK Fst 计算模块

**使用PLINK计算群体间Fst值的专业工具 | Professional Tool for Calculating Fst Using PLINK**

## 📖 功能概述 | Overview

PLINK Fst计算模块是一个专门用于使用PLINK计算群体间遗传分化指数（Fst）的工具。支持全基因组、特定区间和滑动窗口三种计算模式，提供完整的质控过滤和结果分析功能。

**自动两两群体比较**：当有3个或更多群体时，自动计算所有可能的群体对之间的Fst值。

## ✨ 主要特性 | Key Features

- **🎯 三种计算模式**: 全基因组、区间、滑动窗口
- **👥 自动两两比较**: 多群体时自动计算所有群体对的Fst值
- **📋 简化的样本文件**: 只需两列（样品名 分组）
- **🔬 质控过滤**: 支持MAF、GENO、MIND、HWE过滤
- **📊 自动结果解析**: 统计摘要和高Fst SNP提取
- **🪟 滑动窗口**: 支持自定义窗口大小和步长
- **🛡️ 完整验证**: 输入文件和参数验证
- **📝 标准日志**: 遵循项目日志规范

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 全基因组Fst计算
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./output
```

### 高级用法 | Advanced Usage

```bash
# 指定染色体
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./output \
    --mode region \
    --chr 1

# 滑窗计算
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./output \
    --mode window \
    --window-size 1000000 \
    --step-size 500000
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入文件前缀 (.bed/.bim/.fam) 或 VCF文件 (.vcf/.vcf.gz) | `-i mydata` 或 `-i data.vcf` |
| `-s, --sample` | 样本分组文件 | `-s samples.txt` |
| `-o, --output` | 输出目录 | `-o ./output` |

### 样本文件格式 | Sample File Format

**输入格式**（samples.txt）：
```
Sample1  POP1
Sample2  POP1
Sample3  POP2
Sample4  POP2
```

**脚本自动转换为PLINK格式**：
```
Sample1  Sample1  POP1
Sample2  Sample2  POP1
Sample3  Sample3  POP2
Sample4  Sample4  POP2
```

### 计算模式 | Calculation Modes

| 模式 | 说明 | 适用场景 |
|------|------|----------|
| `global` | 全基因组Fst | 整体分化评估 |
| `region` | 区间Fst | 特定区域分析 |
| `window` | 滑窗Fst | 局部分化模式 |

### 区间模式参数 | Region Mode Parameters

| 参数 | 描述 |
|------|------|
| `--chr` | 染色体 |
| `--from-bp` | 起始位置（bp） |
| `--to-bp` | 结束位置（bp） |

### 滑窗模式参数 | Window Mode Parameters

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--window-size` | 窗口大小（bp） | 必需 |
| `--step-size` | 步长（bp） | 必需 |

### 质控参数 | Quality Control Parameters

| 参数 | 描述 | 推荐值 |
|------|------|--------|
| `--maf` | 最小等位基因频率 | 0.05 |
| `--geno` | 基因型缺失率 | 0.1 |
| `--mind` | 个体缺失率 | 0.1 |
| `--hwe` | HWE p值阈值 | 0.0001 |

## 💡 使用示例 | Usage Examples

### 示例1：全基因组Fst计算（PLINK格式）

```bash
biopytools plink-fst \
    -i /data/population \
    -s population_groups.txt \
    -o ./fst_results
```

**输出**：
- `plink_fst.fst` - 所有SNP的Fst值
- `plink_fst_high_fst.txt` - 高分化SNP（Fst > 0.2）

### 示例2：使用VCF文件直接计算

```bash
# 从VCF文件直接计算Fst（无需预先转换）
biopytools plink-fst \
    -i /data/variants.vcf.gz \
    -s population_groups.txt \
    -o ./fst_results

# 或使用未压缩的VCF
biopytools plink-fst \
    -i /data/variants.vcf \
    -s population_groups.txt \
    -o ./fst_results
```

**说明**：
- VCF文件（.vcf或.vcf.gz）会被自动识别
- 程序内部使用PLINK的`--vcf`参数处理
- 对于滑窗模式，会自动转换为PLINK格式以生成窗口

### 示例3：指定染色体计算

```bash
# 计算1号染色体的Fst
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./chr1_fst \
    --mode region \
    --chr 1
```

### 示例4：特定区间Fst

```bash
# 1号染色体1-2Mb区间
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./region_fst \
    --mode region \
    --chr 1 \
    --from-bp 1000000 \
    --to-bp 2000000
```

### 示例5：滑窗Fst分析

```bash
# 1Mb窗口，500K步长
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./window_fst \
    --mode window \
    --window-size 1000000 \
    --step-size 500000
```

### 示例6：质控过滤

```bash
biopytools plink-fst \
    -i mydata \
    -s samples.txt \
    -o ./filtered_fst \
    --maf 0.05 \
    --geno 0.1 \
    --mind 0.1
```

## 📁 输出文件 | Output Files

### 全基因组/区间模式

```
output_dir/
├── plink_fst.fst                       # PLINK原始结果（整体Fst）
├── plink_fst_summary.txt               # 整体Fst统计摘要
├── plink_fst_high_fst.txt              # 高分化SNP列表（整体）
├── plink_fst_pairwise_summary.txt      # 两两比较Fst汇总表
├── pairwise_POP1_vs_POP2.fst           # POP1 vs POP2的Fst值
├── pairwise_POP1_vs_POP3.fst           # POP1 vs POP3的Fst值
├── pairwise_POP2_vs_POP3.fst           # POP2 vs POP3的Fst值
└── samples.cluster                     # 转换后的样本文件
```

### 滑窗模式

```
output_dir/
├── plink_fst_windows.fst              # 所有窗口的Fst值
├── plink_fst_windows_summary.txt
├── plink_fst_windows_high_fst.txt
├── windows.bed                         # 窗口BED文件
└── samples.cluster
```

**注意**：滑窗模式不自动执行两两比较（计算量较大），如需两两比较请使用global或region模式。

### 结果文件格式

**.fst文件**：
```
CHR     SNP      A1  A2  FST
1       rs1001   A   G   0.125
1       rs1002   C   T   0.234
```

**统计摘要**：
- 总SNP数
- 平均/中位数Fst
- 高分化SNP数量和比例

**两两比较汇总表**（plink_fst_pairwise_summary.txt）：
```
Population1  Population2  Mean_Fst  Median_Fst  Min_Fst  Max_Fst  Std_Fst  Num_SNPs  High_Fst_SNPs  High_Fst_Percent
POP1         POP2         0.0523    0.0487      0.0000   0.2456   0.0321   10000     245            2.45
POP1         POP3         0.0891    0.0854      0.0012   0.3567   0.0456   10000     567            5.67
POP2         POP3         0.1024    0.0989      0.0023   0.4123   0.0512   10000     723            7.23
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **PLINK 1.9+**
- **Python** (版本 3.7+)
- **Python包**:
  - `click`
  - `pandas`

## ⚠️ 注意事项 | Important Notes

1. **样本文件格式**: 必须是两列（样品名 分组）
2. **至少两个群体**: Fst计算至少需要两个群体
3. **样本名匹配**: 样本名必须与输入文件中的样品名完全匹配（VCF的列名或.fam的IID）
4. **窗口设置**: 窗口大小应大于步长
5. **质控建议**: 建议在Fst计算前进行质控过滤
6. **VCF输入**: 支持直接使用VCF文件（.vcf或.vcf.gz），滑窗模式会自动转换为PLINK格式
7. **两两比较**:
   - 当有3个或更多群体时，自动进行所有可能的两两比较
   - 2个群体时，整体Fst即为该群体对的Fst
   - 滑窗模式不自动执行两两比较（计算量大）
   - 两两比较结果保存在 `plink_fst_pairwise_summary.txt`

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: 样本文件格式错误**
```
错误 | Error: 跳过无效行
解决 | Fix: 确保每行有两列，用空格或Tab分隔
```

**Q: Fst计算失败**
```
检查 | Check: 
- 是否有至少两个群体
- 样本名是否匹配.fam文件
- 是否有足够的SNP
```

**Q: 滑窗模式很慢**
```
优化 | Optimization:
- 增大步长（如 window_size 的 50%）
- 使用质控过滤减少SNP数量
- 考虑只分析特定染色体
```

## 📚 应用场景 | Application Scenarios

1. **群体遗传学研究**: 评估不同群体间的遗传分化
2. **选择信号检测**: 寻找高Fst区域（可能的受选择位点）
3. **种群结构分析**: 结合PCA、Admixture等分析
4. **局部适应性**: 滑窗Fst识别局部适应区域
5. **种群保护**: 评估亚种群分化程度
6. **多群体比较**: 自动化两两比较，识别高度分化的群体对

## 📄 许可证 | License

MIT License

---

## 🔬 引用信息 | Citation

使用PLINK进行Fst计算，请引用：

```
Purcell S, Neale B, Todd-Brown K et al. PLINK: a tool set for whole-genome 
association and population-based linkage analyses. Am J Hum Genet. 2007 Sep;89(3):559-75.
doi: 10.1086/519795
```
