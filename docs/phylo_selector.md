# 系统发育树样品选择工具

**从系统发育树中选择代表性样品，支持分层抽样和最大间距优化 | Select Representative Samples from Phylogenetic Trees with Stratified Sampling and Maximum Spacing Optimization**

## 功能概述 | Overview

系统发育树样品选择工具是一个专业的生物信息学工具，用于从基于枝长的系统发育树中选择最具代表性的样品子集。工具实现了两种策略：简单选择（最大间距贪心算法）和分层选择（结合分组信息的最大间距算法），适用于三代基因组测序、群体遗传学和进化生物学研究。

## 主要特性 | Key Features

- **两种选择策略**:
  - 简单选择：使用最大间距贪心算法
  - 分层选择：结合分组信息，确保各组代表性
- **智能样品分配**：按比例分配各组样品数，确保平衡
- **间距优化**：最大化样品间最小距离，避免选择相近样品
- **交集处理**：自动处理进化树与分组表的交集
- **详细报告**：生成完整的统计报告和可视化文件
- **灵活输出**：支持多种输出格式（TXT、CSV、报告、可视化）

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 简单选择（不分组）
biopytools phylo-selector -i tree.nwk -o selected_samples

# 分层选择（带分组信息）
biopytools phylo-selector -i tree.nwk -g groups.txt -o selected_samples -n 100

# 自定义参数
biopytools phylo-selector \
    -i tree.nwk \
    -g groups.txt \
    -o results/selected_samples \
    -n 150 \
    --min-distance 0.0005 \
    --min-samples-per-group 2
```

## 命令行参数 | Command Line Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 |
|------|------|
| `-i, --newick-file` | Newick格式系统发育树文件 |
| `-o, --output-prefix` | 输出文件前缀 |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-g, --group-file` | `None` | 样品分组表文件 |
| `-n, --n-samples` | `150` | 选择样品总数 |
| `--min-distance` | `0.001` | 最小枝长间距阈值 |
| `--min-samples-per-group` | `1` | 每组最小样品数 |
| `--max-ratio-per-group` | `0.8` | 每组最大样品比例 |

### 输出控制 | Output Control

| 参数 | 描述 |
|------|------|
| `--no-report` | 不生成详细报告 |
| `--no-csv` | 不生成CSV文件 |
| `--no-visualization` | 不生成间距可视化 |

## 输入文件格式 | Input File Formats

### Newick系统发育树文件 | Newick Phylogenetic Tree File

标准Newick格式的系统发育树文件：

```newick
((Sample1:0.005,Sample2:0.006):0.003,Sample3:0.012,Sample4:0.015);
```

**文件要求**:
- 标准Newick格式
- 样品名使用字母、数字、下划线、连字符
- 枝长使用科学计数法或小数格式

### 样品分组表文件 | Sample Group Table File

两列格式：样品名、分组名（支持CSV、TSV、空格分隔）：

```csv
Sample_ID    Group
W580A        Group1
GDW012A      Group1
W215A        Group2
DB198A       Group3
```

**文件要求**:
- 第一列：样品名称（必须与Newick树中的名称一致）
- 第二列：分组名称
- 支持制表符、逗号或空格分隔
- 跳过以#开头的注释行

## 选择策略说明 | Selection Strategy Description

### 策略1：简单选择（最大间距贪心算法）

**适用场景**：无分组信息，仅需基于枝长选择代表性样品

**算法步骤**：
1. 选择枝长最小和最大的样品作为锚点
2. 迭代选择与所有已选样品距离最远的样品
3. 最大化样品间的最小距离
4. 确保样品均匀分布在整个枝长范围

**优势**：
- 避免选择枝长过于相近的样品
- 样品在枝长范围内均匀分布
- 最大化覆盖整个变异范围

### 策略2：分层选择（按比例分配 + 组内最大间距）

**适用场景**：有分组信息（如地理来源、品种、表型），需确保各组代表性

**算法步骤**：
1. 求进化树与分组表的样品交集
2. 统计各组样品数，按比例分配应选样品数
3. 确保每组至少选择`min_samples_per_group`个样品
4. 每组最多不超过该组样品数的`max_ratio_per_group`比例
5. 在每个组内使用最大间距算法选择样品
6. 合并所有组的样品

**分配策略**：
- **第一轮**：确保每组至少有最小样品数
- **第二轮**：按比例分配剩余样品
- **第三轮**：处理剩余未分配的样品

**优势**：
- 保证每个分组都有代表
- 避免某些组被过度代表或遗漏
- 组内样品间距优化
- 结合先验知识和数据驱动

## 输出文件说明 | Output Files Description

### 1. 样品列表文件 | Sample List File

**文件名**: `{prefix}.txt`

**内容**：
```
# 系统发育树样品选择结果
# 策略: 分层选择
# 总样品数: 806
# 选择数量: 150

Sample1
Sample2
Sample3
...
```

### 2. 详细报告 | Detailed Report

**文件名**: `{prefix}_report.txt`

**内容**：
- 项目概况（样品数、选择比例）
- 选择策略说明
- 样品间距统计（最小、最大、平均、中位数）
- 间距质量分析
- 已选样品列表（按枝长排序，含间距信息）
- 统计对比（全部样品 vs 选中样品）

### 3. CSV文件 | CSV File

**文件名**: `{prefix}.csv`

**内容**：
```csv
序号,样品编号,枝长,与前一样品间距
1,Sample1,0.000144,-
2,Sample2,0.001732,0.001588
3,Sample3,0.002340,0.000608
...
```

### 4. 间距可视化 | Spacing Visualization

**文件名**: `{prefix}_visualization.txt`

**内容**：间距分布直方图
```
间距分布直方图 (每个 * 代表1个样品对):
<0.0005         [  0]
0.0005-0.001    [ 83] *******************
0.001-0.002     [ 62] ***************
...
```

### 5. 日志文件 | Log File

**文件名**: `{prefix}.log`

**内容**：详细的分析过程记录

## 使用示例 | Usage Examples

### 示例1：三代基因组测序样品选择

```bash
# 从806个样品中选择150个代表
biopytools phylo-selector \
    -i phylogenetic_tree.nwk \
    -o results/selected_samples \
    -n 150

# 输出文件：
# - results/selected_samples.txt (样品列表)
# - results/selected_samples_report.txt (详细报告)
# - results/selected_samples.csv (CSV文件)
# - results/selected_samples_visualization.txt (可视化)
# - results/selected_samples.log (日志文件)
```

### 示例2：考虑地理来源的分层选择

```bash
# 按地理分组选择样品
biopytools phylo-selector \
    -i population_tree.nwk \
    -g geographic_groups.txt \
    -o results/geographic_representatives \
    -n 200 \
    --min-samples-per-group 5

# geographic_groups.txt内容：
# Sample1    North
# Sample2    North
# Sample3    South
# Sample4    South
# Sample5    East
# ...
```

### 示例3：品种选择（平衡各品种代表）

```bash
# 选择不同品种的代表
biopytools phylo-selector \
    -i cultivar_tree.nwk \
    -g cultivar_groups.txt \
    -o results/cultivar_selection \
    -n 100 \
    --min-samples-per-group 3 \
    --max-ratio-per-group 0.7

# 确保每个品种至少3个样品
# 任何品种不超过总样品数的70%
```

### 示例4：快速测试（仅生成样品列表）

```bash
# 不生成报告和可视化，加快速度
biopytools phylo-selector \
    -i test_tree.nwk \
    -o test_selected \
    -n 50 \
    --no-report \
    --no-csv \
    --no-visualization
```

## 交集处理说明 | Intersection Processing

工具会自动处理进化树和分组表的交集：

**示例场景**：
- 进化树包含：806个样品
- 分组表包含：1000个样品
- 交集（有效样品）：750个样品

**日志输出**：
```
样品交集统计|Sample intersection statistics:
  进化树独有|Tree only: 56 个样品|samples
  分组表独有|Group table only: 250 个样品|samples
  交集（使用）|Intersection (used): 750 个样品|samples
```

**注意**：只选择在两个集合中都存在的样品

## 性能建议 | Performance Recommendations

### 数据规模建议

| 样品数 | 推荐选择数 | 最小间距建议 | 运行时间 |
|--------|-----------|-------------|----------|
| 100-500 | 50-100 | 0.001 | <1分钟 |
| 500-1000 | 100-200 | 0.001 | 1-3分钟 |
| 1000-5000 | 200-500 | 0.0005 | 3-10分钟 |
| 5000+ | 500+ | 0.0005 | 10分钟+ |

### 参数调优建议

**提高代表性**：
- 降低`--min-distance`（如0.0005）
- 增加`--n-samples`
- 使用分组信息

**提高多样性**：
- 增加`--min-samples-per-group`
- 降低`--max-ratio-per-group`

**加快运行**：
- 使用`--no-report`、`--no-csv`、`--no-visualization`
- 减少`--n-samples`

## 常见问题 | FAQ

### Q1: 如何确定选择多少个样品？

**A**: 根据研究目的和预算：
- 初步探索：选择10-20%
- 深入研究：选择20-30%
- 全覆盖：如预算允许，可考虑50%+

### Q2: 分组文件格式错误怎么办？

**A**: 检查：
1. 文件编码是否为UTF-8
2. 分隔符是否正确（自动检测制表符、逗号、空格）
3. 样品名是否与Newick树中完全一致
4. 是否有注释行（以#开头）

### Q3: 为什么某些组没有样品被选中？

**A**: 可能原因：
1. 该组样品不在进化树中（查看交集统计）
2. 该组样品太少，被分配数量为0
3. 增加`--min-samples-per-group`可避免

### Q4: 如何验证选择结果？

**A**: 方法：
1. 查看`_report.txt`中的间距统计
2. 查看`_visualization.txt`中的直方图
3. 使用`_csv`文件在Excel中分析
4. 在系统发育树中标记选中样品

### Q5: 能否同时使用多个分组文件？

**A**: 当前版本不支持。建议：
1. 预先合并分组信息
2. 或分别运行多次选择
3. 或创建组合分组列

## 最佳实践 | Best Practices

1. **准备阶段**：
   - 确保Newick树格式正确
   - 样品名在两个文件中完全一致
   - 预览分组表，检查分组合理性

2. **参数选择**：
   - 从默认参数开始
   - 根据结果调整`--n-samples`
   - 根据样品分布调整`--min-distance`

3. **结果验证**：
   - 检查交集统计，确保没有遗漏重要样品
   - 查看间距统计，评估样品多样性
   - 在系统发育树中可视化选中样品

4. **结果使用**：
   - 保存`.txt`文件作为送测清单
   - 保存`_report.txt`用于项目存档
   - 保存`.csv`文件用于进一步分析

## 相关工具 | Related Tools

- [VCF系统发育分析](vcf_phylo.md)：从VCF构建系统发育树
- [MSA可视化](msaviz.md)：多序列比对可视化
- [IQ-TREE](iqtree.md)：最大似然系统发育树
- [RAxML](raxml.md)：RAxML系统发育分析

## 技术支持 | Support

如有问题或建议，请联系：
- 邮箱: yzwl_lixg@outlook.com
- 项目地址: https://github.com/your-org/biopytools

---

**版本**: 1.0.0
**作者**: Xiang LI
**最后更新**: 2026-01-30
