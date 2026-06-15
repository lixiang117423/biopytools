# 🌳 VCF系统发育分析模块

**基于VCF变异数据的快速系统发育树构建工具 | Fast Phylogenetic Tree Construction Tool Based on VCF Variant Data**

## 📖 功能概述 | Overview

VCF系统发育分析模块是一个专业的生物信息学工具，用于从VCF格式的SNP数据直接构建系统发育树。集成了VCF2Dis遗传距离计算和邻接法(NJ)树构建算法，支持完整的从原始变异数据到高质量系统发育树的分析流程，适用于群体遗传学、进化生物学和系统发育学研究。

## ✨ 主要特性 | Key Features

- **🔄 完整流程**: VCF文件→遗传距离矩阵→系统发育树的一体化分析
- **⚡ 快速算法**: 基于邻接法(NJ)的高效树构建算法，适合大样本量分析
- **🎯 灵活输入**: 支持VCF文件或已有距离矩阵作为输入
- **📊 质量评估**: 内置结果验证和树拓扑质量检查
- **📝 详细报告**: 生成完整的分析报告和统计信息
- **🛠️ 工具集成**: 无缝集成VCF2Dis遗传距离计算工具
- **📁 灵活输出**: 支持自定义输出路径和文件命名
- **🔍 跳步功能**: 支持跳过距离计算步骤，直接从矩阵构建树

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 从VCF文件构建系统发育树
biopytools vcf-nj-tree -i population_snps.vcf -o phylo_results

# 从已有距离矩阵构建树
biopytools vcf-nj-tree --distance-matrix genetic_dist.txt -o tree_only

# 完整路径指定
biopytools vcf-nj-tree \
    -i /data/cohort_variants.vcf \
    -o /results/phylogeny \
    -t /results/cohort_tree.nwk
```

### 高级用法 | Advanced Usage

```bash
# 大规模群体分析
biopytools vcf-nj-tree \
    -i large_population.vcf.gz \
    -o comprehensive_analysis \
    -w ./phylo_working_dir/ \
    -t final_tree.newick

# 使用自定义工具路径
biopytools vcf-nj-tree \
    -i wild_species.vcf \
    -o wild_phylogeny \
    --vcf2dis-path /opt/local/VCF2Dis \
    --working-dir /scratch/phylo_work/

# 分步分析（从已有距离矩阵）
biopytools vcf-nj-tree \
    --distance-matrix precomputed_dist.txt \
    -o final_analysis \
    --skip-vcf2dis
```

## 📋 命令行参数 | Command Line Parameters

### 输入文件参数 | Input File Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | `None` | 📁 输入VCF文件路径 |
| `-d, --distance-matrix` | `None` | 📊 已有距离矩阵文件路径 |

### 输出配置参数 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `phylo_analysis` | 📝 输出文件前缀 |
| `-t, --tree-output` | `None` | 🌳 系统发育树输出文件路径 |

### 工具路径参数 | Tool Path Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--vcf2dis-path` | `VCF2Dis` | 🛠️ VCF2Dis程序路径 |
| `-w, --working-dir` | `.` | 📂 工作目录路径 |

### 行为控制参数 | Behavior Control

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--skip-vcf2dis` | `False` | ⏭️ 跳过VCF2Dis步骤，直接从距离矩阵构建树 |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **0** | 🔍 依赖检查 | 验证VCF2Dis工具和系统环境 |
| **1** | 📊 距离计算 | 从VCF计算样本间遗传距离矩阵 |
| **2** | 🌳 树构建 | 使用邻接法构建系统发育树 |
| **3** | ✅ 结果验证 | 检查树结构和距离矩阵质量 |
| **4** | 📝 报告生成 | 生成详细的分析报告和统计信息 |

## 📁 输入文件格式 | Input File Formats

### VCF变异文件 | VCF Variant File

标准VCF格式的群体变异数据：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1  Sample2  Sample3  Sample4
chr1    10177   .       A       AC      100     PASS    DP=20   GT      0/1     1/1     0/0     0/1
chr1    10352   .       T       TA      90      PASS    DP=18   GT      0/1     0/0     0/1     1/1
chr2    54678   .       G       GA      85      PASS    DP=22   GT      1/1     0/1     0/0     0/1
```

**文件要求**:
- 标准VCF 4.0+格式
- 包含多个样本的基因型数据
- 建议预先进行质量控制（MAF、缺失率过滤）
- 样本名称不含特殊字符
- 基因型字段完整（GT格式）

**数据质量建议**:
- SNP位点MAF > 0.05
- 样本缺失率 < 20%
- 位点缺失率 < 10%
- 高质量基因型调用

### 遗传距离矩阵文件 | Genetic Distance Matrix File

方阵格式的样本间距离数据：

```
        Sample1  Sample2  Sample3  Sample4
Sample1  0.0000   0.1234   0.2345   0.3456
Sample2  0.1234   0.0000   0.2567   0.3678
Sample3  0.2345   0.2567   0.0000   0.2789
Sample4  0.3456   0.3678   0.2789   0.0000
```

**文件要求**:
- 方阵格式（行数=列数）
- 第一行和第一列为样本名称
- 距离值对称，对角线为0
- 制表符或空格分隔
- 距离值非负实数

## 📂 输出文件说明 | Output Files Description

### 系统发育树文件 | Phylogenetic Tree File

Newick格式的系统发育树：

```newick
((Sample1:0.056,Sample2:0.067):0.034,Sample3:0.123,Sample4:0.145);
```

**文件特点**:
- 标准Newick格式，兼容大部分系统发育软件
- 枝长代表进化距离
- 括号表示拓扑结构
- 分号结尾

### 遗传距离矩阵 | Genetic Distance Matrix

从VCF计算的距离矩阵：

```
        Sample1  Sample2  Sample3  Sample4  Sample5
Sample1  0.0000   0.0234   0.0345   0.0456   0.0567
Sample2  0.0234   0.0000   0.0467   0.0578   0.0689
Sample3  0.0345   0.0467   0.0000   0.0689   0.0791
Sample4  0.0456   0.0578   0.0689   0.0000   0.0812
Sample5  0.0567   0.0689   0.0791   0.0812   0.0000
```

### 分析报告文件 | Analysis Report File

详细的分析统计信息：

```
========================================
VCF系统发育分析报告 | VCF Phylogenetic Analysis Report
========================================

分析时间 | Analysis Time: 2024-12-17 15:30:45
输入文件 | Input File: population_snps.vcf

样本信息 | Sample Information:
  - 总样本数 | Total Samples: 25
  - 样本名称 | Sample Names: [Sample1, Sample2, ..., Sample25]

SNP信息 | SNP Information:
  - 总SNP数 | Total SNPs: 15,678
  - 过滤后SNPs | Filtered SNPs: 12,345
  - 使用率 | Usage Rate: 78.8%

距离矩阵统计 | Distance Matrix Statistics:
  - 最小距离 | Minimum Distance: 0.0012
  - 最大距离 | Maximum Distance: 0.2345
  - 平均距离 | Mean Distance: 0.0891
  - 标准差 | Standard Deviation: 0.0456

系统发育树统计 | Phylogenetic Tree Statistics:
  - 树节点数 | Tree Nodes: 49
  - 树枝数 | Tree Branches: 48
  - 最大枝长 | Maximum Branch Length: 0.1456
  - 平均枝长 | Mean Branch Length: 0.0345

质量检查 | Quality Check:
  - 距离矩阵对称性 | Distance Matrix Symmetry: ✅ PASS
  - 对角线有效性 | Diagonal Validity: ✅ PASS
  - 树拓扑有效性 | Tree Topology Validity: ✅ PASS
```

### 日志文件 | Log File

详细的分析过程记录：

```
2024-12-17 15:30:00 INFO  开始VCF系统发育分析 | Starting VCF phylogenetic analysis
2024-12-17 15:30:01 INFO  检查依赖软件 | Checking dependencies
2024-12-17 15:30:02 INFO  VCF2Dis工具路径 | VCF2Dis tool path: /usr/local/bin/VCF2Dis
2024-12-17 15:30:05 INFO  解析VCF文件 | Parsing VCF file
2024-12-17 15:30:06 INFO  提取样本信息 | Extracting sample information
2024-12-17 15:30:07 INFO  发现25个样本 | Found 25 samples
2024-12-17 15:30:10 INFO  计算遗传距离矩阵 | Calculating genetic distance matrix
2024-12-17 15:32:45 INFO  距离矩阵计算完成 | Distance matrix calculation completed
2024-12-17 15:32:46 INFO  构建系统发育树 | Building phylogenetic tree
2024-12-17 15:32:48 INFO  邻接法算法执行 | Neighbor-Joining algorithm execution
2024-12-17 15:32:50 INFO  树构建完成 | Tree construction completed
2024-12-17 15:32:51 INFO  验证结果 | Validating results
2024-12-17 15:32:52 INFO  生成分析报告 | Generating analysis report
2024-12-17 15:32:55 INFO  分析完成 | Analysis completed successfully
```

## 🔧 使用示例 | Usage Examples

### 人类群体遗传学分析

```bash
# 人类群体SNP系统发育分析
biopytools vcf-nj-tree \
    -i human_population.vcf.gz \
    -o human_phylogeny \
    -w ./human_analysis/ \
    -t human_evolution_tree.nwk

# 输出文件：
# - human_phylogeny.dis (距离矩阵)
# - human_evolution_tree.nwk (系统发育树)
# - human_phylogeny_report.txt (分析报告)
# - human_phylogeny.log (日志文件)
```

### 野生植物群体研究

```bash
# 野生植物遗传多样性分析
biopytools vcf-nj-tree \
    -i wild_plant_snps.vcf \
    -o wild_plant_diversity \
    --vcf2dis-path /opt/bioinfo/VCF2Dis \
    --working-dir /scratch/plant_work/

# 后续可视化分析
Rscript plot_tree.R wild_plant_diversity.nwk
```

### 微生物菌株分类

```bash
# 细菌菌株系统发育分类
biopytools vcf-nj-tree \
    -i bacterial_strains.vcf \
    -o bacterial_classification \
    -t bacterial_tree.newick

# 基于已有距离矩阵
biopytools vcf-nj-tree \
    --distance-matrix bacterial_distance.txt \
    -o bacterial_classification \
    --skip-vcf2dis
```

### 大规模群体分析

```bash
# 千人级别大规模分析
biopytools vcf-nj-tree \
    -i large_cohort.vcf.gz \
    -o large_scale_phylogeny \
    -w ./large_analysis/ \
    --threads 88
```

## ⚡ 性能优化 | Performance Optimization

### 硬件配置建议

| 样本数 | SNP数 | 内存需求 | CPU时间 | 存储空间 |
|--------|-------|----------|---------|----------|
| 50 | 10K | <4GB | <5分钟 | <100MB |
| 200 | 50K | 8-16GB | 15-30分钟 | 500MB-1GB |
| 500 | 100K | 32-64GB | 1-2小时 | 2-5GB |
| 1000+ | 100K+ | 64GB+ | 2-4小时 | 5-10GB+ |

### 优化策略

```bash
# 1. 预处理VCF文件
bcftools view -i 'MAF>0.05 && INFO/DP>10' raw.vcf -o filtered.vcf

# 2. 使用高速存储
export TMPDIR=/ssd/temp
biopytools vcf-nj-tree -i data.vcf -o results -w /ssd/analysis/

# 3. 并行处理（如果支持）
# VCF2Dis工具本身的并行设置

# 4. 分段分析（超大数据集）
# 按染色体分别分析，再合并结果
```

### 内存管理

```bash
# 监控内存使用
watch -n 1 'free -h && ps aux | grep VCF2Dis'

# 设置内存限制（可选）
ulimit -v 32768  # 限制32GB虚拟内存
biopytools vcf-nj-tree -i large_data.vcf -o output
```

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**VCF2Dis工具找不到**
```bash
# 检查工具是否安装
which VCF2Dis
VCF2Dis -h

# 指定工具路径
biopytools vcf-nj-tree -i data.vcf -o output --vcf2dis-path /custom/path/VCF2Dis

# 添加到PATH
export PATH=$PATH:/opt/bioinfo/bin
```

**VCF文件格式错误**
```bash
# 验证VCF格式
bcftools view -H input.vcf | head -5
bcftools query -l input.vcf | wc -l

# 检查基因型字段
bcftools view -H input.vcf | cut -f10 | head -5

# 标准化VCF格式
bcftools norm -f reference.fa input.vcf -o normalized.vcf
```

**内存不足**
```bash
# 1. 预过滤数据
bcftools view -i 'MAC>1 && DP>5' large.vcf -o reduced.vcf

# 2. 减少样本数量
# 只保留关键样本进行初步分析

# 3. 使用swap空间
sudo swapon /swapfile
```

**距离矩阵异常**
```bash
# 检查矩阵对称性
python -c "
import numpy as np
dist = np.loadtxt('distance.txt')
print('对称性检查:', np.allclose(dist, dist.T))
print('对角线检查:', np.all(np.diag(dist) == 0))
"
```

### 性能调试 | Performance Debugging

```bash
# 监控系统资源
htop          # CPU和内存
iotop         # I/O性能
nvidia-smi    # GPU使用（如适用）

# 测试不同数据量
for n in 10 50 100 200; do
    echo "测试样本数: $n"
    biopytools vcf-nj-tree -i test_${n}.vcf -o test_${n}_output
done

# 时间分析
time biopytools vcf-nj-tree -i benchmark.vcf -o benchmark_output
```

### 结果验证 | Result Validation

```bash
# 验证Newick格式
cat output.nwk

# 使用R包ape验证
Rscript -e "library(ape); tree <- read.tree('output.nwk'); plot(tree); nodelabels()"

# 计算树统计信息
Rscript -e "library(ape); tree <- read.tree('output.nwk'); print(tree); print(Ntip(tree)); print(Nnode(tree))"

# 距离矩阵一致性检查
# 比较VCF2Dis输出和手动计算的距离
```

## 📊 质量控制建议 | Quality Control Recommendations

### 输入数据质量控制

```bash
# 1. 基本质量过滤
bcftools filter -e 'QUAL<30 || DP<10' input.vcf -o filtered1.vcf

# 2. MAF过滤
bcftools view -i 'MAF>0.05' filtered1.vcf -o filtered2.vcf

# 3. 缺失率过滤
vcftools --vcf filtered2.vcf --max-missing 0.8 --recode --out filtered3

# 4. 样本质量检查
vcftools --vcf filtered3.vcf --missing-indv --out missing_report
# 移除缺失率高的样本

# 5. 连锁不平衡修剪（可选）
plink --vcf filtered3.vcf --indep-pairwise 50 5 0.2 --out pruning
plink --vcf filtered3.vcf --extract pruning.prune.in --make-bed --out final_data
```

### 分析参数优化

```bash
# 保守设置（高质量数据）
biopytools vcf-nj-tree -i high_quality.vcf -o conservative_result

# 宽松设置（数据稀疏时）
# 预先进行较宽松的过滤，然后分析

# 中等设置（常规推荐）
biopytools vcf-nj-tree -i standard_data.vcf -o standard_result
```

### 结果验证策略

```bash
# 1. Bootstrap验证（如果有）
# 多次随机采样构建树，检查拓扑一致性

# 2. 不同方法比较
# 与最大似然法、贝叶斯法结果比较

# 3. 距离方法比较
# 尝试不同的距离计算方法
# 如：欧氏距离、曼哈顿距离等

# 4. 外群检验
# 包含已知的外群样本，检验根位置是否合理
```

## 🔗 相关文档 | Related Documentation

- [VCF2Dis工具文档](https://github.com/BGI-shenzhen/VCF2Dis)
- [邻接法算法原理](https://en.wikipedia.org/wiki/Neighbour_joining)
- [Newick格式说明](https://en.wikipedia.org/wiki/Newick_format)
- [系统发育分析最佳实践](https://academic.oup.com/mbe/article/35/5/1303/4997218)
- [biopytools VCF转换工具](vcf2phylip.md)
- [IQ-TREE系统发育软件](https://iqtree.org/)
- [RAxML-NG软件手册](https://github.com/amkozlov/raxml-ng)
- [FigTree树可视化工具](http://tree.bio.ed.ac.uk/software/figtree/)

## 📄 许可证 | License

本模块遵循MIT许可证。详细信息请参见LICENSE文件。

## 🤝 贡献指南 | Contributing

欢迎提交Issue和Pull Request来改进本模块。

## 📞 技术支持 | Support

如有技术问题，请联系：
- 邮箱: yzwl_lixg@outlook.com
- 项目地址: https://github.com/your-org/biopytools

---

**最后更新**: 2024年12月17日
**版本**: 2.0.0
**作者**: biopytools开发团队