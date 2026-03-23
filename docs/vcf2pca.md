# VCF2PCA 主成分分析工具

**VCF文件主成分分析工具，支持VCF2PCACluster和PLINK双后端 | VCF Principal Component Analysis Tool with VCF2PCACluster and PLINK Backends**

## 功能概述 | Overview

VCF2PCA是一个功能强大的VCF文件主成分分析工具，支持两种分析后端：**VCF2PCACluster**（默认）和**PLINK**。它能够从VCF格式的基因型数据中提取主成分，帮助研究者识别样本群体结构、亲缘关系和种群分层，广泛应用于群体遗传学、进化生物学和生物医学研究。

## 主要特性 | Key Features

### 核心功能 | Core Features
- **双后端支持**: 集成VCF2PCACluster和PLINK两种分析引擎
- **多种聚类算法**: V2P后端支持K-means、DBSCAN、EM高斯聚类
- **质控过滤**: PLINK后端提供完整的质量控制流程（MAF、缺失率、HWE检验）
- **Kinship矩阵**: V2P后端可计算样本间亲缘关系矩阵
- **灵活配置**: 丰富的参数选项，满足不同分析需求

### 技术特点 | Technical Highlights
- **高效处理**: 多线程支持，快速处理大规模数据
- **内存优化**: V2P后端内存占用低，支持百万级SNP
- **标准化输出**: 结果与tassel、gapit、gcta完全一致
- **日志追踪**: 详细的运行日志，便于问题排查
- **可视化支持**: 支持PCA结果可视化（PLINK后端）

## 快速开始 | Quick Start

### 使用VCF2PCACluster后端（推荐）| Using VCF2PCACluster Backend (Recommended)

```bash
# 基本PCA分析
biopytools vcf2pca \
    -i variants.vcf \
    -o pca_results

# 启用聚类分析
biopytools vcf2pca \
    -i variants.vcf \
    -o pca_results \
    --cluster \
    --cluster-method kmeans \
    --cluster-k 3
```

### 使用PLINK后端 | Using PLINK Backend

```bash
# 带质控的PCA分析
biopytools vcf2pca \
    -i variants.vcf \
    -o pca_results \
    --backend plink \
    --maf 0.05 \
    --missing 0.1 \
    --plot

# 跳过质控
biopytools vcf2pca \
    -i variants.vcf \
    -o pca_results \
    --backend plink \
    --skip-qc
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入VCF文件路径（支持.vcf.gz）| `-i variants.vcf` |
| `-o, --output` | 输出目录路径 | `-o pca_output` |

### 后端选择 | Backend Selection

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-b, --backend` | `v2p` | 分析后端：`v2p` (VCF2PCACluster) 或 `plink` |

### PCA参数 | PCA Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-c, --components` | `10` | 主成分数量 |
| `-t, --threads` | `12` | 线程数 |

### 质控参数（仅PLINK后端）| Quality Control (PLINK Backend Only)

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--maf` | `0.05` | 最小等位基因频率阈值 |
| `--missing` | `0.1` | 最大缺失率阈值 |
| `--hwe` | `1e-6` | Hardy-Weinberg平衡p值阈值 |
| `--skip-qc` | `false` | 跳过质量控制过滤 |

### 聚类参数（仅V2P后端）| Clustering (V2P Backend Only)

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--cluster` | `false` | 启用聚类分析 |
| `--cluster-method` | `kmeans` | 聚类方法：`kmeans`、`dbscan`、`em` |
| `--cluster-k` | `3` | K-means聚类数量 |

### 可视化参数 | Visualization Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-P, --plot` | `false` | 生成PCA可视化图表（PLINK后端） |

### 工具路径 | Tool Paths

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--vcf2pca-path` | `~/software/VCF2PCACluster-1.42/bin/VCF2PCACluster` | VCF2PCACluster程序路径 |
| `--plink-path` | `plink` | PLINK程序路径 |
| `--bcftools-path` | `bcftools` | BCFtools程序路径 |

### 其他参数 | Other Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-s, --sample-info` | 样本信息文件路径 | `-s samples.txt` |

## 后端对比 | Backend Comparison

### VCF2PCACluster后端（V2P）| VCF2PCACluster Backend

**优点：**
- ✅ 一键完成，无需中间文件转换
- ✅ 速度快，多线程优化
- ✅ 内存效率高
- ✅ 内置多种聚类算法（K-means、DBSCAN、EM）
- ✅ 输出Kinship矩阵
- ✅ 结果与tassel/gapit/gcta一致

**适用场景：**
- 大规模数据集（百万级SNP）
- 需要聚类分析
- 需要Kinship矩阵
- 内存有限的环境

### PLINK后端

**优点：**
- ✅ 成熟稳定，行业标准
- ✅ 质量控制流程完善
- ✅ 生态系统兼容性好
- ✅ 支持Python可视化（matplotlib/seaborn）
- ✅ 输出结果详细

**适用场景：**
- 需要严格的质量控制
- 需要与PLINK生态系统兼容
- 需要Python可视化
- 小到中等规模数据集

## 输出结果 | Output Results

### V2P后端输出 | V2P Backend Output

```
pca_output/
├── vcf2pca.log                      # 运行日志
├── pca_explained_variance.txt       # 解释方差表
├── vcf2pca.PCA.eig                  # 特征值文件
├── vcf2pca.PCA.eigvec               # 特征向量文件
├── vcf2pca.Kinship                  # Kinship矩阵
└── vcf2pca.PCA.Cluster              # 聚类结果（如果启用）
```

### PLINK后端输出 | PLINK Backend Output

```
pca_output/
├── vcf2pca.log                      # 运行日志
├── pca.eigenval                     # 特征值
├── pca.eigenvec                     # 特征向量
├── pca_summary.txt                  # 分析总结
└── plots/                           # 可视化图表（如果启用）
    ├── PCA_PC1_vs_PC2.pdf
    ├── PCA_PC1_vs_PC3.pdf
    └── ...
```

## 输出文件格式 | Output File Formats

### 特征向量文件 | Eigenvector File

```
Sample1   PC1:0.0234   PC2:-0.0156   PC3:0.0089   ...
Sample2   PC1:-0.0189  PC2:0.0223    PC3:-0.0123  ...
Sample3   PC1:0.0087   PC2:-0.0091   PC3:0.0156    ...
```

### 解释方差表 | Explained Variance Table

```
PC      Explained_Variance_Ratio    Cumulative_Variance_Ratio
PC1     0.045678                    0.045678
PC2     0.032345                    0.078023
PC3     0.028912                    0.106935
...
```

## 使用示例 | Usage Examples

### 示例1：基础PCA分析 | Example 1: Basic PCA Analysis

```bash
biopytools vcf2pca \
    -i population.vcf.gz \
    -o pca_basic
```

**适用场景：** 快速查看样本群体结构

### 示例2：带聚类的完整分析 | Example 2: Analysis with Clustering

```bash
biopytools vcf2pca \
    -i population.vcf.gz \
    -o pca_cluster \
    --cluster \
    --cluster-method kmeans \
    --cluster-k 3 \
    --components 20 \
    --threads 16
```

**适用场景：** 需要识别亚群分组

### 示例3：严格质控的PLINK分析 | Example 3: PLINK Analysis with Strict QC

```bash
biopytools vcf2pca \
    -i population.vcf.gz \
    -o pca_plink \
    --backend plink \
    --maf 0.01 \
    --missing 0.05 \
    --hwe 1e-10 \
    --plot \
    --components 15
```

**适用场景：** 需要严格质控和可视化

### 示例4：使用自定义工具路径 | Example 4: Custom Tool Paths

```bash
biopytools vcf2pca \
    -i variants.vcf \
    -o pca_custom \
    --vcf2pca-path /opt/VCF2PCACluster/bin/VCF2PCACluster \
    --threads 8
```

**适用场景：** 工具安装在非默认路径

## 常见问题 | FAQ

### Q1: 如何选择合适的后端？
**A:**
- 数据量大（>50万SNP）：推荐V2P后端
- 需要聚类分析：推荐V2P后端
- 需要严格质控：推荐PLINK后端
- 不确定时：先尝试V2P后端（默认）

### Q2: 聚类方法如何选择？
**A:**
- `kmeans`: 最常用，适合球形簇，需要预先指定K值
- `dbscan`: 自动发现簇数，适合密度不均的数据
- `em`: 基于概率模型，适合重叠簇

### Q3: 如何确定主成分数量？
**A:**
- 查看解释方差表，选择累积解释方差达到70-90%的PC数
- 观察scree plot的拐点
- 通常前10-20个PC已足够

### Q4: V2P和PLINK结果一致吗？
**A:** 是的，V2P后端的PCA结果与PLINK、tassel、gapit、gcta完全一致。

### Q5: 如何解读聚类结果？
**A:** 聚类结果文件包含每个样本的簇标签，可以：
- 结合地理信息分析群体结构
- 识别离群样本
- 验证样本分组

## 引用 | Citation

如果使用VCF2PCACluster后端，请引用：

```
doi: https://doi.org/10.1186/s12859-024-05770-1
```

## 更新日志 | Changelog

### v2.0.0 (2026-03-23)
- 重命名模块：vcf_pca → vcf2pca
- 集成VCF2PCACluster后端
- 支持双后端切换
- 新增聚类分析功能
- 统一命名规范

### v1.0.0
- 初始版本，基于PLINK的PCA分析

## 相关工具 | Related Tools

- `vcf2genotype` - VCF基因型提取
- `vcf2nj` - VCF邻接树构建
- `vcf2phylip` - VCF转phylip格式
- `vcf2gene` - VCF变异基因注释

## 技术支持 | Support

如有问题或建议，请联系：xiang.li@yourlab.edu
