# VCF2PCA 主成分分析工具

**VCF文件主成分分析工具，支持VCF2PCACluster和PLINK双后端 | VCF Principal Component Analysis Tool with VCF2PCACluster and PLINK Backends**

## 功能概述 | Overview

VCF2PCA是一个功能强大的VCF文件主成分分析工具，支持两种分析后端：**PLINK**（默认，支持SNP/INDEL等所有变异类型）和**VCF2PCACluster**（仅SNP，带kinship/聚类）。它能够从VCF格式的基因型数据中提取主成分，帮助研究者识别样本群体结构、亲缘关系和种群分层，广泛应用于群体遗传学、进化生物学和生物医学研究。

> ⚠️ **默认行为**：默认使用 PLINK 后端且**不过滤 VCF**（不跑 MAF/缺失率/HWE）。需要质控时加 `--apply-qc`；需要 VCF2PCACluster 的 kinship/聚类时加 `-b v2p`（注意 v2p 仅支持 SNP，会跳过 Indel）。

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

### 默认用法（PLINK后端，不过滤，支持SNP/INDEL）| Default (PLINK backend, no filtering, SNP/INDEL)

```bash
# 基本PCA分析(默认不过滤,自动剔除零基因型样本避免崩溃)
biopytools vcf2pca -i variants.vcf -o pca_results

# INDEL数据同样适用
biopytools vcf2pca -i indels.vcf.gz -o pca_indel

# 分组着色绘图(配合样本信息文件)
biopytools vcf2pca -i variants.vcf -o pca_results -s sample_info.txt -g population -P
```

### 启用质控过滤 | Enable QC Filtering

```bash
# 加 --apply-qc 才跑 MAF/缺失率/HWE 过滤(默认不过滤)
biopytools vcf2pca -i variants.vcf -o pca_results --apply-qc --maf 0.05 --missing 0.1
```

### 使用VCF2PCACluster后端（仅SNP，带kinship/聚类）| Using VCF2PCACluster Backend (SNP-only)

```bash
# 切换v2p后端(注意:仅支持SNP,会跳过Indel)
biopytools vcf2pca -i snp_variants.vcf -o pca_results -b v2p

# 启用聚类分析
biopytools vcf2pca -i snp_variants.vcf -o pca_results -b v2p --cluster --cluster-method kmeans --cluster-k 3
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
| `-b, --backend` | `plink` | 分析后端：`plink` (默认,支持SNP/INDEL) 或 `v2p` (VCF2PCACluster,仅SNP) |

### PCA参数 | PCA Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-c, --components` | `10` | 主成分数量 |
| `-t, --threads` | `12` | 线程数 |

### 质控参数（仅PLINK后端）| Quality Control (PLINK Backend Only)

> 默认不过滤；加 `--apply-qc` 启用过滤后才使用以下阈值。| No filtering by default; thresholds below apply only when `--apply-qc` is set.

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--apply-qc` | `false` | 启用质控过滤（默认不过滤）|Enable QC filtering (default: no filtering) |
| `--maf` | `0.05` | 最小等位基因频率阈值 |
| `--missing` | `0.1` | 最大缺失率阈值 |
| `--hwe` | `1e-6` | Hardy-Weinberg平衡p值阈值 |

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
| `-g, --group-column` | `None` | 分组列名（配合`-s`样本信息文件，按分组着色）|

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
    --apply-qc \
    --maf 0.01 \
    --missing 0.05 \
    --hwe 1e-10 \
    --plot \
    --components 15
```

**适用场景：** 需要严格质控和可视化（默认后端即 PLINK，加 `--apply-qc` 启用过滤）

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
- **默认用 PLINK 后端**：支持 SNP/INDEL 等所有变异类型，默认不过滤
- INDEL 数据：必须用 PLINK 后端（v2p 会跳过 Indel）
- 需要 kinship 矩阵/聚类分析：用 V2P 后端（`-b v2p`，仅 SNP）
- 数据量极大（>50万 SNP）且只要 SNP：V2P 后端更快

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

### v3.0.0 (2026-08-05)
- **合并 vcf_pca 模块**：独立 `biopytools vcf-pca` 已删除，功能并入 `vcf2pca`（核心逻辑本就相同）
- **默认后端改为 PLINK**（原 v2p）：支持 SNP/INDEL 等所有变异类型；v2p 降为可选（`-b v2p`）
- **默认不过滤 VCF**：`--skip-qc` 改为 `--apply-qc`（过滤从默认变为 opt-in）
- **自动剔除零基因型样本**：默认不过滤时，自动剔除 100% 缺失样本以避免 PLINK 崩溃
- 新增 `-g/--group-column` 分组着色参数（原 vcf_pca 功能）

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
