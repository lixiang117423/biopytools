# VCF2PCA 主成分分析 | VCF Principal Component Analysis

一句话理解：**从 VCF 基因型数据算出主成分(PCA)，把样本在遗传上的亲疏远近压到一张二维图**，用于看群体结构、样本分层、揪出离群样本。

## 功能概述 | Overview

- 双后端：`plink`(默认，支持 SNP/INDEL 等所有变异类型)与 `v2p`(VCF2PCACluster，仅 SNP，带 kinship/聚类)
- 默认**不过滤** VCF，并自动剔除零基因型样本，避免 PLINK 崩溃
- plink 后端：VCF→PLINK→(可选质控)→PCA→结果/可视化
- v2p 后端：一键 PCA + Kinship 矩阵 + 聚类(kmeans/dbscan/em)
- 输出特征值、特征向量、解释方差、合并 Excel 与可视化图

## 快速开始 | Quick Start

```bash
biopytools vcf2pca -i variants.vcf -o pca_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| PCA(主成分分析) | 把成千上万个位点的信息压缩成几个「主方向」，样本的远近关系一目了然 |
| 主成分(PC) | 压缩后的「坐标轴」；PC1 通常解释最多差异，PC2 次之 |
| 特征值 | 每个主成分「解释了多少差异」的大小；越大越重要 |
| 特征向量 | 每个样本在各主成分上的「坐标」，画散点图就靠它 |
| 解释方差比 | 每个主成分解释的差异占比；前几个加起来越高，压缩越「保真」 |
| 群体结构 | 样本按遗传背景自然分成几堆，常对应地理/谱系 |
| Kinship(亲缘矩阵) | 两两样本之间的亲缘关系打分，找重复样本/近亲 |
| 聚类 | 自动把样本分成若干组，不依赖人工标签 |

## 输入 | Input

### VCF 文件

标准 VCF(支持 .vcf/.vcf.gz)，含基因型。

### 样本信息文件(可选，`-s`)

带表头的 TSV，第一列样本名需与 VCF 一致；配合 `-g 列名` 按分组给散点图着色：

```text
sample	population
sample1	popA
sample2	popA
sample3	popB
```

## 分析流程 | Pipeline

plink 后端(默认)：

```text
VCF
  │ 步骤1: VCF → PLINK(自动剔除零基因型样本)
  ▼
步骤2: 质控(默认跳过；--apply-qc 才做 MAF/缺失/HWE 过滤)
  │
  ▼
步骤3: PLINK --pca 计算主成分
  │
  ▼
步骤4: 结果处理(特征值/特征向量/解释方差/Excel)
  │
  └─ 步骤5: 可视化(--plot 时生成散点图/碎石图)
```

v2p 后端：VCF2PCACluster 一步完成 PCA + Kinship + 聚类。

## 输出 | Output

plink 后端(默认)：

```text
pca_output/
├── pca_results.eigenval           # 特征值(PLINK 原始)
├── pca_results.eigenvec           # 特征向量(PLINK 原始)
├── pca_eigenvectors_formatted.txt # 格式化的特征向量(FID/IID/PC1..PCn)
├── pca_explained_variance.txt     # 解释方差表(PC/特征值/方差比/累计)
├── eigenvalues_explanation.txt    # 特征值说明文档
├── pca_results.xlsx               # 合并 Excel(特征向量+解释方差两个sheet)
├── pca_with_sample_info.txt       # PCA 结果与样本信息合并(配合 -s)
├── pca_analysis_summary.txt       # 分析总结报告
├── pca_scree_plot.png             # 碎石图(--plot 时)
├── pca_scatter_plots.png          # PCA 散点图(--plot 时)
├── pca_pairs_plot.png             # 成对散点图(--plot 时)
├── vcf2pca.log                    # 运行日志
└── <样本名前缀>.bed/bim/fam       # 中间 PLINK 文件
```

v2p 后端：

```text
pca_output/
├── vcf2pca.eigenval                # 特征值
├── vcf2pca.eigenvec                # 特征向量
├── pca_sample_point.txt            # 重命名后的特征向量
├── vcf2pca.Normalized_IBS.Kinship  # Kinship 矩阵
├── vcf2pca.cluster                 # 聚类结果(启用 --cluster 时)
├── pca_explained_variance.txt      # 解释方差表
└── vcf2pca.log                     # 运行日志
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 把 `pca_eigenvectors_formatted.txt` 的 PC1、PC2 画成散点图，样本聚成几堆=存在群体结构；某样本孤零零=可能是离群样本或污染。

- **PC1/PC2 散点图**：样本距离越近遗传越相似；按分组着色后，能直观看到群体是否分开
- **pca_explained_variance.txt**：`Explained_Variance_Ratio` 越大该 PC 越重要；前 2-3 个 PC 累计解释 70% 以上，二维图就比较可信
- **Kinship 矩阵(v2p)**：两样本值接近 1 提示重复样本或近亲，值接近 0 提示无关
- **聚类结果(v2p)**：每个样本一个簇标签，可与已知分组交叉验证

## 参数选择建议 | Parameter Guidance

- **后端选择**：默认 `plink`；INDEL 数据必须用 plink(v2p 会跳过 Indel)；要 kinship/聚类用 `-b v2p`(仅 SNP)
- **质控默认关闭**：加 `--apply-qc` 才做 MAF/缺失/HWE 过滤，阈值 `--maf 0.05`/`--missing 0.1`/`--hwe 1e-6` 一般不用动
- **`-c/--components`**：默认 10 个主成分，一般够用；要画 3D 或后续分析可加大
- **`-P/--plot`**：生成散点图/碎石图；配 `-s` + `-g` 按分组着色更直观
- **聚类(v2p)**：`--cluster` 开启，`--cluster-method kmeans` 最常用，`--cluster-k` 设为预期群体数


<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入VCF文件｜Input VCF file |
| `--output, -o` | `./pca_output` | Path | 输出目录｜Output directory |
| `--backend, -b` | `plink` | v2p/plink | 分析后端｜Analysis backend: plink (default, supports SNP/INDEL) or v2p (VCF2PCACluster, SNP-only) |
| `--sample-info, -s` | — | Path | 样本信息文件｜Sample information file |
| `--components, -c` | `10` | int | 主成分数量｜Number of principal components |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold (PLINK backend only) |
| `--missing` | `0.1` | float | 最大缺失率阈值｜Maximum missing rate threshold (PLINK backend only) |
| `--hwe` | `1e-06` | float | Hardy-Weinberg平衡p值阈值｜Hardy-Weinberg equilibrium p-value (PLINK backend only) |
| `--apply-qc` | — |  | 启用质控过滤(MAF/缺失率/HWE,默认不过滤)｜Enable QC filtering (default: no filtering) (PLINK backend only) |
| `--cluster` | — |  | 启用聚类分析｜Enable clustering analysis (V2P backend only) |
| `--cluster-method` | `kmeans` | kmeans/dbscan/em | 聚类方法｜Clustering method: kmeans, dbscan, em (V2P backend only) |
| `--cluster-k` | `3` | int | K-means聚类数｜Number of clusters for K-means (V2P backend only) |
| `--plot, -P` | — |  | 生成PCA可视化图表｜Generate PCA visualization plots |
| `--group-column, -g` | — | str | 分组列名(配合-s样本信息文件,按分组着色)｜Column name for grouping (with sample info file) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--vcf2pca-path` | — | str | VCF2PCACluster路径｜VCF2PCACluster path |
| `--plink-path` | `plink` | str | PLINK软件路径｜PLINK software path |
| `--bcftools-path` | `bcftools` | str | BCFtools软件路径｜BCFtools software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入VCF文件路径｜Input VCF file path |
| `-v, --vcf` | — |  |  |
| `-o, --output` | `./pca_output` |  | 输出目录｜Output directory |
| `-s, --sample-info` | — |  | 样本信息文件｜Sample information file |
| `-b, --backend` | `plink` | v2p/plink | 分析后端｜Analysis backend: plink (default, supports SNP/INDEL) or v2p (VCF2PCACluster, SNP-only) |
| `-c, --components` | `10` | int | 主成分数量｜Number of principal components |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold (PLINK backend only) |
| `--missing` | `0.1` | float | 最大缺失率阈值｜Maximum missing rate threshold (PLINK backend only) |
| `--hwe` | `1e-06` | float | Hardy-Weinberg平衡p值阈值｜Hardy-Weinberg equilibrium p-value (PLINK backend only) |
| `--apply-qc` | — | store_true | 启用质控过滤(MAF/缺失率/HWE,默认不过滤)｜Enable QC filtering (default: no filtering) (PLINK backend only) |
| `--cluster` | — | store_true | 启用聚类分析｜Enable clustering analysis (V2P backend only) |
| `--cluster-method` | `kmeans` | kmeans/dbscan/em | 聚类方法｜Clustering method: kmeans, dbscan, em (V2P backend only) |
| `--cluster-k` | `3` | int | K-means聚类数｜Number of clusters for K-means (V2P backend only) |
| `-P, --plot` | — | store_true | 生成PCA可视化图表｜Generate PCA visualization plots |
| `-g, --group-column` | — |  | 分组列名(配合-s样本信息文件,按分组着色)｜Column name for grouping (with sample info file) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--vcf2pca-path` | — |  | VCF2PCACluster路径｜VCF2PCACluster path |
| `--plink-path` | `plink` |  | PLINK软件路径｜PLINK software path |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- PLINK(默认 `plink`，可 `--plink-path` 覆盖)
- bcftools(默认 `bcftools`)
- VCF2PCACluster(仅 v2p 后端，默认 `~/software/VCF2PCACluster-1.42/bin/VCF2PCACluster`)
- matplotlib/seaborn(生成可视化图时)

## 常见问题 | FAQ

**Q1：默认会过滤 VCF 吗？**
不会。默认**不做任何过滤**，只自动剔除 100% 缺失(零基因型)的样本以避免 PLINK 崩溃。要 MAF/缺失/HWE 过滤需加 `--apply-qc`。

**Q2：INDEL 数据能用吗？**
能用，但必须用默认的 plink 后端。v2p 后端只支持 SNP，会跳过 Indel。

**Q3：plink 和 v2p 结果一样吗？**
PCA 主成分结果一致(v2p 与 PLINK、GCTA 等主流工具一致)；差异在 v2p 额外提供 Kinship 和聚类，plink 后端额外提供质控与 Python 可视化。

**Q4：换参数重跑要删旧文件吗？**
本模块无断点续传，每次运行重算并覆盖同名输出。换 `-c`、`--apply-qc` 等重跑同一 `-o` 即可；想保留多组结果请换输出目录。