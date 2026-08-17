# JanusX GWAS 与基因组选择 | JanusX GWAS and Genomic Selection

一句话理解：**一个「统一入口」做全基因组关联分析(GWAS)和基因组选择(GS)两类统计，外加 PCA 和结果可视化**——输入基因型(VCF 或 PLINK)+ 表型，就能跑多种模型、出关联结果和图。

## 功能概述 | Overview

- 封装 JanusX(Rust 高性能二进制)的统一接口，4 个子命令：`gwas`(关联分析)、`gs`(基因组选择)、`pca`(主成分分析)、`postgwas`(结果可视化+注释)
- GWAS 支持 4 种模型：`lm`(一般线性)、`lmm`(混合线性)、`fastlmm`、`farmcpu`
- GS 支持 5 种模型：`GBLUP`、`rrBLUP`、`BayesA`、`BayesB`、`BayesCpi`
- 基因型支持 VCF(.vcf/.vcf.gz) 或 PLINK binary(.bed/.bim/.fam 前缀)
- 底层是「拼一条命令 → 调 JanusX 二进制跑一次 → 列出输出」，**Python 层没有多步流程、没有断点续传**

## 快速开始 | Quick Start

```bash
biopytools janusx gwas -i data.vcf.gz -p pheno.txt -m lmm -o gwas_out
```

最小输入：一个 VCF + 一个 TSV 表型文件(第一列样本名，后面是性状值)，用默认 lmm 模型跑关联分析。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GWAS(关联分析) | 全基因组扫描，找出「哪些位点与某个性状(如株高、抗病)显著相关」 |
| GS(基因组选择) | 不找位点，而是用全基因组信息**预测**一个个体未来会表现多好(像「基因算命」) |
| 模型(model) | 统计上怎么算的方法；lm 最快但不考虑亲缘，lmm/farmcpu 校正亲缘关系更稳 |
| 亲缘关系矩阵(GRM) | 个体间基因相似度表格，用于排除「有亲戚关系导致的假关联」 |
| 主成分分析(PCA) | 把高维基因型压缩成几个「主轴」，常用来画样本的分群散点图 |
| 曼哈顿图(Manhattan) | 把每个位点的显著性(-log10(P))按染色体位置画成柱状，像城市天际线，冒尖的塔就是候选位点 |
| QQ图(QQ plot) | 检验「整体信号是否正常」的诊断图：点都贴在对角线上=正常，右上角翘起=有真信号 |
| MAF | 「少数派碱基」在群体里占多少；太低=位点没信息量，默认过滤 <0.05 |
| GEBV | 基因组估计育种值，GS 输出的「预测成绩」 |

## 输入 | Input

### 基因型文件

- VCF：`.vcf` 或 `.vcf.gz`(默认 `-t vcf`)；或 PLINK binary 前缀(`-t bfile`，即 .bed/.bim/.fam 三件套去掉扩展名的前缀)。
- 样本名须与表型文件第一列一致。

### 表型文件

TSV 制表符分隔，第一列样本名，其余列是性状值：

```text
sample1    10.5    0.85
sample2    12.3    0.92
sample3    11.8    0.78
```

- 多个性状列时，用 `-n` 指定分析哪一列(零基，从 0 数)。
- 可选协变量文件 `-c`；可选 Q 矩阵/PCA `-q`(传数量) 或 `--grm` 传预计算 GRM。

## 参数说明 | Parameters

### gwas - 关联分析 | GWAS

**通俗理解|In plain words:** `-m` 选模型——没把握就先默认 `lmm`(稳)，样本巨大想快用 `farmcpu`。`--maf`/`--geno` 是位点过滤(默认 0.05/0，一般不用动)。`-q` 给「群体结构」几个主成分做校正，`-k` 控制亲缘矩阵算法(默认 1，一般不用动)。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --genotype` | 必填 | 基因型(VCF 或 PLINK 前缀) |
| `-p, --pheno` | 必填 | 表型文件 |
| `-t, --type` | `vcf` | 基因型类型 vcf/bfile |
| `-m, --models` | `lmm` | 模型 lm/lmm/fastlmm/farmcpu(可多选) |
| `--maf` | `0.05` | 最小等位基因频率阈值 |
| `--geno` | `0.0` | 缺失率阈值(0=不过滤) |
| `-k, --grm` | `1` | 亲缘矩阵方法 |
| `-q, --qcov` | `0` | PCA 数量或 Q 矩阵文件(0=不用) |
| `-c, --cov` | 无 | 协变量文件 |
| `-n, --ncol` | 无 | 表型列索引(零基) |
| `--plot` | 关 | 生成图表 |
| `-th, --threads` | `12` | 线程数 |
| `-o, --output-dir` | `./janusx_gwas_output` | 输出目录 |
| `--prefix` | 自动 | 输出前缀(默认=基因型文件名) |

### gs - 基因组选择 | Genomic Selection

**通俗理解|In plain words:** 用全基因组预测个体的「成绩」。`-m` 选预测模型(默认 `GBLUP` 最常用)。`--cv` 是交叉验证折数——给了它会评估预测准不准，适合先跑一遍验证模型可信度；`--pcd` 开启 PCA 降维(位点极多时可提速)。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --genotype` | 必填 | 基因型 |
| `-p, --pheno` | 必填 | 表型文件 |
| `-m, --models` | `GBLUP` | 模型 GBLUP/rrBLUP/BayesA/BayesB/BayesCpi(可多选) |
| `--cv` | 无 | 交叉验证折数(评估预测精度) |
| `--pcd` | 关 | 启用 PCA 降维 |
| `-n, --ncol` | 无 | 表型列索引(零基) |
| `-o, --output-dir` | `./janusx_gs_output` | 输出目录 |

### pca - 主成分分析 | PCA

**通俗理解|In plain words:** 画「样本分群图」。`-d` 控制输出几个主成分；`--plot` 出 2D 散点，`--plot3d` 出 3D 旋转 GIF；`-g` 给个分组文件就能把不同群体染成不同颜色。**注意输入优先级：给了 `--grm` 或 `--pcfile` 就用预计算结果，否则才读基因型重新算。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --genotype` | 必填 | 基因型(无 --grm/--pcfile 时) |
| `--grm` / `--pcfile` | 无 | 预计算 GRM/PCA 前缀(优先使用) |
| `-d, --dim` | `3` | 输出主成分数量 |
| `--plot` / `--plot3d` | 关 | 2D 散点 / 3D 旋转 GIF |
| `-g, --group` | 无 | 分组文件(染色用) |
| `--color` | `1` | 调色板索引(0-6) |
| `-o, --output-dir` | `./janusx_pca_output` | 输出目录 |

### postgwas - 结果可视化与注释 | PostGWAS

**通俗理解|In plain words:** 把 GWAS 结果(TSV)画成曼哈顿图/QQ图，并按阈值圈显著位点做基因注释。`--chr-col`/`--pos-col`/`--pvalue-col` 告诉它结果文件里三列分别叫什么；`--threshold` 设显著性线；`-a` + `--ab` 给注释文件和窗口做候选基因标注。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-f, --files` | 必填 | GWAS 结果文件(可多个) |
| `--chr-col` / `--pos-col` / `--pvalue-col` | `#CHROM`/`POS`/`p` | 染色体/位置/P值列名 |
| `--threshold` | 无 | 显著性阈值 |
| `--format` | `png` | 输出格式 pdf/png/svg/tif |
| `-a, --anno` / `--ab` | 无 | 注释文件 / 注释窗口(kb) |
| `-o, --output-dir` | `./janusx_postgwas_output` | 输出目录 |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 四个子命令各自「独立成章」，没有跨子命令的自动流水线。每个子命令就是把参数拼成一条 JanusX 命令交给二进制执行，内部过滤、建矩阵、跑模型都由 JanusX 完成。

```text
gwas:    基因型 + 表型  → JanusX gwas → 关联结果 TSV + 图(可选)
gs:      基因型 + 表型  → JanusX gs   → 预测结果(GEBV) + 图(可选)
pca:     基因型(或GRM) → JanusX pca  → 主成分坐标 + 散点图(可选)
postgwas:GWAS结果TSV   → JanusX postGWAS → 曼哈顿图/QQ图 + 显著位点注释
```

> 典型串联用法：先 `pca` 看群体结构 → 再 `gwas` 时用 `-q 3` 把前 3 个主成分当协变量校正 → 最后 `postgwas` 可视化。

## 输出 | Output

### gwas / gs / pca

结果文件名由 JanusX 二进制决定(源码只做 `glob {prefix}.*` 列清单)，一般按 `{prefix}.{trait}.{model}.tsv` 之类命名，`--plot` 时附带 `.png`/`.svg` 图。日志在输出目录根下：`janusx_gwas.log` / `janusx_gs.log` / `janusx_pca.log`。

### postgwas

```text
janusx_postgwas_output/
├── {输入文件名}.manh.png      # 曼哈顿图(每个输入结果一张)
├── {输入文件名}.qq.png        # QQ图
├── {阈值}.anno.tsv            # 显著位点注释(仅给 --threshold 时)
└── janusx_postgwas.log        # 日志
```

## 结果解读 | Interpreting Results

### GWAS 结果(TSV)

**通俗理解|In plain words:** 一行一个位点，看 P 值列——P 越小越显著。画成曼哈顿图，`-log10(P)` 越高(塔越高)越像真的关联位点。

- 关注 `P值` 或 `-log10(P)` 列；P ≤ 1e-5 或 1e-6 是常见候选线。
- QQ 图里右上角「翘尾巴」说明有超出随机预期的强信号；整体都贴线说明没有显著信号。

### GS 结果(GEBV)

**通俗理解|In plain words:** 每个个体一个「预测成绩」。用了 `--cv` 会额外给出预测精度(相关性)，精度越高(越接近 1)说明模型越可信。

### PCA 结果

**通俗理解|In plain words:** 每个样本在主成分坐标轴上的位置。散点图上靠得近的样本基因越像；用 `-g` 分组染色后能一眼看出群体结构。

### postgwas 注释(anno.tsv)

**通俗理解|In plain words:** 把超过阈值的显著位点「翻译」成附近有哪些基因，方便直接看候选基因。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规 GWAS | `-m lmm`(默认)，样本特别大改用 `-m farmcpu` |
| 有群体结构 | 先 `pca` 看结构，`gwas` 时加 `-q 3` 校正 |
| 多模型对比 | `-m lm lmm farmcpu` 一次跑多个，对比结果一致性 |
| 基因组选择 | 默认 `GBLUP`；先 `--cv 5` 验证精度再上正式数据 |
| 位点极多(>千万) | `gs` 加 `--pcd` 降维提速 |
| 只想画图/注释 | `postgwas -f result.tsv --threshold 1e-5 -a gene.gff` |
| 基因型已是 PLINK | `-t bfile` 并用 .bed/.bim/.fam 前缀 |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genotype` | 必填 |  | 基因型文件(VCF或PLINK前缀)｜Genotype file (VCF or PLINK prefix) |
| `-p, --pheno` | 必填 |  | 表型文件｜Phenotype file |
| `-t, --type` | `vcf` | vcf/bfile | 基因型类型｜Genotype type |
| `-m, --models` | `['lmm']` | lm/lmm/fastlmm/farmcpu | GWAS模型｜GWAS models |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold |
| `--geno` | `0.0` | float | 缺失率阈值(0=不过滤｜Missing rate threshold, 0=no filtering) |
| `-k, --grm` | `1` |  | 亲缘关系矩阵方法｜GRM method |
| `-q, --qcov` | `0` | int | PCA数量或Q矩阵文件路径｜Number of PCs or path to Q matrix file (0=不使用Q矩阵｜no Q matrix) |
| `-c, --cov` | — |  | 协变量文件｜Covariate file |
| `-n, --ncol` | — | int | 表型列索引(零基)｜Phenotype column indices (zero-based) |
| `--plot` | — |  | 生成图表｜Generate plots |
| `--chunksize` | `100000` | int | SNP分块大小｜SNP chunk size |
| `--mmap-limit` | — | int | 内存映射限制｜Memory map limit |
| `-th, --threads` | `12` | int | 线程数｜Number of threads |
| `-o, --output-dir` | `./janusx_gwas_output` | Path | 输出目录｜Output directory |
| `--prefix` | — |  | 输出文件前缀｜Output file prefix |
| `--janusx-path` | — | Path | JanusX可执行文件路径｜JanusX executable path |
| `--pcd` | — |  | 启用PCA降维｜Enable PCA-based dimensionality reduction |
| `--cv` | — | int | 交叉验证折数｜K-fold cross-validation |
| `--grm` | — |  | 预计算的GRM前缀｜Precomputed GRM prefix |
| `--pcfile` | — |  | 预计算的PCA文件｜Precomputed PCA file |
| `-d, --dim` | `3` | int | 输出主成分数量｜Number of PCs to output |
| `--plot3d` | — |  | 生成3D旋转GIF｜Generate 3D rotating GIF |
| `-g, --group` | — |  | 分组文件路径｜Group file path |
| `--color` | `1` | int | 调色板索引(0-6)｜Color palette index |
| `-f, --files` | 必填 |  | GWAS结果文件列表｜List of GWAS result files |
| `--chr-col` | `#CHROM` |  | 染色体列名｜Chromosome column name |
| `--pos-col` | `POS` |  | 位置列名｜Position column name |
| `--pvalue-col` | `p` |  | P值列名｜P-value column name |
| `--threshold` | — | float | 显著性阈值｜Significance threshold |
| `--noplot` | — |  | 禁用绘图｜Disable plotting |
| `--hl, --highlight` | — |  | 高亮区域BED文件｜Highlight regions BED file |
| `--format` | `png` | pdf/png/svg/tif | 输出格式｜Output format |
| `-a, --anno` | — |  | 注释文件路径｜Annotation file path |
| `--ab, --anno-broaden` | — | int | 注释窗口｜Annotation window (kb) |
| `--desc-item` | `description` |  | GFF描述键｜GFF description key |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genotype` | 必填 |  | 基因型文件(VCF或PLINK前缀)｜Genotype file (VCF or PLINK prefix) |
| `-p, --pheno` | 必填 |  | 表型文件｜Phenotype file |
| `-t, --type` | `vcf` | vcf/bfile | 基因型类型｜Genotype type (default: vcf) |
| `-m, --models` | `['lmm']` | lm/lmm/fastlmm/farmcpu | GWAS模型｜GWAS models (default: lmm) |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold (default: 0.05) |
| `--geno` | `0.0` | float | 缺失率阈值(0=不过滤｜Missing rate threshold, 0=no filtering) |
| `-k, --grm` | `1` |  | 亲缘关系矩阵方法｜GRM method (default: 1) |
| `-q, --qcov` | `0` | int | PCA数量或Q矩阵文件路径｜Number of PCs or path to Q matrix file (default: 0, no Q matrix) |
| `-c, --cov` | — |  | 协变量文件｜Covariate file |
| `-n, --ncol` | — | int | 表型列索引(零基)｜Phenotype column indices (zero-based) |
| `--plot` | — | store_true | 生成图表｜Generate plots |
| `--chunksize` | `100000` | int | SNP分块大小｜SNP chunk size (default: 100000) |
| `--mmap-limit` | — | int | 内存映射限制｜Memory map limit |
| `-th, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-o, --output-dir` | `./janusx_gwas_output` |  | 输出目录｜Output directory (default: ./janusx_gwas_output) |
| `--prefix` | — |  | 输出文件前缀｜Output file prefix |
| `--janusx-path` | — |  | JanusX可执行文件路径｜JanusX executable path |
| `--pcd` | — | store_true | 启用PCA降维｜Enable PCA-based dimensionality reduction |
| `--cv` | — | int | 交叉验证折数｜K-fold cross-validation |
| `--grm` | — |  | 预计算的GRM前缀｜Precomputed GRM prefix |
| `--pcfile` | — |  | 预计算的PCA文件｜Precomputed PCA file |
| `-d, --dim` | `3` | int | 输出主成分数量｜Number of PCs to output (default: 3) |
| `--plot3d` | — | store_true | 生成3D旋转GIF｜Generate 3D rotating GIF |
| `-g, --group` | — |  | 分组文件路径｜Group file path |
| `--color` | `1` | int | 调色板索引(0-6)｜Color palette index (default: 1) |
| `-f, --files` | 必填 |  | GWAS结果文件列表｜List of GWAS result files |
| `--chr-col` | `#CHROM` |  | 染色体列名｜Chromosome column name (default: #CHROM) |
| `--pos-col` | `POS` |  | 位置列名｜Position column name (default: POS) |
| `--pvalue-col` | `p` |  | P值列名｜P-value column name (default: p) |
| `--threshold` | — | float | 显著性阈值｜Significance threshold |
| `--noplot` | — | store_true | 禁用绘图｜Disable plotting |
| `--hl, --highlight` | — |  | 高亮区域BED文件｜Highlight regions BED file |
| `--format` | `png` | pdf/png/svg/tif | 输出格式｜Output format (default: png) |
| `-a, --anno` | — |  | 注释文件路径｜Annotation file path |
| `--ab, --anno-broaden` | — | int | 注释窗口｜Annotation window (kb) |
| `--desc-item` | `description` |  | GFF描述键｜GFF description key (default: description) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- JanusX 二进制(非 conda)：默认 `~/software/JanusX/JanusX.bin`，可用环境变量 `JANUSX_PATH` 或 `--janusx-path` 覆盖
- Python 3，无第三方库依赖(纯 subprocess 封装)

## 常见问题 | FAQ

**Q1：`--maf`/`--geno` 默认值到底是多少？**
源码里是 `--maf 0.05`、`--geno 0.0`(即不过滤缺失)。旧文档写的 0.02/0.05 已过时。

**Q2：输出文件叫什么名字？**
gwas/gs/pca 的结果文件名由 JanusX 二进制产出，源码只是列出 `{prefix}.*` 文件；跑完看输出目录和日志里的「Output files」清单即可。postgwas 的图名是确定的(`.manh.<格式>`、`.qq.<格式>`)。

**Q3：能断点续传吗？**
不能。每个子命令都是「跑一次 JanusX」，重跑会整体重来。

**Q4：`--grm` 和 `--pcfile` 什么时候用？**
pca 子命令里，给了 `--grm`(或 `--pcfile`) 就跳过重新算，直接用预计算结果画图，省时间。

**Q5：表型文件多列怎么选？**
用 `-n` 指定列索引(零基)，例如 `-n 1` 表示第 2 列；GWAS 支持多个 `-n` 一次分析多列。
