# K-mer GWAS分析模块

**K-mer based Genome-Wide Association Study Analysis Module**

## 功能概述|Overview

K-mer GWAS分析模块是一个完整的基于k-mer的全基因组关联分析工具包，基于KMERIA软件构建，提供从k-mer计数到关联分析的完整流程。支持完整流程一键式运行和各步骤独立运行，适用于二倍体和多倍体物种的关联分析。

## 主要特性|Key Features

- **完整分析流程**: k-mer计数 → 矩阵构建 → 过滤 → 格式转换 → 关联分析
- **一键式运行**: 支持完整流程自动化运行，适合大规模样本分析
- **断点续传**: 支持从任意步骤继续运行
- **质控统计**: 自动生成各步骤QC报告
- **可视化**: 自动生成曼哈顿图和QQ图
- **k-mer注释**: 支持关联k-mer的基因组定位和功能注释
- **灵活参数**: 支持自定义k-mer大小、丰度阈值、倍性等参数
- **高效处理**: 优化的批处理和并行计算

## 安装|Installation

### 1. KMERIA软件安装

```bash
# KMERIA已安装在
/share/org/YZWL/yzwl_lixg/software/kmeria

# Conda环境
/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmeriaenv
```

### 2. 环境变量配置

```bash
# 添加到 ~/.zshrc 或 ~/.bashrc
export KMERIA_HOME=/share/org/YZWL/yzwl_lixg/software/kmeria
export KMERIA_CONDA=/share/org/YZWL/yzwl_lixg/miniforge3/envs/kmeriaenv
export LD_LIBRARY_PATH=${KMERIA_HOME}/lib:${KMERIA_CONDA}/lib:$LD_LIBRARY_PATH
export PATH=${KMERIA_HOME}/bin:${KMERIA_HOME}/bimbamAsso:${KMERIA_HOME}/external_tools:$PATH
```

## 快速开始|Quick Start

### 准备数据 - 输入文件格式详解

#### 1. FASTQ文件目录 (`--fastq-dir`)

**目录结构示例**：
```
data/fastq/
├── sample1_R1.fq.gz
├── sample1_R2.fq.gz
├── sample2_R1.fq.gz
├── sample2_R2.fq.gz
├── sample3_1.fq.gz
├── sample3_2.fq.gz
└── ...
```

**文件命名规则**：
- 支持双端测序命名：`sample_R1.fq.gz` + `sample_R2.fq.gz`
- 或：`sample_1.fq.gz` + `sample_2.fq.gz`
- 自动识别配对关系
- 支持格式：`.fq.gz`, `.fastq.gz`, `.fq`, `.fastq`

---

#### 2. 样本列表文件 (`--samples`)

**格式**: 纯文本，每行一个样本名

**示例** (`samples.txt`):
```text
sample1
sample2
sample3
SRR28578485
SRR28578484
...
```

**说明**：
- 样本名必须与FASTQ文件的前缀一致
- 例如：`sample1_R1.fq.gz` 的样本名为 `sample1`
- 不要包含文件扩展名
- 每行一个样本，无表头
- 空行和以 `#` 开头的行会被忽略

---

#### 3. 测序深度文件 (`--depth-file`)

**格式**: Tab分隔 (TSV)

**示例** (`depth.txt`):
```text
sample1	45.2
sample2	52.8
sample3	38.9
SRR28578485	143.2531311
SRR28578484	151.7451801
```

**说明**：
- 第1列：样本名（与样本列表和FASTQ文件一致）
- 第2列：测序深度（浮点数，单位：X）
- 无表头行
- Tab分隔 (`\t`)
- 深度值通常从比对结果统计得到（如 `samtools depth`）

**混合倍性物种格式**（可选）：
```text
sample1	45.2	4
sample2	52.8	8
sample3	38.9	4
```
- 第3列：倍性（可选，如果不提供则使用 `--ploidy` 参数的值）

---

#### 4. 表型文件 (`--pheno-file`)

**格式**: Tab或空格分隔

**示例 1 - 简单格式** (`pheno.txt`):
```text
sample1	1.5
sample2	2.3
sample3	1.8
SRR28578485	0.015350
SRR28578484	-0.152834
SRR28578303	0.375435
```

**示例 2 - 带协变量的格式**:
```text
sample1	1.5	0	0.5	male
sample2	2.3	1	-0.2	female
sample3	1.8	0	0.8	male
```

**说明**：
- 第1列：样本名（必须与样本列表一致）
- 第2列：表型值（数值型，连续或分类变量）
  - 连续型性状：如株高、产量等（小数）
  - 分类性状：0/1 或其他编码
- 第3列及以后：可选的协变量（如群体结构、性别等）
- 默认使用第2列作为表型，可通过 `--pheno-col` 参数指定
- 支持 Tab 或空格分隔
- 无表头行（或有表头行时会被第一行数据覆盖）

**表型值示例**：
```text
# 连续型性状（如产量）
sample001	25.6
sample002	30.2
sample003	28.7

# 二分类性状（如抗病/感病）
sample001	0
sample001	1
sample002	1

# 标准化后的表型值（mean=0, sd=1）
sample001	0.015350
sample002	-0.152834
sample003	0.375435
```

---

### 输入文件格式验证

**样本名一致性检查**：
```bash
# 检查样本列表中的所有样本名是否在深度文件中
awk 'NR==FNR{a[$1]=1; next} !($1 in a){print "Missing in depth:", $1}' samples.txt depth.txt

# 检查深度文件中的所有样本名是否在表型文件中
awk 'NR==FNR{a[$1]=1; next} !($1 in a){print "Missing in pheno:", $1}' depth.txt pheno.txt
```

**统计信息**：
```bash
# 统计样本数量
wc -l samples.txt

# 查看深度分布
awk '{print $2}' depth.txt | sort -n | head -20
awk '{print $2}' depth.txt | sort -n | tail -20

# 查看表型值分布
awk '{print $2}' pheno.txt | sort -n | head -20
awk '{print $2}' pheno.txt | sort -n | tail -20
```

### 完整流程运行

```bash
# 基本用法
biopytools kmeria pipeline \
    -i /data/fastq \
    --samples samples.txt \
    -d depth.txt \
    -p pheno.txt \
    -o /data/kmeria_results \
    -t 24

# 自定义参数
biopytools kmeria pipeline \
    -i /data/fastq \
    --samples samples.txt \
    -d depth.txt \
    -p pheno.txt \
    -o results \
    -k 31 \
    --min-abund 5 \
    --max-abund 1000 \
    --missing-ratio 0.6 \
    --ploidy 2 \
    -t 32 \
    --enable-qc \
    --enable-visualization
```

### 分步运行

```bash
# Step 1: k-mer计数
biopytools kmeria count \
    -i /data/fastq \
    --samples samples.txt \
    -o 01_kmer_counts \
    -k 31 \
    -t 24

# Step 2: 矩阵构建
biopytools kmeria kctm \
    -i 01_kmer_counts \
    -o 02_kmer_matrices \
    -t 24

# Step 3: 过滤
biopytools kmeria filter \
    -i 02_kmer_matrices \
    -o 03_filtered_matrices \
    -d depth.txt \
    -p 2 \
    -t 24

# Step 4: 转换为BIMBAM格式
biopytools kmeria m2b \
    -i 03_filtered_matrices \
    -o 04_bimbam \
    -t 24

# Step 5: 关联分析
biopytools kmeria asso \
    -i 04_bimbam \
    -p pheno.txt \
    -o 05_association \
    -t 64
```

## 参数说明|Parameters

### pipeline - 完整流程

#### 必需参数|Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --fastq-dir` | FASTQ文件目录| `-i /data/fastq` |
| `--samples` | 样本列表文件| `--samples samples.txt` |
| `-d, --depth-file` | 测序深度文件| `-d depth.txt` |
| `-p, --pheno-file` | 表型文件| `-p pheno.txt` |

#### k-mer参数|K-mer Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-k, --kmer-size` | 31 | K-mer大小 (2-31) |
| `--min-abund` | 5 | 最小k-mer丰度 |
| `--max-abund` | 1000 | 最大k-mer丰度 |

#### 过滤参数|Filter Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--missing-ratio` | 0.6 | 缺失率阈值 (0-1) |
| `--ploidy` | 2 | 基因组倍性 |

#### 性能参数|Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | 24 | 线程数 |
| `--batch-size` | 4 | 批处理大小 |

#### 流程控制|Pipeline Control

| 参数 | 描述 |
|------|------|
| `--step` | 从指定步骤开始 (count/kctm/filter/m2b/asso) |

#### 可选功能|Optional Features

| 参数 | 描述 |
|------|------|
| `--enable-qc` | 启用质控统计 (默认启用) |
| `--enable-visualization` | 启用可视化 (默认启用) |
| `--enable-annotation` | 启用k-mer注释 |
| `--genome-file` | 参考基因组 (注释用) |
| `--gff-file` | GFF注释文件 (注释用) |

### count - k-mer计数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --fastq-dir` | *必需* | FASTQ文件目录 |
| `--samples` | *必需* | 样本列表文件 |
| `-o, --output-dir` | `./01_kmer_counts` | 输出目录 |
| `-k, --kmer-size` | 31 | K-mer大小 |
| `-t, --threads` | 24 | 线程数 |
| `-b, --batch-size` | 4 | 批处理大小 |
| `-C, --count-separate-strands` | False | 分别计数链 |
| `-T, --text-output` | False | 文本输出 |

### kctm - 矩阵构建

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | *必需* | 输入目录 |
| `-o, --output-dir` | `./02_kmer_matrices` | 输出目录 |
| `-t, --threads` | 24 | 线程数 |

### filter - 过滤

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | *必需* | 输入目录 |
| `-o, --output-dir` | `./03_filtered_matrices` | 输出目录 |
| `-d, --depth-file` | *必需* | 测序深度文件 |
| `-c, --max-abund` | 1000 | 最大丰度 |
| `-s, --missing-ratio` | 0.6 | 缺失率 |
| `-p, --ploidy` | 2 | 倍性 |
| `-t, --threads` | 24 | 线程数 |

### m2b - 格式转换

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --in` | *必需* | 输入目录 |
| `-o, --out` | `./04_bimbam` | 输出目录 |
| `-t, --threads` | 24 | 线程数 |
| `--no-normalize` | False | 不归一化 |
| `--quantile-norm` | False | 分位数归一化 |

### asso - 关联分析

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | *必需* | 输入目录 |
| `-p, --pheno-file` | *必需* | 表型文件 |
| `-o, --output-dir` | `./05_association` | 输出目录 |
| `-n, --pheno-col` | 1 | 表型列 |
| `-c, --covar-file` | - | 协变量文件 |
| `-k, --kinship-file` | - | 亲缘关系矩阵 |
| `-t, --threads` | 64 | 线程数 |

## 输出结果|Output

### 目录结构

```
kmeria_results/
├── 01_kmer_counts/          # k-mer计数结果
│   ├── sample1_k31.bin
│   └── sample2_k31.bin
├── 02_kmer_matrices/        # k-mer矩阵
│   ├── kmer_matrix.0001.bin
│   └── kmer_matrix.0002.bin
├── 03_filtered_matrices/    # 过滤后的矩阵
│   ├── filtered_0001.txt
│   └── filtered_0002.txt
├── 04_bimbam/              # BIMBAM格式文件
│   ├── *.bimbam.gz
│   └── sample_list.txt
├── 05_association/         # 关联分析结果
│   ├── association_results.txt
│   └── significant_kmers.txt
├── qc_reports/             # QC报告
│   └── qc_report.json
├── visualization/          # 可视化结果
│   ├── manhattan_plot.pdf
│   └── qq_plot.pdf
├── annotation/             # k-mer注释
│   └── kmer_annotations.json
└── commands_log.json       # 命令日志
```

## 使用示例|Examples

### 示例1: 1000样本完整分析

```bash
biopytools kmeria pipeline \
    -i /data/1000_samples/fastq \
    --samples /data/1000_samples/samples.txt \
    -d /data/1000_samples/depth.txt \
    -p /data/1000_samples/trait_phenotype.txt \
    -o /data/1000_samples/kmeria_gwas \
    -k 31 \
    -t 32 \
    --batch-size 10 \
    --ploidy 2 \
    --enable-qc \
    --enable-visualization
```

### 示例2: 从过滤步骤继续

```bash
biopytools kmeria pipeline \
    -i /data/fastq \
    --samples samples.txt \
    -d depth.txt \
    -p pheno.txt \
    -o results \
    --step filter \
    -t 32
```

### 示例3: 只运行k-mer计数

```bash
biopytools kmeria count \
    -i /data/fastq \
    --samples samples.txt \
    -o 01_kmer_counts \
    -k 31 \
    -t 32
```

### 示例4: 自定义参数的关联分析

```bash
biopytools kmeria asso \
    -i 04_bimbam \
    -p pheno.txt \
    -o 05_association \
    -c covariates.txt \
    -k kinship_matrix.txt \
    -t 64
```

## 注意事项|Important Notes

1. **测序深度**: 建议平均深度 ≥ 30X
2. **样本数量**: 建议至少100个样本以获得可靠的关联结果
3. **k-mer大小**:
   - 较小k-mer (15-21): 更敏感，计算更快，但假阳性更高
   - 较大k-mer (25-31): 更特异，假阳性更少，但计算资源需求更大
4. **倍性设置**:
   - 二倍体: `--ploidy 2`
   - 四倍体: `--ploidy 4`
   - 六倍体: `--ploidy 6`
5. **内存需求**:
   - 1000样本分析建议至少500GB内存
   - 使用`--batch-size`参数控制内存使用

## 故障排除|Troubleshooting

### 问题1: "kmeria: command not found"

**解决方案**:
```bash
# 检查环境变量
echo $PATH | grep kmeria
echo $LD_LIBRARY_PATH | grep kmeria

# 重新加载环境变量
source ~/.zshrc  # 或 source ~/.bashrc
```

### 问题2: 内存不足

**解决方案**:
```bash
# 减小批处理大小
--batch-size 2

# 减小k-mer大小
-k 21

# 增加过滤阈值
--missing-ratio 0.8
```

### 问题3: 某步骤失败

**解决方案**:
```bash
# 从失败的步骤继续
--step filter  # 从filter步骤继续
```

## 参考文献|References

```bibtex
@article{chen2025kmeria,
  title={A k-mer-based GWAS approach empowering gene mining in polyploids},
  author={Chen, Shuai and others},
  journal={Research Square},
  year={2025},
  doi={10.21203/rs.3.rs-7347406/v1}
}
```

## 相关链接|Related Links

- **KMERIA GitHub**: https://github.com/Sh1ne111/KMERIA
- **KMERIA Wiki**: https://github.com/Sh1ne111/KMERIA/wiki
- **biopytools文档**: https://github.com/your-org/biopytools

## 许可证|License

MIT License

---

**版本**: 1.0.0
**更新日期**: 2026-01-13
**作者**: BioPyTools Team
