# 📝 HiTE 转座子检测与注释模块

**专业的基因组转座子检测与注释工具 | Professional Genome Transposon Detection and Annotation Tool**

基于 HiTE v3.3.3，通过 Singularity 容器化部署，提供单基因组和泛基因组转座子分析功能

---

## 📖 功能概述 | Overview

HiTE 转座子检测与注释模块是基于 **HiTE v3.3.3** 构建的强大工具，通过 Singularity 容器化部署，提供两种主要分析模式：

1. **HiTE 单基因组分析** (`biopytools hite`): 针对单个基因组的转座子检测与注释
2. **panHiTE 泛基因组分析** (`biopytools panhite`): 针对群体基因组的转座子检测与比较分析

### 核心特性

- 🧬 **全长转座子检测**: 使用动态边界调整方法，准确识别全长转座子
- 🎯 **多种 TE 类型支持**: LTR、TIR、Helitron、Non-LTR
- 📊 **泛基因组分析**: 构建泛基因组转座子库，进行群体比较
- 🔄 **断点续跑**: 支持从中断点恢复分析
- 🛡️ **容器化部署**: Singularity 容器确保环境一致性
- 📝 **详细日志**: 完整的分析过程记录和错误追踪

---

## ✨ 主要特性 | Key Features

### HiTE 单基因组分析

- ✅ **动态边界调整**: 高精度的全长 TE 检测
- ✅ **基因组注释**: 使用检测到的 TE 库注释基因组
- ✅ **蛋白结构域预测**: 预测 TE 中的保守蛋白结构域
- ✅ **嵌套 TE 处理**: 自动移除嵌套转座子
- ✅ **植物/动物模式**: 针对不同生物类型优化

### panHiTE 泛基因组分析

- ✅ **泛基因组 TE 库构建**: 整合多个基因组的转座子
- ✅ **低拷贝 TE 恢复**: 重新比对低拷贝 TE 到泛基因组
- ✅ **基因表达分析**: 整合 RNA-seq 数据进行 TIDELs 分析
- ✅ **批量处理**: 高效处理数百个基因组

---

## 🚀 快速开始 | Quick Start

### HiTE 单基因组分析

#### 基本用法

```bash
# 分析植物基因组
biopytools hite \
    -i genome.fa \
    -o hite_results \
    -t 12

# 分析动物基因组
biopytools hite \
    -g animal_genome.fa \
    --no-plant \
    -o results \
    -t 12
```

#### 完整分析（包含注释）

```bash
biopytools hite \
    -g plant_genome.fa \
    -o annotated_results \
    -t 12 \
    --plant \
    --annotate \
    --domain
```

#### 恢复中断的分析

```bash
biopytools hite \
    -i genome.fa \
    -o hite_results \
    -t 12 \
    --recover
```

---

### panHiTE 泛基因组分析

#### 基本用法（仅 TE 检测）

```bash
# 准备基因组列表文件 genome_list.txt
# 格式: 每行一个基因组文件名
genome1.fa
genome2.fa
genome3.fa

# 运行 panHiTE
biopytools panhite \
    -p pan_genomes/ \
    -i genome_list.txt \
    -o panhite_results \
    -t 12
```

#### 包含基因注释

```bash
# 准备基因组列表文件 genome_list_with_genes.txt
# 格式: 基因组文件    基因注释文件
genome1.fa    genome1.gff
genome2.fa    genome2.gff
genome3.fa    genome3.gff

# 运行 panHiTE
biopytools panhite \
    -p pan_genomes/ \
    -i genome_list_with_genes.txt \
    --genes-dir annotations/ \
    -o panhite_results \
    -t 12
```

#### 包含 RNA-seq 分析（TIDELs）

```bash
# 准备基因组列表文件 genome_list_full.txt
# 格式: 基因组    基因注释    是否配对    RNA-seq文件
genome1.fa    genome1.gff    1    rna1_1.fq.gz    rna1_2.fq.gz
genome2.fa    genome2.gff    1    rna2_1.fq.gz    rna2_2.fq.gz

# 运行 panHiTE
biopytools panhite \
    -p pan_genomes/ \
    -i genome_list_full.txt \
    --genes-dir annotations/ \
    --rna-dir rnaseq/ \
    -o panhite_results \
    -t 12
```

#### 仅生成 panTE 库（跳过分析）

```bash
biopytools panhite \
    -p pan_genomes/ \
    -i genome_list.txt \
    -o panhite_results \
    --skip-analyze
```

---

## 📋 参数说明 | Parameters

### HiTE 参数 (`biopytools hite`)

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 基因组 FASTA 文件路径 | `-i genome.fa` |

#### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./hite_output` | 📁 输出目录路径 |
| `-t, --threads` | `12` | 🔧 线程数 |

#### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--plant/--no-plant` | `--plant` | 🌱 是否为植物基因组 |
| `--annotate` | `False` | 📝 是否注释基因组 |
| `--domain` | `False` | 🔬 是否预测蛋白结构域 |
| `--te-type` | `all` | 🎯 TE 类型 [ltr\|tir\|helitron\|non-ltr\|all] |
| `--recover` | `False` | 🔄 是否启用断点续跑 |

#### 高级参数 | Advanced Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--chunk-size` | `400` | 基因组分块大小 (MB) |
| `--miu` | `1.3e-8` | 中性突变率 |
| `--min-te-len` | `80` | 最小 TE 长度 (bp) |

#### 容器配置 | Container Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--singularity-cmd` | `~/miniforge3/envs/...` | Singularity 命令路径 |
| `--sif-file` | `/share/.../hite_3.3.3.sif` | HiTE SIF 镜像路径 |

---

### panHiTE 参数 (`biopytools panhite`)

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-p, --pan-genomes-dir` | 泛基因组目录路径 | `-p pan_genomes/` |
| `-i, --input` | 基因组列表文件路径 | `-i genome_list.txt` |

#### 可选路径 | Optional Paths

| 参数 | 描述 | 示例 |
|------|------|------|
| `--genes-dir` | 基因注释文件目录 | `--genes-dir annotations/` |
| `--rna-dir` | RNA-seq 数据目录 | `--rna-dir rnaseq/` |

#### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./panhite_output` | 📁 输出目录路径 |
| `-t, --threads` | `12` | 🔧 线程数 |

#### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--te-type` | `all` | 🎯 TE 类型 |
| `--miu` | `1.3e-8` | 中性突变率 |
| `--skip-analyze` | `False` | 仅生成 panTE 库 |
| `--recover` | `False` | 断点续跑 |
| `--debug` | `False` | 调试模式 |

---

## 📂 输出文件 | Output Files

### HiTE 输出文件

#### 主要输出文件

| 文件名 | 描述 |
|--------|------|
| `confident_TE.cons.fa` | 高置信度 TE 序列（推荐使用） |
| `confident_ltr_cut.fa.cons` | LTR 转座子 |
| `confident_tir_*.fa` | TIR 转座子 |
| `confident_helitron_*.fa` | Helitron 转座子 |
| `confident_non_ltr_*.fa` | Non-LTR 转座子 |
| `all_TE.fa` | 所有 TE 序列 |
| `low_confident_TE.cons.fa` | 低置信度 TE 候选 |

#### 注释输出文件（使用 `--annotate` 时）

| 文件名 | 描述 |
|--------|------|
| `HiTE.out` | HiTE 注释结果 |
| `HiTE.gff` | GFF 格式注释（可用于 IGV 可视化） |
| `HiTE.tbl` | TE 统计信息表 |

#### 结构域输出文件（使用 `--domain` 时）

| 文件名 | 描述 |
|--------|------|
| `confident_TE.cons.fa.domain` | TE 蛋白结构域预测结果 |

#### 结果汇总文件

| 文件名 | 描述 |
|--------|------|
| `hite_summary.json` | JSON 格式的结果汇总 |
| `hite_*.log` | 详细的分析日志 |

---

### panHiTE 输出文件

#### 主要输出文件

| 文件名 | 描述 |
|--------|------|
| `panTE_library.fa` | 泛基因组 TE 库 |
| `panTE_annotation.gff` | panTE 注释文件 |
| `*_stats.txt` | 统计信息 |
| `*_cluster.fa` | 聚类结果 |
| `presence_absence*.txt` | 存在-缺失矩阵 |

#### 结果汇总文件

| 文件名 | 描述 |
|--------|------|
| `panhite_summary.json` | JSON 格式的结果汇总 |
| `panhite_*.log` | 详细的分析日志 |

---

## 🛠️ 环境要求 | Requirements

### 系统要求

- **操作系统**: Linux (推荐 CentOS 7+, Ubuntu 16.04+)
- **内存**: 推荐 128 GB RAM
- **CPU**: 推荐 40+ 核心数
- **磁盘空间**: 根据基因组大小，至少 100GB 可用空间

### 软件依赖

- **Singularity**: 3.8+ (通过 conda 安装)
- **HiTE SIF 镜像**: `/share/org/YZWL/yzwl_lixg/software/singularity/hite_3.3.3.sif`

### 安装 Singularity

```bash
# 使用 conda 安装
conda create -n singularity_v.3.8.7 -c conda-forge singularity=3.8.7
conda activate singularity_v.3.8.7

# 验证安装
singularity --version
```

---

## 🔧 故障排除 | Troubleshooting

### 常见问题

#### 1. Singularity 未找到

**错误信息**:
```
Singularity不存在|Singularity not found: ...
```

**解决方案**:
```bash
# 检查 Singularity 安装
which singularity

# 或指定 Singularity 路径
biopytools hite -i genome.fa --singularity-cmd /path/to/singularity
```

#### 2. SIF 文件未找到

**错误信息**:
```
HiTE SIF文件不存在|HiTE SIF file not found: ...
```

**解决方案**:
```bash
# 检查 SIF 文件
ls -lh /share/org/YZWL/yzwl_lixg/software/singularity/hite_3.3.3.sif

# 或指定 SIF 文件路径
biopytools hite -i genome.fa --sif-file /path/to/hite.sif
```

#### 3. 容器启动失败

**错误信息**:
```
FATAL:   While initializing...
```

**解决方案**:
```bash
# 检查 Singularity 配置
singularity exec hite_3.3.3.sif python --version

# 尝试使用 conda 安装的 Singularity
~/miniforge3/envs/singularity_v.3.8.7/bin/singularity --version
```

#### 4. 权限问题

**错误信息**:
```
Permission denied
```

**解决方案**:
```bash
# 检查输出目录权限
ls -ld output_dir/

# 修改权限
chmod 755 output_dir/
```

#### 5. 内存不足

**错误信息**:
```
Killed
```

**解决方案**:
```bash
# 减少线程数
biopytools hite -i genome.fa -t 20

# 或增加系统交换空间
```

---

## 📚 更多示例 | More Examples

### 示例 1: 检测特定类型 TE

```bash
# 只检测 LTR 转座子
biopytools hite \
    -i genome.fa \
    -o ltr_results \
    -t 12 \
    --te-type ltr

# 只检测 TIR 转座子
biopytools hite \
    -i genome.fa \
    -o tir_results \
    -t 12 \
    --te-type tir
```

### 示例 2: 调整参数优化

```bash
# 小基因组
biopytools hite \
    -g small_genome.fa \
    -o results \
    -t 12 \
    --chunk-size 200 \
    --min-te-len 100

# 大基因组
biopytools hite \
    -g large_genome.fa \
    -o results \
    -t 24 \
    --chunk-size 500 \
    --miu 1.5e-8
```

### 示例 3: 批量分析多个基因组

```bash
# 使用脚本批量分析
for genome in genomes/*.fa; do
    name=$(basename $genome .fa)
    biopytools hite \
        -g $genome \
        -o results/${name} \
        -t 12 \
        --plant
done
```

### 示例 4: panHiTE 仅生成 TE 库

```bash
# 快速生成 panTE 库，不进行深入分析
biopytools panhite \
    -p pan_genomes/ \
    -i genome_list.txt \
    -o pante_library \
    -t 12 \
    --skip-analyze
```

---

## 📖 参考文献 | References

### HiTE

Hu, K., Ni, P., Xu, M. et al. HiTE: a fast and accurate dynamic boundary adjustment approach for full-length transposon detection and annotation. *Nat Commun* **15**, 5573 (2024).

https://doi.org/10.1038/s41467-024-49912-8

### panHiTE

Hu, K. et al. HiTE: accurate transposable element identification and pan-genome analysis. *bioRxiv* preprint (2025).

https://doi.org/10.1101/2025.02.15.638472v1

### 官方文档

- HiTE GitHub: https://github.com/CSU-KangHu/HiTE
- HiTE Wiki: https://github.com/CSU-KangHu/HiTE/wiki

---

## 📞 技术支持 | Support

如有问题或建议，请联系：
- BioPyTools 项目组
- GitHub Issues: https://github.com/your-org/biopytools/issues

---

## 📝 版本历史 | Version History

- **v1.0.0** (2026-01-18): 初始版本
  - ✅ HiTE 单基因组分析
  - ✅ panHiTE 泛基因组分析
  - ✅ Singularity 容器化部署
  - ✅ 完整的日志和结果处理
