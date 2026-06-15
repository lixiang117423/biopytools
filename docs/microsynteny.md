# 📝 Microsynteny 微观共线性分析模块

**基于JCVI的自动化微观共线性分析和可视化工具 | Automated Microsynteny Analysis and Visualization Tool Based on JCVI**

## 📖 功能概述 | Overview

Microsynteny 微观共线性分析模块是一个专业的基因级别共线性分析和可视化工具，基于JCVI (Python Comparative Genomics Toolkit) 构建，提供从基因组数据准备到最终微观共线性图绘制的完整自动化流程。支持多物种比较、智能区块提取、灵活的参数配置和专业的出版级可视化，适用于比较基因组学、进化生物学和基因家族研究。

## ✨ 主要特性 | Key Features

- **🔄 完整自动化流程**: GFF转换→序列比对→共线性计算→区块提取→自动绘图五步骤全自动化
- **🎯 断点续传功能**: 支持步骤级别的断点续传，已完成的步骤自动跳过
- **🛡️ 智能数据预处理**: 自动GFF→BED转换、文件格式验证和物种识别
- **🔍 灵活共线性分析**: 支持两两物种比较、可配置的共线性分数阈值
- **📊 智能区块提取**: 根据目标基因自动提取相关共线性区块
- **🎨 自动布局生成**: 动态生成最优的物种排列布局
- **📂 多物种支持**: 支持2个及以上物种的微观共线性比较
- **⚙️ 高度可配置**: 自定义JCVI路径、线程数、延伸基因数等参数
- **📊 详细日志记录**: 完整的处理过程日志和错误追踪
- **🚀 专业可视化**: 生成出版级的PDF/SVG/PNG格式共线性图

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 准备数据文件夹
# genome_folder/
# ├── A.fa, A.gff
# ├── B.fa, B.gff
# └── C.fa, C.gff

# 准备目标基因列表文件
# genes.txt (两列：物种ID  基因ID)
# A    GeneA_001
# B    GeneB_001
# C    GeneC_001

# 运行微观共线性分析
biopytools microsynteny \
    -i ./genome_folder \
    -g genes.txt \
    -o output_results
```

### 高级用法 | Advanced Usage

```bash
# 自定义JCVI路径和分析参数
biopytools microsynteny \
    -i ./genome_folder \
    -g genes.txt \
    -j ~/miniforge3/envs/jcvi_v.1.5.7 \
    --cscore 0.95 \
    --extend-genes 50 \
    -t 24 \
    -o custom_output

# 只运行特定步骤
biopytools microsynteny \
    -i ./genome_folder \
    -g genes.txt \
    --step 1 \
    -o step1_preprocess
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --genome-folder` | 基因组文件夹路径（包含A.fa, A.gff等文件）| `-i ./genome_data` |
| `-g, --gene-list` | 目标基因列表文件（两列：物种ID 基因ID）| `-g genes.txt` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-j, --jcvi-path` | `~/miniforge3/envs/jcvi_v.1.5.7` | 🛠️ JCVI环境路径 |
| `-o, --output-dir` | `./microsynteny_output` | 📁 输出目录路径 |
| `-t, --threads` | `12` | ⚙️ 线程数 |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--extend-genes` | `30` | 📏 延伸基因数（在目标基因每侧延伸的基因数）|
| `--cscore` | `0.99` | 🎯 共线性分数阈值（0-1，值越高越严格）|

### 处理控制选项 | Processing Control Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--step` | `全部` | 🎯 只运行指定步骤 (1/2/3/4) |
| `--log-level` | `INFO` | 📊 日志级别 (DEBUG/INFO/WARNING/ERROR/CRITICAL) |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 🔄 数据预处理 | GFF转BED、提取CDS序列、合并BED文件 |
| **2** | 🔍 共线性分析 | 两两物种序列比对、共线性计算 |
| **3** | 📦 区块提取 | 提取微观共线性区块、目标基因相关区块 |
| **4** | 🎨 绘图 | 自动生成layout文件、绘制共线性图 |

## 📁 输入文件格式 | Input File Formats

### 基因组文件夹 | Genome Folder

基因组文件夹应包含以下格式的文件：

```
genome_folder/
├── A.fa      # 物种A的基因组序列（FASTA格式）
├── A.gff     # 物种A的基因注释（GFF3格式）
├── B.fa      # 物种B的基因组序列
├── B.gff     # 物种B的基因注释
└── ...
```

**文件命名规则：**
- 序列文件：`{物种ID}.fa` 或 `{物种ID}.fasta`
- 注释文件：`{物种ID}.gff` 或 `{物种ID}.gff3`

### 基因组序列文件 | Genome Sequence File (FASTA)

标准FASTA格式的基因组序列：

```fasta
>chr1
ATCGATCGATCGATCGATCGATCGATCG...
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### 基因注释文件 | Gene Annotation File (GFF3)

标准GFF3格式的基因注释：

```gff3
##gff-version 3
##sequence-region chr1 1 10000000
chr1    RefSeq gene    1000   5000   .   +   .   ID=gene001;Name=GeneA
chr1    RefSeq mRNA    1000   5000   .   +   .   ID=mrna001;Parent=gene001
chr1    RefSeq CDS     1200   1400   .   +   0   ID=cds001;Parent=mrna001
```

### 目标基因列表文件 | Target Gene List File

两列文本文件（制表符或空格分隔）：

```text
A    GeneA_001
A    GeneA_002
B    GeneB_001
C    GeneC_001
```

**格式要求：**
- 第一列：物种ID（与文件名前缀一致）
- 第二列：基因ID（与GFF文件中的ID一致）
- 支持制表符或空格分隔
- 支持注释行（以#开头）

## 💡 使用示例 | Usage Examples

### 示例1：双物种NLR基因簇分析 | Example 1: Two-species NLR gene cluster analysis

```bash
# 分析两个物种的NLR基因簇共线性
biopytools microsynteny \
    -i ./nlr_genomes \
    -g nlr_genes.txt \
    -o nlr_synteny \
    --extend-genes 20 \
    --cscore 0.99
```

### 示例2：多物种全基因组比较 | Example 2: Multi-species whole-genome comparison

```bash
# 五个物种的全基因组微观共线性比较
biopytools microsynteny \
    -i ./plant_genomes \
    -g target_genes.txt \
    -o multi_species_synteny \
    -t 24 \
    --log-level DEBUG
```

### 示例3：分步骤处理 | Example 3: Step-by-step processing

```bash
# 第一步：数据预处理
biopytools microsynteny \
    -i ./genome_data \
    -g genes.txt \
    --step 1 \
    -o preprocess_output

# 第二步：共线性分析（如果第一步成功完成）
biopytools microsynteny \
    -i ./genome_data \
    -g genes.txt \
    --step 2 \
    -o preprocess_output

# 第三步：区块提取
biopytools microsynteny \
    -i ./genome_data \
    -g genes.txt \
    --step 3 \
    -o preprocess_output

# 第四步：绘图
biopytools microsynteny \
    -i ./genome_data \
    -g genes.txt \
    --step 4 \
    -o preprocess_output
```

### 示例4：使用自定义JCVI环境 | Example 4: Using custom JCVI environment

```bash
# 指定自定义JCVI环境路径
biopytools microsynteny \
    -i ./genome_data \
    -g genes.txt \
    -j ~/miniforge3/envs/jcvi_custom \
    -o custom_env_output
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
microsynteny_output/
├── 1_preprocess/              # 数据预处理结果
│   ├── A.bed                  # 物种A的BED文件
│   ├── B.bed                  # 物种B的BED文件
│   └── all_species.bed        # 合并的BED文件
├── 2_synteny/                 # 共线性分析结果
│   ├── A.B.last               # LAST比对结果
│   ├── A.B.anchors            # 共线性锚点
│   └── A.B.lifted.anchors     # 扩展的共线性区块
├── 3_blocks/                  # 区块提取结果
│   ├── A.B.blocks             # 微观共线性区块
│   └── merged.blocks          # 合并的区块文件
├── 4_plot/                    # 绘图结果
│   ├── synteny.layout         # 自动生成的layout文件
│   ├── microsynteny.pdf       # PDF格式共线性图
│   ├── microsynteny.svg       # SVG格式共线性图（可编辑）
│   └── microsynteny.png       # PNG格式共线性图
├── logs/                      # 日志文件
│   └── microsynteny.log       # 完整日志记录
├── .step_1.done               # 步骤1完成标记
├── .step_2.done               # 步骤2完成标记
├── .step_3.done               # 步骤3完成标记
└── .step_4.done               # 步骤4完成标记
```

### 关键输出文件说明 | Key Output Files Description

**预处理文件 (1_preprocess/)：**
- `*.bed`: 基因位置信息（染色体、起止位置、方向）
- `all_species.bed`: 所有物种的合并BED文件

**共线性文件 (2_synteny/)：**
- `*.last`: LAST序列比对原始结果
- `*.anchors`: 高质量共线性锚点
- `*.lifted.anchors`: 扩展后的完整共线性区块

**区块文件 (3_blocks/)：**
- `*.blocks`: 微观共线性区块文件
- `merged.blocks`: 合并后的区块文件（用于绘图）

**可视化文件 (4_plot/)：**
- `*.pdf`: 矢量PDF图（适合论文发表）
- `*.svg`: 可编辑SVG图（可用Inkscape/AI编辑）
- `*.png`: 位图PNG图（适合快速预览）

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **JCVI** (版本 1.5.7 或更新)
  - 安装方式: `conda install -c bioconda jcvi` 或 `pip install jcvi`
  - Python环境: 需要3.7+
  - 下载地址: https://github.com/tanghaibao/jcvi

- **LAST** (序列比对工具)
  - JCVI会自动调用LAST进行序列比对
  - 下载地址: https://gitlab.com/mcfrith/last

### Python包依赖 | Python Package Dependencies

- `click`: 命令行界面框架
- `pathlib`: 路径处理
- `logging`: 日志记录
- `dataclasses`: 配置数据类
- `subprocess`: 命令执行

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐4核以上，多物种比较推荐8核+）
- **RAM**: 最少4GB（大基因组推荐16GB以上）
- **存储**: 预留基因组文件大小10倍的磁盘空间
- **网络**: 如需下载依赖软件，建议稳定网络

## ⚠️ 注意事项 | Important Notes

1. **文件命名规范**: 基因组文件必须按 `{物种ID}.fa` 和 `{物种ID}.gff` 格式命名
2. **基因ID一致性**: 基因列表文件中的基因ID必须与GFF文件中的ID完全一致
3. **JCVI环境**: 确保JCVI环境正确安装并可用
4. **内存需求**: 大基因组或多物种比较需要较多内存
5. **运行时间**: 序列比对和共线性计算可能需要较长时间
6. **断点续传**: 已完成的步骤会自动跳过，无需手动删除中间文件

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "找不到JCVI Python" 错误**
```bash
# 检查JCVI安装
ls ~/miniforge3/envs/jcvi_v.1.5.7/bin/python

# 如路径不同，使用-j参数指定正确路径
biopytools microsynteny ... -j /path/to/jcvi
```

**Q: "基因组文件不存在" 错误**
```bash
# 检查文件命名是否正确
ls ./genome_folder/
# 应该看到: A.fa, A.gff, B.fa, B.gff 等

# 检查genes.txt中的物种ID是否与文件名前缀一致
cat genes.txt
```

**Q: "基因ID在GFF文件中找不到" 错误**
```bash
# 检查GFF文件中的基因ID格式
grep "GeneA_001" ./genome_folder/A.gff

# 确保genes.txt中的ID与GFF中的ID完全一致（包括大小写）
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools microsynteny ... -t 4

# 或分步运行
biopytools microsynteny ... --step 1
biopytools microsynteny ... --step 2
```

**Q: 共线性分析失败（LAST错误）**
```bash
# 检查LAST是否安装
which lastal lastdb

# 如果未安装，安装LAST
conda install -c bioconda last
```

## 📚 相关资源 | Related Resources

- [JCVI官方文档](https://github.com/tanghaibao/jcvi)
- [JCVI Wiki](https://github.com/tanghaibao/jcvi/wiki)
- [LAST序列比对工具](https://gitlab.com/mcfrith/last)
- [GFF3格式规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [比较基因组学最佳实践](https://www.nature.com/articles/s41592-021-01158-9)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: JCVI和LAST软件本身有各自的许可证，请遵守相关软件的使用条款。

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用JCVI相关文献：

```
Tang, H., Zhang, X., Miao, R. et al.
JCVI: a graphical toolkit for comparative genomics.
BMC Bioinformatics 22, 361 (2021).
https://doi.org/10.1186/s12859-021-04315-x
```

---

## 💻 高级技巧 | Advanced Tips

### 自定义布局文件 | Custom Layout File

如果自动生成的布局不理想，可以手动编辑 `synteny.layout` 文件：

```text
# x,   y, rotation,   ha,     va,   color, ratio, label
0.5, 0.6,        0, center,    top,      ,     1, Species A
0.5, 0.4,        0, center, bottom,      ,    .5, Species B
0.5, 0.2,        0, center, bottom,      ,    .5, Species C
# edges
e, 0, 1
e, 0, 2
```

### 批量处理多个基因家族 | Batch Processing Multiple Gene Families

```bash
# 创建多个genes.txt文件
# nlr_genes.txt, disease_genes.txt, tf_genes.txt

# 批量运行
for genes in nlr disease tf; do
    biopytools microsynteny \
        -i ./genome_data \
        -g ${genes}_genes.txt \
        -o ${genes}_synteny
done
```

### 结果后期处理 | Post-processing

```bash
# 使用Inkscape或Adobe Illustrator编辑SVG图
inkscape microsynteny_output/4_plot/microsynteny.svg

# 高亮特定基因、调整颜色、添加注释等
```
