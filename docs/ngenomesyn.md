# 🧬 NGenomeSyn 基因组共线性可视化工具

**使用NGenomeSyn进行基因组共线性可视化和结构变异分析的专业工具**

## 📖 功能概述 | Overview

NGenomeSyn工具基于NGenomeSyn软件，提供基因组共线性可视化功能。支持：
- **两种比对工具**：Minimap2（快速）和MUMmer（精确）
- **两种分析模式**：基础可视化和SyRI深度结构变异分析
- **自动化流程**：从比对到可视化的完整pipeline

## ✨ 主要特性 | Key Features

- **🔧 双比对引擎**：支持Minimap2和MUMmer两种比对工具
- **🧬 结构变异识别**：SyRI模式自动识别倒位、易位、重复等结构变异
- **🎨 智能可视化**：结构变异通过坐标方向编码，NGenomeSyn自动识别并着色
- **📋 简化的样本文件**：只需两列（基因组文件  基因组名称）
- **🔬 完整工作流**：自动完成比对→过滤→SyRI分析→转换→可视化
- **📝 标准日志**：遵循项目日志规范

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# Minimap2模式（快速，默认）
biopytools ngenomesyn \
    -s samples.txt \
    -o ./output

# MUMmer模式（精确）
biopytools ngenomesyn \
    -s samples.txt \
    -o ./output \
    --aligner mummer

# SyRI结构变异分析（推荐）
biopytools ngenomesyn \
    -s samples.txt \
    -o ./output \
    --use-syri
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-s, --sample-map` | 样本映射文件 | `-s samples.txt` |
| `-o, --output` | 输出目录 | `-o ./output` |

### 样本文件格式 | Sample File Format

**格式**：
```
genome_file.fa    GenomeName
reference.fa     Ref
query.fa         Query
```

**示例**：
```
/share/data/ref.fa    Reference
/share/data/query.fa  Sample1
```

### 分析模式参数 | Analysis Mode Parameters

| 参数 | 描述 | 默认值 | 可选值 |
|------|------|--------|--------|
| `--aligner` | 比对工具选择 | minimap2 | minimap2, mummer |
| `--use-syri` | 启用SyRI结构变异分析 | False | - |

### Minimap2参数 | Minimap2 Parameters

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--preset` | Minimap2预设模式 | asm5 |

### MUMmer参数 | MUMmer Parameters

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--mummer-match` | MUMmer匹配选项 | --maxmatch |
| `--mummer-min-len` | MUMmer最小匹配长度 | 100 |
| `--delta-identity` | Delta过滤最小一致性(%) | 90.0 |
| `--delta-min-len` | Delta过滤最小比对长度 | 100 |

### 通用参数 | General Parameters

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--min-length` | LINK文件最小比对长度(bp) | 5000 |

### 染色体过滤 | Chromosome Filtering

| 参数 | 描述 | 示例 |
|------|------|------|
| `--chromosome` | 指定染色体 | `"1,2,3"` 或 `"1-5"` |

### 可视化参数 | Visualization Parameters

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--format` | 输出格式 | svg, png |
| `--threads` | 线程数 | 12 |

## 💡 使用示例 | Usage Examples

### 示例1：Minimap2快速可视化

```bash
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./minimap2_output
```

**输出**：
- `.len`文件（基因组长度信息）
- `.paf`文件（Minimap2比对结果）
- `.link`文件（NGenomeSyn输入）
- `.svg`和`.png`可视化图

### 示例2：MUMmer精确比对

```bash
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./mummer_output \
    --aligner mummer
```

**输出**：
- `.len`文件
- `.delta`文件（MUMmer比对结果）
- `.filtered.delta`（过滤后的delta）
- `.coords`文件（坐标格式）
- `.link`文件
- 可视化图

### 示例3：SyRI结构变异分析（推荐）

```bash
# Minimap2 + SyRI
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./syri_output \
    --use-syri

# MUMmer + SyRI
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./syri_output \
    --aligner mummer \
    --use-syri
```

**输出**：
- 比对文件（`.paf`或`.coords`）
- `*_syri_results/`目录：
  - `syri.out`（结构变异详情）
  - `syri.vcf`（VCF格式）
- `.link`文件（包含结构变异信息）
- 彩色可视化图

**SyRI识别的结构变异类型**：
- ✅ **INV**：倒位（Inversion）
- ✅ **TRA**：易位（Translocation）
- ✅ **DUP**：重复（Duplication）
- ✅ **INVTR**：倒位+易位
- ✅ **INVDP**：倒位+重复
- ✅ 及其他复杂结构变异

### 示例4：指定染色体分析

```bash
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./chr1_3_output \
    --chromosome "1,2,3"
```

### 示例5：自定义Minimap2参数

```bash
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./high_quality_output \
    --use-syri \
    --preset asm20 \
    --min-length 100000 \
    --threads 24
```

### 示例6：自定义MUMmer参数

```bash
biopytools ngenomesyn \
    -s genome_samples.txt \
    -o ./custom_mummer_output \
    --aligner mummer \
    --use-syri \
    --mummer-match --mum \
    --mummer-min-len 200 \
    --delta-identity 95 \
    --delta-min-len 500 \
    --min-length 10000
```

## 📁 输出文件 | Output Files

### Minimap2模式（不含SyRI）

```
output_dir/
├── Ref.len                    # 参考基因组长度文件
├── Query.len                  # 查询基因组长度文件
├── Ref_vs_Query.paf          # Minimap2比对结果
├── Ref_vs_Query.link         # LINK文件（NGenomeSyn输入）
├── ngenomesyn.conf           # NGenomeSyn配置文件
├── ngenomesyn.svg            # SVG可视化图
└── ngenomesyn.png            # PNG可视化图
```

### Minimap2 + SyRI模式

```
output_dir/
├── Ref.len                    # 参考基因组长度文件
├── Query.len                  # 查询基因组长度文件
├── Ref_vs_Query.paf          # Minimap2比对结果
├── Ref_vs_Query.link         # LINK文件（含结构变异信息）
├── Ref_vs_Query_syri_results/ # SyRI分析结果目录
│   ├── syri.out              # 结构变异详情
│   ├── syri.vcf              # VCF格式变异
│   └── syri.log              # SyRI运行日志
├── ngenomesyn.conf           # NGenomeSyn配置文件
├── ngenomesyn.svg            # 彩色可视化图
└── ngenomesyn.png            # 彩色可视化图
```

### MUMmer模式（不含SyRI）

```
output_dir/
├── Ref.len                    # 参考基因组长度文件
├── Query.len                  # 查询基因组长度文件
├── Ref_vs_Query.delta        # MUMmer比对结果
├── Ref_vs_Query.filtered.delta # 过滤后的delta
├── Ref_vs_Query.coords       # 坐标格式
├── Ref_vs_Query.link         # LINK文件
├── ngenomesyn.conf           # NGenomeSyn配置文件
├── ngenomesyn.svg            # 可视化图
└── ngenomesyn.png            # 可视化图
```

### MUMmer + SyRI模式

```
output_dir/
├── Ref.len
├── Query.len
├── Ref_vs_Query.delta
├── Ref_vs_Query.filtered.delta
├── Ref_vs_Query.coords
├── Ref_vs_Query.link         # 含结构变异信息
├── Ref_vs_Query_syri_results/
│   ├── syri.out
│   ├── syri.vcf
│   └── syri.log
├── ngenomesyn.conf
├── ngenomesyn.svg
└── ngenomesyn.png
```

## 🎨 可视化说明 | Visualization Guide

### 结构变异在LINK文件中的编码

**重要**：Syri2Link转换时，结构变异信息通过**坐标方向**编码，而非显式标签：

- **正向共线性**：`chrA startA endA chrB startB endB`
  - 特征：startB < endB（正向坐标）
  - 颜色：通常为绿色

- **倒位（INV）**：`chrA startA endA chrB endB startB`
  - 特征：startB > endB（坐标反向）
  - 颜色：通常为红色

- **易位（TRA）**：涉及不同染色体的共线性关系
  - 颜色：通常为蓝色

NGenomeSyn会自动识别这些坐标模式，并使用不同颜色显示。

### SyRI输出文件格式

SyRI的`syri.out`文件包含11列：
```
chrA  startA  endA  seqA  seqB  chrB  startB  endB  ID  groupID  type
```

- **type列**：结构变异类型
  - `SYN`：正向共线性
  - `INV`：倒位
  - `TRA`：易位
  - `INVDP`：倒位+重复
  - `INVTR`：倒位+易位
  - `DUP`：重复
  - `NOTAL`：无共线性
  - `SNP/INS/DEL`：小变异（Syri2Link会过滤掉）

Syri2Link转换时会：
- 过滤NOTAL和小变异
- 过滤长度 < min_alignment_length的比对
- 保留SYN、INV、TRA等大尺度结构变异
- 将坐标信息写入LINK文件

## ⚙️ 系统要求 | System Requirements

### 依赖软件 | Dependencies

**必需软件**：
- **Python** 3.7+
- **NGenomeSyn** 1.43+
- **Perl** (GetTwoGenomeSyn.pl脚本需要)

**比对工具（二选一）**：
- **Minimap2** 2.0+ 或
- **MUMmer** 4.0+ (包含mummer, delta-filter, show-coords)

**可选软件**：
- **SyRI** (使用--use-syri时需要)

**Python包**：
- `click`
- `biopython`

## 📚 工作流程 | Workflow

### Minimap2流程（不含SyRI）

```
1. 读取样本文件
2. 计算染色体长度 → 生成.len文件
3. Minimap2比对 → 生成.paf文件
4. PAF转LINK → GetTwoGenomeSyn.pl Paf2Link
5. 生成NGenomeSyn配置
6. 运行NGenomeSyn → 生成可视化图
```

### Minimap2 + SyRI流程

```
1. 读取样本文件
2. 计算染色体长度 → 生成.len文件
3. Minimap2比对 → 生成.paf文件
4. SyRI分析：
   - 读取.paf
   - 识别结构变异
   - 生成syri.out和syri.vcf
5. SyRI结果转LINK → GetTwoGenomeSyn.pl Syri2Link
6. 生成NGenomeSyn配置
7. 运行NGenomeSyn → 生成彩色可视化图
```

### MUMmer流程（不含SyRI）

```
1. 读取样本文件
2. 计算染色体长度 → 生成.len文件
3. MUMmer比对 → 生成.delta文件
4. Delta过滤 → delta-filter
5. 转换为coords → show-coords
6. Coords转LINK → GetTwoGenomeSyn.pl Coords2Link
7. 生成NGenomeSyn配置
8. 运行NGenomeSyn → 生成可视化图
```

### MUMmer + SyRI流程

```
1. 读取样本文件
2. 计算染色体长度 → 生成.len文件
3. MUMmer比对 → 生成.delta文件
4. Delta过滤 → delta-filter
5. 转换为coords → show-coords
6. SyRI分析：
   - 读取.coords和.filtered.delta
   - 识别结构变异
   - 生成syri.out和syri.vcf
7. SyRI结果转LINK → GetTwoGenomeSyn.pl Syri2Link
8. 生成NGenomeSyn配置
9. 运行NGenomeSyn → 生成彩色可视化图
```

## ⚠️ 注意事项 | Important Notes

1. **样本文件格式**：必须是两列（基因组文件  基因组名称），用Tab分隔
2. **基因组文件**：支持FASTA和FASTA.GZ格式
3. **至少两个基因组**：共线性分析至少需要两个基因组
4. **SyRI模式**：
   - 需要更多时间和计算资源
   - 结果更丰富，能识别详细的结构变异
   - 推荐用于研究结构变异
5. **染色体编号**：染色体编号从1开始，连续编号
6. **MUMmer vs Minimap2**：
   - Minimap2：速度快，适合快速可视化
   - MUMmer：更精确，适合精细分析

## 📄 许可证 | License

MIT License

---

## 🔬 引用信息 | Citation

使用本工具时，请引用以下软件：

**NGenomeSyn**:
```
Sun Haoran et al. NGenomeSyn: a three-dimensional genome visualization tool for multiple genomes.
Bioinformatics, 2023.
```

**SyRI** (如果使用):
```
Akhtar et al. SyRI: finding genomic rearrangements and structural variants from whole-genome assemblies.
Genome Biology, 2019.
```

**Minimap2** (如果使用):
```
Li H. Minimap2: pairwise alignment for nucleotide sequences.
Bioinformatics, 2018.
```

**MUMmer** (如果使用):
```
Kurtz S et al. Versatile and open software for comparing large genomes.
Genome Biology, 2004.
```
