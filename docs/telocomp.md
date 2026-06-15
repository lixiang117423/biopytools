# TeloComp 端粒鉴定分析模块

**端粒序列识别与过滤工具 | Telomere Sequence Identification and Filtering Tool**

## 功能概述 | Overview

TeloComp 端粒鉴定分析模块是一个基于 ONT/HiFi 长读长数据的端粒序列识别和过滤工具。专注于从原始数据中提取含端粒重复序列的reads，为后续端粒组装和补全提供高质量的输入数据。

**当前模块功能范围**：
- ✅ **Filter_1**: 端粒reads识别和过滤（支持ONT/HiFi单种或两种数据）
- ✅ **Filter_2**: 端粒序列提取和修剪（需要同时提供ONT和HiFi数据）
- ✅ **结果统计**: 自动生成端粒reads数量统计和汇总报告
- ❌ **Assembly**: 端粒序列组装（**未实现**）
- ❌ **Complement**: 端粒补全和可视化图表生成（**未实现**）

**重要说明**：
- **Filter_1**: 支持单种或两种数据类型（ONT/HiFi），完成端粒reads过滤
- **Filter_2**: 需要同时提供 ONT 和 HiFi 数据。如果只有单种数据类型，会自动跳过Filter_2步骤
- 对于只有单种数据类型的场景，Filter_1 的输出（`filter1_ont.bam` / `filter1_hifi.bam`）已经包含了端粒鉴定结果
- **完整的可视化图表（telomere_plots）**需要使用TeloComp完整流程的Complement步骤，该步骤需要额外的WGS数据和NextPolish工具

## 主要特性 | Key Features

- **端粒序列鉴定**: 基于长读长数据识别含端粒重复序列的 reads
- **智能过滤系统**: 多级过滤策略，准确提取端粒序列
- **单/双数据类型支持**: 智能检测数据类型，自动调整分析流程
- **自动化流程**: 从 reads 比对到端粒分类的全流程自动化
- **灵活参数配置**: 支持植物、动物等不同物种的端粒序列模式
- **结果汇总统计**: 自动生成端粒鉴定结果汇总报告，展示端粒reads数量
- **详细日志记录**: 完整的分析过程日志和错误追踪
- **Conda环境管理**: 独立的 conda 环境，避免依赖冲突

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本端粒鉴定分析（同时使用ONT和HiFi数据）
biopytools telocomp \
    -g genome.fa \
    --ont ont_data.fastq.gz \
    --hifi hifi_data.fastq.gz \
    -o telomere_results

# 仅使用ONT数据
biopytools telocomp \
    -g genome.fa \
    --ont ont_data.fastq.gz \
    -o telomere_results

# 仅使用HiFi数据
biopytools telocomp \
    -g genome.fa \
    --hifi hifi_data.fastq.gz \
    -o telomere_results
```

### 高级用法 | Advanced Usage

```bash
# 自定义参数的端粒鉴定
biopytools telocomp \
    -g genome.fa \
    --ont ont_data.fastq.gz \
    --hifi hifi_data.fastq.gz \
    -m TTAGGG \
    -M 6 \
    -t 24 \
    -c 80 \
    -o results
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 基因组FASTA文件路径（系统会自动创建索引）| `-g genome.fa` |
| `-o, --output-dir` | 输出目录路径 | `-o telomere_output` |

### 输入数据配置 | Input Data Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--ont` | `None` | ONT长读长数据文件路径 |
| `--hifi` | `None` | HiFi长读长数据文件路径 |

**数据类型说明**：
- **必须至少提供 ONT 或 HiFi 数据其中一种**
- **仅提供一种数据类型**（如只有ONT或只有HiFi）：
  - Filter_1 正常执行，完成端粒reads过滤
  - Filter_2 自动跳过（因为需要两种数据类型）
  - 主要输出：`filter1_ont.bam` 或 `filter1_hifi.bam`（包含端粒reads）
- **同时提供两种数据类型**：
  - Filter_1 和 Filter_2 都正常执行
  - 完整的端粒提取、修剪和合并流程
  - 额外输出：`filter2_output/trim_L/` 和 `trim_R/` 目录

### 处理配置 | Processing Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `-c, --coverage` | `100` | 覆盖度参数 (0-100)，用于修剪 reads |
| `--motifs` | `[TTAGGG, CCCTAAA]` | 端粒重复序列模式列表 |

### 端粒配置 | Telomere Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --motif` | `CCCTAAA` | 端粒重复序列（植物默认） |
| `-M, --motif-num` | `7` | 端粒重复序列碱基数 |

**常见物种端粒序列**:
- 植物: `CCCTAAA` / `TTTAGGG`
- 动物: `TTAGGG` / `CCCTAA`
- 其他: 请根据物种调整

### 流程控制选项 | Pipeline Control Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--skip-filter` | `False` | 跳过Filter步骤 |
| `--no-visualization` | `False` | 跳过可视化步骤 |

## 输入文件格式 | Input File Formats

### 基因组文件 | Genome File

标准FASTA格式的基因组序列，**系统会自动创建索引**：

```bash
# 手动创建索引（可选，系统会自动创建）
samtools faidx genome.fa
```

```fasta
>chromosome1
ATCGATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### ONT/HiFi数据文件 | ONT/HiFi Data Files

标准 FASTQ 格式的长读长测序数据（支持 gzip 压缩）：

```fastq
@read_id
ATCGATCGATCGATCGATCGATCGATCG...
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
telomere_output/
├── telomere_positions.txt     # 端粒位置检测报告（重要！）
├── telomere_summary.txt       # 端粒鉴定结果汇总
├── telocomp.log               # 完整运行日志
├── filter1_ont.bam            # Filter_1过滤后的ONT BAM文件
├── filter1_hifi.bam           # Filter_1过滤后的HiFi BAM文件
└── filter2_output/            # Filter_2输出目录（如果运行）
    ├── trim_L/                # 左端端粒序列
    │   ├── seq1.fa
    │   ├── seq2.fa
    │   └── ...
    └── trim_R/                # 右端端粒序列
        ├── seq1.fa
        ├── seq2.fa
        └── ...
```

### 关键输出文件说明 | Key Output Files Description

**telomere_positions.txt** - 端粒位置检测报告（**最重要的输出**）
  - 告诉你在哪些染色体的哪些位置检测到了端粒
  - 左端端粒：染色体起始位置（1-5000 bp区域）
  - 右端端粒：染色体末端位置（染色体长度-5000 bp到末端）
  - 包含每条染色体左端和右端的端粒reads数量统计
  - 示例：
    ```
    染色体|Chromosome: Chr10 (长度|Length: 133,184,359 bp)
    --------------------------------------------------------------------------------
      左端端粒|Left telomere: 1-5000 bp区域
        ONT reads: 0
        HiFi reads: 3
        总计|Total: 3

      右端端粒|Right telomere: 133,179,359-133,184,359 bp区域
        ONT reads: 0
        HiFi reads: 5
        总计|Total: 5
    ```

**telomere_summary.txt** - 端粒鉴定结果汇总报告
  - 包含Filter_1和Filter_2的统计信息
  - 端粒reads数量统计
  - 关于可视化功能的说明

**filter1_ont.bam / filter1_hifi.bam** - Filter_1 步骤输出的 BAM 文件
  - 包含含端粒重复序列的比对结果
  - 可用于后续分析或可视化

**trim_L/**: 左端端粒序列目录（仅在Filter_2运行时生成）
  - 包含识别出的染色体左端端粒序列
  - FASTA 格式，每条序列一个文件

**trim_R/**: 右端端粒序列目录（仅在Filter_2运行时生成）
  - 包含识别出的染色体右端端粒序列
  - FASTA 格式，每条序列一个文件

## 分析流程 | Analysis Pipeline

### 完整流程步骤 | Complete Pipeline Steps

```
1. 检查基因组索引
   ↓
2. Filter_1: 过滤含端粒的reads
   - 使用 minimap2 比对 ONT/HiFi 数据到基因组
   - 使用 teloclip 处理 SAM 文件
   - 过滤含端粒重复序列的 reads
   - 输出 filter1_ont.bam 和 filter1_hifi.bam
   ↓
3. Filter_2: 检测和提取端粒序列（需要同时有ONT和HiFi数据）
   - 根据覆盖度修剪 reads
   - 分类左侧和右侧端粒序列
   - 输出到 trim_L 和 trim_R 目录
   ↓
4. 端粒位置分析（✅ 新增功能）
   - 分析 filter1 BAM 文件
   - 识别左端端粒（染色体起始位置附近）
   - 识别右端端粒（染色体末端位置附近）
   - 生成 telomere_positions.txt 报告
   ↓
5. 结果汇总
   - 统计端粒序列数量
   - 生成结果汇总报告
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

**Conda 环境**:
- Python 3.11+
- 所有依赖自动安装在 conda 环境中

**必须预先安装的软件**:
- **TeloComp** (版本 1.0.0)
  - 安装路径: `/share/org/YZWL/yzwl_lixg/software/telocomp/TeloComp-1.0.0/`
- **GenomeSyn**
  - 安装路径: `/share/org/YZWL/yzwl_lixg/software/GenomeSyn/GenomeSyn-main/GenomeSyn-1.2.7/`

**Conda 包**:
- samtools
- minimap2
- bwa
- flye
- Python: numpy, pandas, pysam, biopython, matplotlib, svglib

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐20核以上）
- **RAM**: 最少16GB（大基因组推荐64GB以上）
- **存储**: 预留基因组文件大小5倍的磁盘空间

## 注意事项 | Important Notes

1. **基因组索引**: 系统会自动创建索引（如果不存在）
2. **端粒序列**: 根据物种选择正确的端粒重复序列模式
3. **输入数据要求**:
   - 必须至少提供 ONT 或 HiFi 数据其中一种
   - 如果只有单种数据类型，Filter_2 会自动跳过，Filter_1 的输出仍然包含端粒鉴定结果
   - 同时提供两种数据类型可以获得更完整的端粒提取和修剪结果
4. **磁盘空间**: Filter_2 会生成大量中间文件，确保有足够磁盘空间
5. **运行时间**: 取决于数据量和基因组大小，可能需要数小时到数天

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "基因组索引文件不存在" 错误**
```bash
# 系统会自动创建索引，无需手动操作
# 如需手动创建，可运行：
samtools faidx genome.fa
```

**Q: "ModuleNotFoundError: No module named 'pysam'" 错误**
```bash
# 确保激活了正确的 conda 环境
conda activate telocomp

# 或重新安装依赖
cd /share/org/YZWL/yzwl_lixg/software/telocomp/TeloComp-1.0.0/bin
bash setup.sh
```

**Q: "至少需要提供 ONT 或 HiFi 数据" 错误**
```bash
# 必须提供至少一种长读长数据
biopytools telocomp -g genome.fa --ont data.fastq.gz -o output
```

**Q: 磁盘空间不足**
```bash
# 清理中间文件
rm -rf telomere_output/filter2_output/trim_L/*
rm -rf telomere_output/filter2_output/trim_R/*
```

## 使用示例 | Usage Examples

### 示例1：植物基因组端粒鉴定 | Example 1: Plant Genome Telomere Identification

```bash
# 使用ONT和HiFi数据鉴定植物端粒（默认CCCTAAA）
biopytools telocomp \
    -g plant_genome.fa \
    --ont plant_ont.fastq.gz \
    --hifi plant_hifi.fastq.gz \
    -o plant_telomere
```

### 示例2：动物基因组端粒鉴定 | Example 2: Animal Genome Telomere Identification

```bash
# 使用动物端粒序列（TTAGGG）
biopytools telocomp \
    -g animal_genome.fa \
    --ont animal_ont.fastq.gz \
    -m TTAGGG \
    -M 6 \
    -o animal_telomere
```

### 示例3：低覆盖度数据 | Example 3: Low Coverage Data

```bash
# 降低覆盖度阈值以获取更多端粒序列
biopytools telocomp \
    -g genome.fa \
    --hifi hifi.fastq.gz \
    -c 50 \
    -o low_coverage_telomere
```

## 与其他模块配合 | Integration with Other Modules

TeloComp 可以与以下模块配合使用：

1. **find_telomere**: 端粒序列查找模块
2. **hifiasm**: HiFi 数据基因组组装
3. **ragtag**: 基因组 scaffolding

## 相关资源 | Related Resources

- [TeloComp官方文档](https://github.com/lxie-0709/TeloComp)
- [端粒研究综述](https://www.nature.com/articles/s41576-020-0031-z)
- [T2T基因组组装指南](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8328731/)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用 TeloComp 相关文献：

```
TeloComp: An efficient integrated software package for telomere extraction and complementation
```

---

**最后更新 | Last Updated**: 2026-01-15
