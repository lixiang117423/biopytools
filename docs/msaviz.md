# MSA可视化模块

**自动多序列比对与可视化工具 | Auto MSA Alignment & Visualization Tool**

## 功能概述 | Overview

MSA可视化模块是一个集成了MAFFT多序列比对和pyMSAviz可视化功能的工具。默认自动运行MAFFT进行序列比对，然后生成美观的可视化图像。支持多种颜色方案、丰富的自定义选项和多种输出格式。

## 主要特性 | Key Features

- **自动比对**: 默认自动运行MAFFT进行多序列比对
- **丰富的颜色方案**: 15种内置颜色方案（Clustal、Zappo、Taylor等）
- **多种输出格式**: 支持PNG、JPG、SVG、PDF格式
- **灵活的布局控制**: 支持换行显示、自定义区域显示
- **Consensus序列**: 自动计算并显示consensus序列和一致性
- **标记和注释**: 支持在特定位点添加标记和文本注释
- **NJ树排序**: 支持基于NJ树的序列排序
- **序列类型检测**: 自动检测核苷酸或蛋白质序列

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 自动比对后可视化（默认）
biopytools msaviz -i sequences.fa -o output.png

# 输入已是比对结果，跳过比对
biopytools msaviz -i aligned.fa -o output.png --skip-align
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --infile` | 输入序列文件或MSA文件 | `-i sequences.fa` |
| `-o, --outfile` | 输出可视化文件 | `-o output.png` |

### 比对参数 | Alignment Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--skip-align` | `False` | 跳过MAFFT比对（输入已是比对结果） |
| `--mafft-path` | `mafft` | MAFFT可执行文件路径 |
| `--mafft-params` | `--auto` | MAFFT参数 |
| `--threads` | `4` | MAFFT线程数 |
| `--keep-alignment` | `True` | 保留比对结果文件（默认）|
| `--no-keep-alignment` | - | 不保留比对结果文件 |

### 颜色方案 | Color Schemes

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--color-scheme` | `Zappo` | 颜色方案 |

**可用的颜色方案**:
- `Clustal` - ClustalX风格颜色
- `Zappo` - Zappo颜色方案（氨基酸默认）
- `Taylor` - Taylor颜色方案
- `Flower` - Flower花色方案
- `Blossom` - Blossom花色方案
- `Sunset` - 日落色方案
- `Ocean` - 海洋色方案
- `Hydrophobicity` - 疏水性颜色
- `HelixPropensity` - 螺旋倾向性
- `StrandPropensity` - 折叠倾向性
- `TurnPropensity` - 转角倾向性
- `BuriedIndex` - 埋藏指数
- `Nucleotide` - 核苷酸颜色（核苷酸默认）
- `Purine/Pyrimidine` - 嘌呤/嘧啶
- `Identity` - 一致性颜色
- `None` - 无颜色

### 区域参数 | Region Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--start` | `1` | 起始位置(1-based) |
| `--end` | `MSA长度` | 结束位置(1-based) |

### 显示参数 | Display Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--wrap-length` | `60` | 换行长度 |
| `--wrap-space-size` | `3.0` | 换行间距 |
| `--label-type` | `id` | 标签类型 (id/description) |
| `--show-grid` | `False` | 显示网格 |
| `--show-count` | `False` | 显示字符统计 |
| `--show-consensus` | `False` | 显示consensus序列 |
| `--consensus-color` | `#1f77b4` | Consensus颜色 |
| `--consensus-size` | `2.0` | Consensus大小 |

### 排序选项 | Sorting Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--sort` | `False` | 按NJ树排序 |

### 输出参数 | Output Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--dpi` | `300` | 图像DPI |

## 使用示例 | Usage Examples

### 示例1：蛋白质序列自动比对可视化 | Example 1: Protein Auto-Align & Visualize

```bash
# 自动比对蛋白质序列并可视化
biopytools msaviz \
    -i protein_sequences.fa \
    -o protein_viz.png \
    --color-scheme Taylor \
    --show-consensus \
    --show-count
```

### 示例2：核苷酸序列比对 | Example 2: Nucleotide Alignment

```bash
# 核苷酸序列比对（自动使用Nucleotide颜色方案）
biopytools msaviz \
    -i nucleotide_sequences.fa \
    -o nucleotide_viz.png \
    --wrap-length 100 \
    --show-grid \
    --show-consensus \
    --threads 8
```

### 示例3：不保留比对结果 | Example 3: Do Not Keep Alignment

```bash
# 比对但不保存结果文件（使用临时文件）
biopytools msaviz \
    -i sequences.fa \
    -o output.png \
    --no-keep-alignment
# 输出: 仅 output.png
```

### 示例4：自定义MAFFT参数 | Example 4: Custom MAFFT Parameters

```bash
# 使用精确的MAFFT参数
biopytools msaviz \
    -i sequences.fa \
    -o output.png \
    --mafft-params "--localpair --maxiterate 1000" \
    --threads 16
```

### 示例5：已有比对结果 | Example 5: Already Aligned Input

```bash
# 输入已经是比对结果，跳过比对
biopytools msaviz \
    -i aligned_sequences.fa \
    -o output.png \
    --skip-align \
    --color-scheme Clustal \
    --show-consensus
```

### 示例6：高质量输出用于发表 | Example 6: High-Quality Output

```bash
# 高分辨率PDF输出（用于发表）
biopytools msaviz \
    -i sequences.fa \
    -o publication_viz.pdf \
    --color-scheme Taylor \
    --wrap-length 80 \
    --show-consensus \
    --show-count \
    --dpi 600
```

### 示例7：区域可视化 | Example 7: Region Visualization

```bash
# 只显示特定区域
biopytools msaviz \
    -i aligned.fa \
    -o region_viz.png \
    --skip-align \
    --start 50 \
    --end 150 \
    --color-scheme Clustal \
    --dpi 600
```

## 工作流程 | Workflow

### 默认流程（自动比对）| Default Workflow (Auto-align)

```
输入序列文件 (FASTA)
    ↓
[可选] 检测序列类型 (核苷酸/蛋白质)
    ↓
运行 MAFFT 比对
    ↓
生成比对结果 (默认保存到输出目录)
    ↓
pyMSAviz 可视化
    ↓
输出图像文件 (PNG/JPG/SVG/PDF) + 比对结果文件 (.aligned.fa)
```

### 跳过比对流程 | Skip Alignment Workflow

```
输入已比对文件 (FASTA/PHYLIP/CLUSTAL)
    ↓
pyMSAviz 可视化
    ↓
输出图像文件
```

## 输入文件格式 | Input File Formats

### FASTA格式（序列文件）

```fasta
>seq1
ATCGATCGATCG
>seq2
ATCGATCGATCG
>seq3
ATCGATCGATCG
```

### 已比对FASTA格式

```fasta
>seq1
ATCGATCG---ATCG
>seq2
ATCG---CGATATCG
>seq3
ATCGATCGATCG---
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **MAFFT** (版本 7+)
  - 下载地址: https://mafft.cbrc.jp/alignment/software/
- **Python** (版本 3.9+)
- **Python包**:
  - `matplotlib` >= 3.5.2
  - `biopython` >= 1.79
  - `click` - 命令行界面

### 安装依赖 | Installing Dependencies

```bash
# 安装MAFFT
# Ubuntu/Debian
sudo apt-get install mafft

# CentOS/RHEL
sudo yum install mafft

# macOS
brew install mafft

# 或从源码安装
wget https://mafft.cbrc.jp/alignment/software/mafft-7.520-without-extensions.tar.gz
tar -xzf mfft-7.520-without-extensions.tar.gz
cd mafft-7.520-without-extensions/
make clean
make
sudo make install

# 安装Python包
pip install matplotlib biopython click
```

## 注意事项 | Important Notes

1. **默认行为**: 默认会自动运行MAFFT比对并保留比对结果文件
2. **文件格式**: 确保输入文件是标准FASTA格式
3. **内存使用**: 大型序列文件可能需要较多内存
4. **DPI设置**: 高DPI值会显著增加输出文件大小和处理时间
5. **颜色方案**: 蛋白质序列默认使用Zappo，核苷酸序列默认使用Nucleotide
6. **MAFFT路径**: 如果MAFFT不在PATH中，使用`--mafft-path`指定路径
7. **比对文件**: 默认保存比对结果为`输出文件名.aligned.fa`，使用`--no-keep-alignment`不保留

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "MAFFT not found" 错误**
```bash
# 安装MAFFT
sudo apt-get install mafft  # Ubuntu/Debian

# 或指定MAFFT路径
biopytools msaviz -i seq.fa -o out.png --mafft-path /path/to/mafft
```

**Q: 想要使用已比对的文件**
```bash
# 使用--skip-align跳过比对步骤
biopytools msaviz -i aligned.fa -o output.png --skip-align
```

**Q: 如何加快比对速度**
```bash
# 增加线程数
biopytools msaviz -i seq.fa -o out.png --threads 16

# 使用快速MAFFT参数
biopytools msaviz -i seq.fa -o out.png --mafft-params "--retree 2 --maxiterate 0"
```

**Q: 输出图像太大**
```bash
# 降低DPI或减少wrap-length
biopytools msaviz -i seq.fa -o out.png --dpi 150 --wrap-length 40
```

**Q: 如何获得更精确的比对**
```bash
# 使用精确的MAFFT参数
biopytools msaviz -i seq.fa -o out.png --mafft-params "--localpair --maxiterate 1000"
```

## 相关资源 | Related Resources

- [pyMSAviz官方文档](https://moshi4.github.io/pyMSAviz/)
- [MAFFT官方文档](https://mafft.cbrc.jp/alignment/software/)
- [BioPython AlignIO文档](https://biopython.org/wiki/AlignIO)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

pyMSAviz同样采用MIT许可证，感谢原作者moshi4的贡献。

## 引用信息 | Citation

如果在学术研究中使用此工具，请分别引用MAFFT和pyMSAviz：

```
MAFFT:
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780.

pyMSAviz:
https://github.com/moshi4/pyMSAviz
```
