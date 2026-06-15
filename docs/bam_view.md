# BAM比对可视化分析模块

**专业的BAM文件比对可视化工具 | Professional BAM Alignment Visualization Tool**

## 功能概述 | Overview

BAM比对可视化分析模块是一个强大的BAM文件可视化工具，基于alignoth软件构建，提供从BAM文件生成交互式比对可视化的完整功能。支持多种输出格式、灵活的高亮选项和详细的配置参数，适用于各种基因组比对分析研究。

## 主要特性 | Key Features

- **多种输出格式**: 支持HTML、Vega-Lite JSON、SVG、PDF等多种输出格式
- **交互式可视化**: HTML格式支持交互式探索，可缩放、高亮显示
- **灵活的高亮选项**: 支持VCF变异位点、BED区域、自定义区间的高亮显示
- **详细的配置参数**: 可自定义reads显示深度、图像宽度、错配显示等
- **覆盖度统计**: 自动计算并显示覆盖度信息
- **错配可视化**: 显示序列错配和变异信息
- **高质量输出**: 基于Vega-Lite的专业级可视化效果

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 生成HTML格式的比对可视化
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr1:1000-2000 \
    -o output_dir

# 指定输出格式为JSON
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr1:5000-6000 \
    -f json \
    -o json_output
```

### 高级用法 | Advanced Usage

```bash
# 使用VCF文件高亮变异位点
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr2:10000-11000 \
    -v variants.vcf \
    -o vcf_highlight

# 使用BED文件高亮区域
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr3:20000-21000 \
    --bed regions.bed \
    -o bed_highlight

# 自定义高亮区间和辅助标签
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr4:30000-31000 \
    -H variant1:30500-30510 \
    -H variant2:30700-30705 \
    -x MD \
    -x NM \
    -o custom_highlight
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-b, --bam` | BAM文件路径 | `-b alignments.bam` |
| `-r, --reference` | 参考序列FASTA文件路径 | `-r reference.fa` |
| `-g, --region` | 可视化区域(格式: chr:start-end) | `-g chr1:1000-2000` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--alignoth-path` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/alignoth/bin/alignoth` | alignoth软件路径 |

### 输出配置 | Output Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./bam_view_output` | 输出目录路径 |
| `-f, --output-format` | `html` | 输出格式 (html/json/svg/pdf) |

### 可视化参数 | Visualization Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-d, --max-read-depth` | `500` | 最大reads显示深度 |
| `-w, --max-width` | `1024` | 最大宽度 |
| `--mismatch-display-min-percent` | `1.0` | 显示错配的最小百分比 |

### 高亮选项 | Highlight Options

| 参数 | 描述 |
|------|------|
| `-v, --vcf` | VCF文件路径(高亮变异位点) |
| `--bed` | BED文件路径(高亮区域) |
| `-H, --highlight` | 高亮区间(可多次使用, 格式: name:start-end) |

### 其他选项 | Other Options

| 参数 | 描述 |
|------|------|
| `-x, --aux-tag` | 辅助标签(可多次使用) |
| `--no-embed-js` | 不嵌入JavaScript(仅HTML格式，减小文件大小) |
| `--plot-all` | 绘制所有reads(不推荐大文件使用) |

## 输入文件格式 | Input File Formats

### BAM文件 | BAM File

标准的BAM格式比对文件，需要建立索引：

```bash
# 如果没有索引，需要先建立索引
samtools index alignments.bam
```

**文件要求**:
- 标准BAM格式（已排序）
- 必须有.bai索引文件
- 包含比对信息和序列数据

### 参考序列文件 | Reference Sequence File

标准FASTA格式的参考基因组序列：

```fasta
>chr1
ATCGATCGATCGATCGATCGATCGATCG...
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

### VCF文件 | VCF File (可选)

标准VCF格式的变异文件，用于高亮变异位点：

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000000>
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    1250 .   A    G    45.2  PASS    DP=30
chr1    1340 .   T    C    62.8  PASS    DP=25
```

**文件要求**:
- 标准VCF格式（版本4.2+）
- 必须有.csi或.tbi索引文件
- 使用前需用bcftools index或tabix建立索引

### BED文件 | BED File (可选)

标准BED格式的区域文件：

```bed
chr1    1000    1500    region1
chr1    2000    2500    region2
```

## 使用示例 | Usage Examples

### 示例1：基本HTML可视化 | Example 1: Basic HTML Visualization

```bash
# 生成HTML格式的交互式可视化
biopytools bam-view \
    -b sample1.bam \
    -r genome.fa \
    -g chr1:50000-51000 \
    -o html_output
```

### 示例2：变异位点高亮 | Example 2: Variant Highlighting

```bash
# 使用VCF文件高亮所有变异位点
biopytools bam-view \
    -b tumor_sample.bam \
    -r human_genome.fa \
    -g chr17:7670000-7680000 \
    -v somatic_variants.vcf \
    -o variant_view
```

### 示例3：自定义区间高亮 | Example 3: Custom Interval Highlighting

```bash
# 高亮特定的多个区间
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr5:1000000-1001000 \
    -H SNP1:1000500-1000501 \
    -H SNP2:1000800-1000801 \
    -H deletion:1000900-1000950 \
    -o custom_highlights
```

### 示例4：基因区域详细分析 | Example 4: Gene Region Detailed Analysis

```bash
# 查看特定基因区域的比对情况
biopytools bam-view \
    -b rnaseq_alignments.bam \
    -r transcriptome.fa \
    -g geneX:2000-3000 \
    -d 1000 \
    -w 2048 \
    -x MD \
    -x NM \
    -o gene_detail
```

### 示例5：生成Vega-Lite JSON | Example 5: Generate Vega-Lite JSON

```bash
# 输出Vega-Lite JSON格式，用于自定义可视化
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr1:10000-11000 \
    -f json \
    -o json_output
```

### 示例6：生成PDF报告 | Example 6: Generate PDF Report

```bash
# 生成PDF格式的静态报告（需要vega-cli）
biopytools bam-view \
    -b alignments.bam \
    -r reference.fa \
    -g chr2:20000-21000 \
    -f pdf \
    -o pdf_report
```

## 输出结果 | Output Results

### HTML输出 | HTML Output

交互式HTML文件，支持：
- 缩放和平移
- 点击reads查看详细信息
- 高亮显示变异区域
- 覆盖度统计图表

**文件特点**:
- 完整自包含（嵌入JavaScript）
- 可直接在浏览器中打开
- 支持现代浏览器的所有功能

### JSON输出 | JSON Output

Vega-Lite规范JSON文件，可用于：
- 自定义可视化
- 集成到网页中
- 转换为其他格式（SVG/PDF）

### SVG输出 | SVG Output

矢量图形文件：
- 无损缩放
- 适合出版和打印
- 可用矢量图软件编辑

### PDF输出 | PDF Output

便携式文档文件：
- 适合报告和论文
- 固定布局
- 高质量输出

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **alignoth** (版本 1.4.0 或更新)
  - 安装方式: `conda install -c bioconda alignoth`
  - 或: `pixi global install alignoth`
  - 或: `cargo install alignoth`

- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 可选依赖 | Optional Dependencies

用于SVG/PDF输出：
- **Node.js** 和 **vega-cli**
  ```bash
  npm install -g vega-cli vega-lite-cli
  ```

### 安装依赖软件 | Installing Dependencies

```bash
# 使用conda安装alignoth（推荐）
conda install -c conda-forge -c bioconda alignoth

# 或使用pixi安装
pixi global install alignoth

# 安装Python包
pip install click

# 安装vega-cli（可选，用于SVG/PDF输出）
npm install -g vega-cli vega-lite-cli
```

## 注意事项 | Important Notes

1. **BAM索引**: BAM文件必须有索引（.bai），使用`samtools index`创建
2. **VCF索引**: 如果使用VCF文件，必须有索引（.csi或.tbi），使用`bcftools index`创建
3. **区域格式**: 区域格式为`chr:start-end`，坐标为1-based且完全包含
4. **内存使用**: 大深度或大区域可能需要较多内存
5. **性能**: 建议单个区域不超过10kb，reads深度不超过1000

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "BAM索引文件不存在" 错误**
```bash
# 创建BAM索引
samtools index alignments.bam
```

**Q: "alignoth not found" 错误**
```bash
# 检查alignoth安装
which alignoth

# 如未安装，使用conda安装
conda install -c bioconda alignoth
```

**Q: "vl2vg not found" 错误**
```bash
# 安装vega-cli（用于SVG/PDF输出）
npm install -g vega-cli vega-lite-cli
```

**Q: 可视化显示不完整**
```bash
# 增加max-read-depth参数
biopytools bam-view ... -d 1000

# 或增加max-width参数
biopytools bam-view ... -w 2048
```

**Q: HTML文件太大**
```bash
# 使用--no-embed-js参数
biopytools bam-view ... --no-embed-js

# 注意：需要网络连接来加载JavaScript库
```

## 与其他工具比较 | Comparison with Other Tools

| 特性 | alignoth (bam-view) | IGV | samtools tview |
|------|---------------------|-----|----------------|
| 交互式HTML | 支持 | 支持 | 不支持 |
| 批量处理 | 支持 | 不支持 | 不支持 |
| 命令行 | 原生 | GUI | 原生 |
| 高亮变异 | 支持 | 支持 | 不支持 |
| 覆盖度图 | 自动 | 手动 | 不显示 |
| 导出格式 | 多种 | 截图 | 文本 |

## 相关资源 | Related Resources

- [alignoth官方文档](https://github.com/koesterlab/alignoth)
- [alignoth在线示例](https://alignoth.github.io/preview.html)
- [Vega-Lite文档](https://vega.github.io/vega-lite/)
- [SAM格式规范](https://samtools.github.io/hts-specs/)

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用alignoth相关文献：

```
Wiegand, F., & Köster, J. (2024).
alignoth: publication-ready alignment plots from BAM files.
Bioinformatics, 40(1), btaf663.
```
