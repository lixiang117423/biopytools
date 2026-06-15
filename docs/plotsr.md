# PlotSR 模块使用文档 | PlotSR Module Usage Guide

## 功能概述 | Overview

PlotSR模块封装了完整的minimap2 + SyRI + PlotSR流程，用于多基因组共线性可视化。

## 快速开始 | Quick Start

```bash
# 基本用法
biopytools plotsr -i genome1.fa -i genome2.fa -i genome3.fa -o output/

# 使用map文件（推荐）
biopytools plotsr -i genomes.map -o output/

# 只分析特定染色体
biopytools plotsr -i genomes.map -o output -c 1
biopytools plotsr -i genomes.map -o output -c "1,2,3"
```

## 输入模式 | Input Modes

### 1. 多次使用`-i`参数（Multiple genomes）
```bash
biopytools plotsr \
  -i genome1.fa \
  -i genome2.fa \
  -i genome3.fa \
  -o output/
```

### 2. 使用文件夹（Auto-discovery）
```bash
# 自动发现文件夹中的所有FASTA文件
biopytools plotsr -i ./genome_folder/ -o output/
```

### 3. 使用map文件（Recommended）
map文件格式（tab分隔）：
```
name1	path/to/genome1.fa
name2	path/to/genome2.fa
name3	path/to/genome3.fa
```

示例：
```bash
biopytools plotsr -i genomes.map -o output/
```

**优先级**：map文件名称 > `-n`参数 > 自动提取

## 参数说明 | Parameters

### 基本参数 | Basic Parameters

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i` | 输入基因组（可多次使用）| 必需 |
| `-o` | 输出目录 | 必需 |
| `-n` | 基因组名称（逗号分隔）| 自动提取 |
| `-t` | 线程数 | 12 |

### 比对参数 | Alignment Parameters

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--minimap2-preset` | minimap2预设 | asm5 |

### 结构变异参数 | Structural Variant Parameters

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-s` | 最小结构变异大小 | 10000 |

### 可视化参数 | Visualization Parameters

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-f` | 字体大小 | 6 |
| `-d` | 图片DPI | 300 |
| `--space-ratio` | 同源染色体间距(0.1-0.75) | 0.7 |
| `--output-format` | 输出格式(pdf/png/svg) | pdf |
| `-v` | 垂直排列染色体 | False |
| `--itx` | 染色体间交互模式 | False |

### 过滤参数 | Filtering Parameters

| 参数 | 说明 |
|------|------|
| `--nosyn` | 不绘制同源区域 |
| `--noinv` | 不绘制倒位 |
| `--notr` | 不绘制易位 |
| `--nodup` | 不绘制重复 |

### 染色体过滤 | Chromosome Filtering

```bash
# 单个染色体（数字索引）
biopytools plotsr -i genomes.map -o output -c 1

# 多个染色体（逗号分隔）
biopytools plotsr -i genomes.map -o output -c "1,2,3"

# 使用染色体名称
biopytools plotsr -i genomes.map -o output -c hifiasm_Chr01

# 混合使用
biopytools plotsr -i genomes.map -o output -c "1,2" -c hifiasm_Chr05
```

**注意**：染色体索引从1开始，按照第一个基因组的染色体顺序。

### 流程控制 | Pipeline Control

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--skip-existing` | 跳过已完成的步骤 | True |
| `--force-run` | 强制重新运行所有步骤 | False |

## 输出文件 | Output Files

```
output/
├── alignment/                    # 比对结果
│   ├── genome1_vs_genome2.bam
│   └── genome2_vs_genome3.bam
├── syri/                         # SyRI结构注释
│   ├── genome1_vs_genome2syri.filtered.out
│   └── genome2_vs_genome3syri.filtered.out
├── plotsr/                       # PlotSR配置文件
│   ├── genomes.txt
│   └── *.chrlen
├── plot.pdf                      # 最终可视化结果
└── plotsr.log                    # 运行日志
```

## 使用示例 | Examples

### 示例1：三个基因组全分析
```bash
biopytools plotsr \
  -i genome1.fa \
  -i genome2.fa \
  -i genome3.fa \
  -o output_full \
  -t 24
```

### 示例2：使用map文件并指定名称
```bash
# genomes.map:
# Col-0	/data/Col-0.fa
# Ler	/data/Ler.fa
# Cvi	/data/Cvi.fa

biopytools plotsr -i genomes.map -o output_col/
```

### 示例3：只分析前3条染色体
```bash
biopytools plotsr -i genomes.map -o output_chr3 -c "1,2,3"
```

### 示例4：高质量PNG输出
```bash
biopytools plotsr \
  -i genomes.map \
  -o output_png \
  --output-format png \
  -d 600 \
  -f 8
```

### 示例5：只显示结构变异，不显示同源区域
```bash
biopytools plotsr \
  -i genomes.map \
  -o output_sv \
  --nosyn
```

## 故障排除 | Troubleshooting

### 问题1：PlotSR报错染色体长度不匹配

**原因**：这是正常的，因为PlotSR会验证SyRI文件中所有染色体的长度。

**解决**：这个错误已经在v1.0.1中修复，确保使用最新版本。

### 问题2：分析时间过长

**原因**：全基因组比对和SyRI分析需要较长时间。

**解决**：
- 使用`-c`参数只分析特定染色体
- 增加`-t`参数提高线程数
- 使用`--skip-existing`跳过已完成步骤

### 问题3：输出图片不清晰

**解决**：
- 增加DPI：`-d 600`
- 增加字体：`-f 8`
- 使用PDF格式：`--output-format pdf`

## 版本历史 | Version History

- **v1.0.0** (2026-02-25)
  - 初始版本
  - 支持minimap2 + SyRI + PlotSR完整流程
  - 支持map文件输入
  - 支持自动跳过已完成步骤
  - 支持染色体过滤

- **v1.0.1** (2026-02-25)
  - 修复SyRI文件顺序问题
  - 修复PlotSR染色体长度验证问题

## 相关链接 | Related Links

- PlotSR GitHub: https://github.com/schneebergerlab/plotsr
- SyRI文档: https://github.com/schneebergerlab/syri
- minimap2文档: https://github.com/lh3/minimap2
