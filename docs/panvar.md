# 泛基因组变异分析工具 | Pan-genome Variant Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2026-04-21

## 概述 | Overview

panvar封装了从Cactus泛基因组图提取结构变异并进行统计分析的完整流程，包含vg deconstruct变异提取、变异分类统计、泛基因组增长曲线绘制和Heaps' Law gamma值拟合。

The panvar tool wraps the complete pipeline from Cactus pangenome graph variant extraction to statistical analysis, including vg deconstruct, variant classification, growth curve plotting, and Heaps' Law gamma fitting.

自动识别输入类型: .gbz文件先执行deconstruct再统计，.vcf文件直接统计。
Auto-detects input type: .gbz runs deconstruct first then summary, .vcf runs summary directly.

## 功能特点 | Features

- **vg deconstruct**: 从GBZ泛基因组图提取结构变异VCF (输入为.gbz时自动执行)
- **变异分类统计**: SNP/InDel(2-50bp)/SV(>50bp)分类与等位基因频率统计
- **增长曲线+gamma**: R/ggplot2绘制泛基因组增长曲线，拟合Heaps' Law: y = a * x^gamma
- **断点续传**: 自动跳过已完成的步骤

## 分析流程 | Analysis Pipeline

1. **vg deconstruct** (仅GBZ输入): 从Cactus GBZ图提取变异VCF
2. **变异分类统计**: SNP/InDel/SV分类 + 等位基因频率统计 (Python/pysam)
3. **增长曲线+gamma**: 随机置换 + Heaps' Law拟合 + ggplot2绘图 (R)

## 使用方法 | Usage

### 依赖 | Dependencies

- vg (conda环境: vg_v.1.7.0)
- R + ggplot2
- Python: pysam

### GBZ输入 (完整流程)

```bash
biopytools panvar -i cactus_output.gbz -P T2T -o output/
```

### VCF输入 (跳过deconstruct)

```bash
biopytools panvar -i cactus_SV.vcf -o output/
```

## 命令行参数 | Command Line Options

| 参数 | 说明 |
|------|------|
| `-i, --input` | GBZ图文件或VCF文件 (自动识别) |
| `-o, --output` | 输出目录 |
| `-P, --ref-path` | 参考路径前缀 (GBZ输入时必需, 如T2T) |
| `-t, --threads` | 线程数 (default: 12) |
| `--ref-size` | 参考基因组大小Mb (0=自动推断) |
| `--permutations` | 增长曲线随机置换次数 (default: 100) |
| `--vg-env` | vg conda环境名 (default: vg_v.1.7.0) |
| `--r-path` | Rscript路径 |

## 输出文件 | Output Files

```
output/
├── 01_deconstruct/                       # (仅GBZ输入)
│   └── cactus_output.vcf
├── 02_summarize/
│   ├── pangenome_variants_summary.tsv     # 变异分类统计表
│   ├── pangenome_growth_curve.png         # 增长曲线 (含gamma值和拟合线)
│   └── plot_growth.R                     # 生成的R脚本 (可复现)
└── 99_logs/
    └── panvar.log
```

### gamma值说明

增长曲线图副标题显示拟合结果: `gamma = X.XXXX | R-squared = X.XXXX`

gamma值反映了泛基因组的开放程度:
- gamma接近0: 泛基因组趋于封闭 (核心基因组为主)
- gamma接近1: 泛基因组高度开放 (大量特有变异)
