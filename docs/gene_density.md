# 基因密度计算模块

**按固定窗口统计每条染色体各区间的基因数量与密度（基因/Mb） | Count genes and gene density (genes/Mb) per fixed-size window along each chromosome**

## 功能概述 | Overview

gene_density 模块读取 GFF3 注释文件，按用户指定窗口大小（默认 100 kb）沿每条染色体滑动统计基因（或其他 feature）数量与密度，输出结构化 TSV 供下游绘图/分析，并可选生成每条染色体的密度分布图。纯 Python 实现（无外部生信工具依赖），计算快、可在登录节点运行。支持指定 `.fai`/FASTA 作为染色体长度来源以提升末尾窗口精度。支持断点续传。

## 快速开始 | Quick Start

```bash
# 基本计算（默认100kb窗口，生成TSV+图）
biopytools gene-density -i annotation.gff3 -o output_dir/

# 指定窗口大小与前缀
biopytools gene-density -i annotation.gff3 -o output_dir/ -w 500000 --prefix sample1

# 用.fai提升末尾窗口精度，且不画图
biopytools gene-density -i annotation.gff3 -o output_dir/ -g genome.fa.fai --no-plot
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --gff` | GFF3 注释文件 |
| `-o, --output-dir` | 输出目录（默认 `./gene_density_output`） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-w, --window-size` | `100000` | 窗口大小（bp） |
| `--feature-type` | `gene` | 统计的 GFF feature 类型（GFF 第 3 列，如 `gene`/`CDS`/`mRNA`） |
| `-g, --genome` | 空 | 染色体长度来源（`.fai` 或 FASTA，可选，提升末尾窗口精度） |
| `--prefix` | GFF 文件名 stem | 输出文件前缀 |
| `--no-plot` | `False` | 不绘制密度图 |

（运行 `biopytools gene-density -h` 查看完整参数列表）

## 输出 | Output

```
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml        # 运行环境与参数记录
├── 01_density/
│   └── {prefix}.gene_density.tsv    # 基因密度TSV ⭐
├── 02_plot/
│   └── {prefix}.gene_density.png    # 每条染色体密度分布图
└── 99_logs/
    └── gene_density.log             # 运行日志
```

`{prefix}.gene_density.tsv` 主要列：

| 列名 | 描述 |
|------|------|
| `chrom` | 染色体/序列名 |
| `start` / `end` | 窗口起止（bp） |
| `gene_count` | 窗口内基因（或指定 feature）数量 |
| `density_genes_per_Mb` | 密度（基因/Mb） |

> 末尾不足一个窗口的区间：未提供 `-g` 长度来源时按实际基因覆盖区间统计；提供后按真实染色体长度切窗。

## 依赖 | Dependencies

- Python ≥ 3.10（纯标准库解析 GFF）
- `matplotlib`（绘图；`--no-plot` 时不需要）
- 无需 conda 环境/外部生信工具

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
