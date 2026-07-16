# 基因信息+序列合并表 | Gene Info + Sequence Table

**从基因组和GFF3提取基因信息表, 同时输出基因DNA + CDS + 蛋白序列 | Extract gene info table from genome + GFF3, with gene DNA + CDS + protein sequences**

## 功能概述 | Overview

gene_table 模块生成一份基因信息合并表, 每行一个基因, 同时导出对应的基因 DNA、CDS、蛋白序列。CDS/蛋白通过 gffread 从 GFF3 + 基因组提取(路径自动检测, 也可显式指定)。

支持每基因仅保留最长转录本(`--longest-only`), 或保留全部转录本(默认)。

## 快速开始 | Quick Start

```bash
# 基础用法(全部转录本)
biopytools gene-table -g genome.fa -f annotation.gff3 -o gene_table.tsv

# 每基因仅最长转录本, 指定前缀(作为 Sample 列)
biopytools gene-table -g genome.fa -f anno.gff3 -o out.tsv --longest-only --prefix SAMPLE1
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 基因组 FASTA |
| `-f, --gff` | GFF3 注释(支持 `.gz`) |
| `-o, --output` | 输出表路径(或目录) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--prefix` | GFF 文件名 | 输出文件前缀 + Sample 列 |
| `--longest-only` | `False` | 每基因仅保留最长转录本 |
| `--transcript-types` | `mRNA transcript` | 视为转录本的 feature 类型 |
| `--gene-type` | `gene` | 基因 feature 类型 |
| `--min-length` | `0` | 基因 DNA 最小长度过滤(0=不过滤) |
| `--gffread` | 自动检测 | gffread 路径 |
| `-v, --verbose` | `False` | 详细日志 |

(运行 `biopytools gene-table -h` 查看完整参数列表)

## 输出 | Output

- 信息合并表 TSV(基因坐标 + 序列文件路径)
- `{prefix}.gene.fa`: 基因 DNA 序列
- `{prefix}.cds.fa`: CDS 序列
- `{prefix}.pep.fa`: 蛋白序列

## 依赖 | Dependencies

- **gffread**: 从 GFF3 + 基因组提取 CDS/蛋白 (随 StringTie 发布)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
