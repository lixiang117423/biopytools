# GFF 整理与 ID 重命名 | GFF Converter (renamegff)

**整理 GFF/GFF3 注释文件并为基因/转录本生成统一规范的 ID | Clean up GFF/GFF3 annotations and assign uniform, species-prefixed IDs to genes and transcripts**

> 注意 | Note：CLI 命令名是 `renamegff`（不是 `gffconverter`）。模块目录名为 `gffconverter`，但实际调用入口为 `biopytools renamegff`。

## 功能概述 | Overview

`renamegff` 用于在基因组注释完成后对 GFF/GFF3 文件做统一整理：解析所有 feature 行，按 `gene → mRNA → CDS/exon` 的层级关系重建 ID，并为所有基因生成形如 `{species_prefix}g{number}` 的连续编号（如 `Ovg010`、`Ovg020`），mRNA 与 CDS 则自动派生为 `{gene_id}.t1`、`{mrna_id}.cds` 等规范命名。

工具可设置起始编号和步长（默认从 10 开始，步长 10），便于多物种/多注释合并时预留区间。支持 `.gz` 压缩输入、多线程加速、保留中间文件、转换样本预览等。典型场景包括：多个结构预测软件结果合并前的 ID 统一、公开发布前的 GFF3 规范化、为不同物种批次分配独立 ID 区间（如物种 A 从 1000 开始、物种 B 从 2000 开始）。

## 快速开始 | Quick Start

```bash
# 基本用法：物种 OV53，前缀 Ov
biopytools renamegff -i input.gff -o output.gff -s OV53 -p Ov

# 自定义起始编号和步长，并预览 10 个转换样本
biopytools renamegff -i annotation.gff3 -o standardized.gff \
    -s SpeciesA -p SpA --start-num 100 --step 5 --show-sample 10
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 GFF/GFF3 文件（支持 `.gz`）|
| `-o, --output` | 输出 GFF 文件路径 |
| `-s, --species-name` | 物种名称（如 `OV53`）|
| `-p, --species-prefix` | 物种前缀，用于生成 gene ID（如 `Ov`）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--start-num` | `10` | 基因编号起始值 |
| `--step` | `10` | 基因编号步长 |
| `-t, --threads` | `12` | 线程数 |
| `-v, --verbose` | 关 | 详细输出模式 |
| `--keep-intermediate` | 关 | 保留中间文件（parsed_features、id_mappings 等）|
| `--show-sample N` | 无 | 显示 N 个转换样本后退出（用于预览）|

（运行 `biopytools renamegff -h` 查看完整参数列表）

### ID 生成规则 | ID Naming Convention

- 基因 ID：`{prefix}g{number}`，例如前缀 `Ov`、起始 10、步长 10 → `Ovg010`, `Ovg020`, `Ovg030` ...
- mRNA ID：`{gene_id}.t1`, `{gene_id}.t2`（按转录本顺序）
- CDS ID：`{mrna_id}.cds`
- 蛋白 ID：`{gene_id}.p1`, `{gene_id}.p2`

## 输出 | Output

```
./
├── output.gff               # 整理后的标准 GFF3
├── output.gff.log           # 转换日志
├── conversion_stats.txt     # 统计信息（verbose 模式）
└──（--keep-intermediate 时）
    ├── parsed_features.tmp  # 解析后的特征
    ├── id_mappings.txt      # 新旧 ID 映射表
    └── validation_report.txt # 验证报告
```

## 依赖 | Dependencies

- Python 3.7+（纯 Python 实现，无强制第三方依赖）

## 引用 | Citation

- GFF3 规范：Sequence Ontology Project, Generic Feature Format Version 3
- Wilkinson, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data* 3, 160018 (2016).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
