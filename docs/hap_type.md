# hap-type | 单倍型提取与可视化

**从 VCF 指定区间提取单倍型，输出兼容 geneHapR 的单倍型结果表 | Extract haplotypes from a VCF region and export geneHapR-compatible tables**

## 功能概述 | Overview

本模块将给定基因组区间（或 BED 文件中的一批区间）内的多个 SNP/INDEL 位点组合成"单倍型块"（haplotype block），为每个样本根据其 GT 生成单倍型字符串，并对相同单倍型分组、编号、统计频率。输出格式与 [geneHapR](https://github.com/ZhangRenL/geneHapR) 的 `hapResult` / `hapSummary` / `sampleHap` 表完全兼容，方便直接导入 geneHapR 做下游可视化（如单倍型网络图、频率分布图）。

流程对每个 BED 行（或单区间）独立运行，可选去除杂合位点（无法确定单一单倍型）与缺失位点。单倍型 ID 默认形如 `H001`、`H002`，前缀与位数可调。

## 快速开始 | Quick Start

```bash
# 单个区间
biopytools hap-type -i sample.vcf -r chr1:1000000-1005000 -o result

# 批量处理（BED 文件）
biopytools hap-type -i sample.vcf -r regions.bed -o result
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --vcf` | VCF 变异文件 |
| `-r, --region` | 基因组区间 `chr:start-end` 或 BED 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | 自动生成 | 输出文件前缀（省略则按区间自动生成） |
| `--hetero-remove` | `False` | 去除杂合位点 |
| `--na-drop` | `True` | 去除缺失位点 |
| `--hap-prefix` | `H` | 单倍型 ID 前缀 |
| `--pad` | `3` | 单倍型 ID 位数（如 3 → `H001`） |

（运行 `biopytools hap-type -h` 查看完整参数列表）

## 输出 | Output

单区间模式（输出目录由 `-o` 决定）：

```
{output_dir}/
├── {chrom}_{start}_{end}.hapResult.txt       # 单倍型矩阵（行=单倍型，列=位点+样本）
├── {chrom}_{start}_{end}.hapResult.xlsx      # Excel 版本
├── {chrom}_{start}_{end}.hapSummary.txt      # 单倍型汇总（ID、频率、位点序列）
├── {chrom}_{start}_{end}.hapSummary.xlsx
├── {chrom}_{start}_{end}.sampleHap.txt       # 样本-单倍型映射
├── {chrom}_{start}_{end}.sampleHap.xlsx
└── hap_type.log                              # 日志
```

BED 模式则对每一行生成对应的一组文件。

## 依赖 | Dependencies

- Python：标准库（pandas 用于 Excel 输出时需要 `openpyxl`）
- 可选：[geneHapR](https://github.com/ZhangRenL/geneHapR)（用于下游 R 端可视化）

## 引用 | Citation

- Zhang R, Chen Y, Yu W, et al. geneHapR: an R package for gene haplotype statistics and visualization. *BMC Bioinformatics*, 2023, 24: 244. doi:10.1186/s12859-023-05513-1

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [geneHapR GitHub](https://github.com/ZhangRenL/geneHapR)
