# vcf-filter | VCF 筛选工具

**按区间、质量、深度、样本、ID、变异类型等多维度筛选 VCF，可选 PLINK 格式转换 | Filter a VCF by region, quality, depth, sample, ID, variant type and optionally convert to PLINK format**

## 功能概述 | Overview

本模块是一个综合性的 VCF 筛选工具，整合了常见的筛选维度：染色体与坐标区间提取、MAF / 缺失率 / 质量 / 深度阈值过滤、保留或移除指定样本或变异 ID、只保留双等位点、移除 INDEL 等。所有筛选支持任意组合，匹配任一未通过的位点都会被剔除。

为提高大文件处理速度，默认跳过输入验证（可用 `--force-validation` 强制开启）。需要做 GWAS、PCA 等分析时，可启用 `--convert-format` 通过 PLINK 转换为 `.bed/.bim/.fam` 二进制格式。

## 快速开始 | Quick Start

```bash
# 按区间提取
biopytools vcf-filter -i input.vcf -c chr1 -s 1000 -e 2000 -o chr1_region.vcf

# 综合过滤：MAF、缺失率、只保留双等位点，并转 PLINK
biopytools vcf-filter -i input.vcf -c chr1,chr2,chr3 \
    --maf 0.05 --max-missing 0.8 --biallelic-only --convert-format
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 VCF 文件路径 |
| `-c, --chr` | 染色体名称（支持逗号分隔多个） |

### 区间与输出 | Region & Output

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | 自动 | 输出 VCF 文件路径 |
| `-s, --start` | - | 起始位置 |
| `-e, --end` | - | 结束位置 |

### 质量与频率过滤 | Quality & Frequency

| 参数 | 描述 |
|------|------|
| `--maf` | 最小等位基因频率 |
| `--max-missing` | 最大缺失率 |
| `--quality-threshold` | 质量阈值 |
| `--min-depth` | 最小深度 |
| `--max-depth` | 最大深度 |

### 样本与位点选择 | Sample & Variant Selection

| 参数 | 描述 |
|------|------|
| `--keep-samples` | 保留样本名（逗号分隔） |
| `--remove-samples` | 移除样本名（逗号分隔） |
| `--keep-ids` | 保留变异位点 ID（逗号分隔） |
| `--remove-ids` | 移除变异位点 ID（逗号分隔） |
| `--biallelic-only` | 只保留双等位基因位点 |
| `--remove-indels` | 移除插入缺失变异 |

### 格式转换与验证 | Conversion & Validation

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--convert-format` | 关闭 | 使用 PLINK 进行格式转换 |
| `--plink-path` | `plink` | PLINK 可执行文件路径 |
| `--allow-extra-chr` | 开启 | 允许额外染色体（非人类命名） |
| `--skip-validation` | 开启 | 跳过输入验证以提高速度 |
| `--force-validation` | 关闭 | 强制执行输入验证 |
| `-v, --verbose` | 关闭 | 显示详细信息 |

（运行 `biopytools vcf-filter -h` 查看完整参数列表）

## 输出 | Output

```
output.vcf                     # 筛选后的 VCF（默认）
output.bed / output.bim / output.fam   # 启用 --convert-format 时的 PLINK 二进制格式
```

## 依赖 | Dependencies

- Python：标准库（实现核心筛选）
- [PLINK 1.9](https://www.cog-genomics.org/plink/)（仅当启用 `--convert-format` 时需要）

## 引用 | Citation

- Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 2015, 4: 7. doi:10.1186/s13742-015-0047-8
- Danecek P et al. The variant call format and VCFtools. *Bioinformatics*, 2011, 27(15): 2156-2158. doi:10.1093/bioinformatics/btr330

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [PLINK 1.9 文档](https://www.cog-genomics.org/plink/1.9/)
