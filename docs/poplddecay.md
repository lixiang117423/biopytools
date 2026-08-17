# LD衰减分析 | Linkage Disequilibrium Decay Analysis (PopLDdecay)

**基于PopLDdecay计算连锁不平衡(LD)衰减, 支持R^2/D'度量和多种绘图方法 | LD decay analysis based on PopLDdecay**

## 功能概述 | Overview

poplddecay 模块封装了 PopLDdecay 工具, 用于计算群体连锁不平衡(LD)衰减。LD衰减反映重组率历史和群体遗传漂变模式, 是群体遗传学、GWAS标记密度评估、关联分析有效性判断的重要指标。本工具支持 VCF 和 Genotype 两种输入格式, 可指定子群体、自定义bin大小、自动推荐LD阈值并绘制衰减曲线图。

## 快速开始 | Quick Start

```bash
# VCF标准分析
biopytools poplddecay -i variants.vcf -o output_prefix

# 按子群体分析
biopytools poplddecay -i variants.vcf -o output_prefix -s subpop_list.txt
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入VCF或Genotype文件 |
| `-o, --output` | 输出文件前缀 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --type` | `vcf` | 输入文件类型(vcf/genotype) |
| `-d, --max-dist` | `300` | 最大距离(kb) |
| `-m, --min-maf` | `0.005` | 最小等位基因频率 |
| `--max-het` | `0.88` | 最大杂合率 |
| `--max-miss` | `0.25` | 最大缺失率 |
| `-s, --subpop` | `None` | 子群体样本列表文件 |
| `--out-type` | `1` | 输出类型(1:R^2, 2:R^2&D', 3:Pairwise LD) |
| `--bin1` | `10` | 短距离bin大小 |
| `--bin2` | `100` | 长距离bin大小 |
| `--break-point` | `100` | 短/长距离分界点 |
| `--measure` | `r2` | LD度量(r2/D/both) |
| `--method` | `MeanBin` | 绘图方法(MeanBin/HW/MedianBin/PercentileBin) |
| `--no-plot` | `False` | 不绘制图像 |
| `--no-recommend-threshold` | `False` | 不推荐LD阈值 |

(运行 `biopytools poplddecay -h` 查看完整参数列表)

## 输出 | Output

- `{prefix}.LD.gz`: LD统计原始结果
- `{prefix}.LD.decay.gz`: LD衰减数据
- `{prefix}.LD.decay.png/pdf`: LD衰减曲线图

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 | Path | 输入VCF或Genotype文件｜Input VCF or Genotype file |
| `-o, --output` | 必填 | Path | 输出文件前缀｜Output file prefix |
| `-t, --type` | `vcf` | vcf/genotype | 输入文件类型｜Input file type (vcf/genotype) [default: vcf] |
| `-d, --max-dist` | `300` | int | 最大距离(kb)｜Max distance in kb [default: 300] |
| `-m, --min-maf` | `0.005` | float | 最小等位基因频率｜Min minor allele frequency [default: 0.005] |
| `--max-het` | `0.88` | float | 最大杂合率｜Max heterozygous rate [default: 0.88] |
| `--max-miss` | `0.25` | float | 最大缺失率｜Max missing rate [default: 0.25] |
| `-s, --subpop` | — | Path | 子群体样本列表文件｜Subpopulation sample list file |
| `--out-type` | `1` | int | 输出类型(1:R^2, 2:R^2&D', 3:Pairwise LD)｜Output type [default: 1] |
| `--bin1` | `10` | int | 短距离bin大小｜Bin size for short distance [default: 10] |
| `--bin2` | `100` | int | 长距离bin大小｜Bin size for long distance [default: 100] |
| `--break-point` | `100` | int | 短/长距离分界点｜Break point [default: 100] |
| `--max-x` | — | int | 最大X坐标(kb)｜Max X coordinate in kb |
| `--measure` | `r2` | r2/D/both | LD度量方法｜LD measure method [default: r2] |
| `--method` | `MeanBin` | MeanBin/HW/MedianBin/PercentileBin | 绘图方法｜Plotting method [default: MeanBin] |
| `--percentile` | `0.5` | float | 百分位数｜Percentile for PercentileBin [default: 0.5] |
| `--no-plot` | — |  | 不绘制图像｜Do not plot figure |
| `--no-recommend-threshold` | — |  | 不推荐LD阈值｜Do not recommend LD threshold |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **PopLDdecay**: LD衰减分析工具 (https://github.com/BGI-shenzhen/PopLDdecay)

## 引用 | Citation

- Zhang C., Dong S.S., Xu J.Y., et al. (2019) PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. Bioinformatics. 35(10):1786-1788.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
