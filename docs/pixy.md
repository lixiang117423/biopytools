# Pixy群体遗传统计 | Pixy Population Genetics Statistics

**基于pixy计算pi/dxy/fst/Watterson theta/Tajima's D, 无偏且支持不变位点 | Unbiased population genetics statistics supporting invariant sites**

## 功能概述 | Overview

pixy 模块封装了 pixy 工具, 用于计算群体遗传学统计量。相比传统工具(如VCFtools), pixy 通过包含不变位点(invariant sites)避免"除以零"导致的fst高估问题, 是当前计算pi、dxy、fst的金标准工具。支持滑窗和全基因组两种模式, 适用于群体遗传结构分析、群体分化评估、选择信号检测等场景。

## 快速开始 | Quick Start

```bash
# 滑窗分析(100kb窗口)
biopytools pixy -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000

# 指定统计量
biopytools pixy -i variants.vcf.gz -p pops.txt -o out --stats pi,fst,dxy
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --vcf-file` | VCF文件(需bgzip压缩并建立tabix索引) |
| `-p, --pop-file` | 群体文件(两列:样本ID 群体名) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./pixy_output` | 输出目录 |
| `--stats` | `pi,dxy,fst,watterson,tajima` | 要计算的统计量,逗号分隔 |
| `-w, --window-size` | `None` | 窗口大小bp(不设则全基因组) |
| `-b, --bed-file` | `None` | BED文件定义窗口 |
| `-s, --sites-file` | `None` | 位点文件(只计算特定位点) |
| `--min-samples` | `0` | 每个群体最小样本数 |
| `--max-missing` | `1.0` | 最大缺失率 |
| `--min-maf` | `0.0` | 最小等位基因频率 |
| `-c, --chromosomes` | `None` | 指定染色体列表,逗号分隔 |
| `-t, --threads` | `12` | 线程数 |
| `--pixy-path` | `pixy` | pixy可执行文件路径 |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` | conda环境路径 |
| `--bypass-invariant-check` | `False` | 强制绕过不变位点检查(默认自动检测) |
| `--keep-intermediate` | `False` | 保留中间文件 |

(运行 `biopytools pixy -h` 查看完整参数列表)

> **注意**: pixy必须指定`-w`(窗口)、`-b`(BED文件)或`-s`(位点文件)之一。

## 输出 | Output

- `pixy_{stat}_{pop}_{window}.txt`: 每个统计量的滑窗/全基因组结果
- `pixy.log`: 运行日志

## 依赖 | Dependencies

- **pixy**: 无偏群体遗传统计工具 (https://github.com/ksamuk/pixy)
- **tabix/bgzip**: VCF索引 (htslib)
- **conda环境**: 推荐在pixy专用环境运行

## 引用 | Citation

- Korunes K.L., Samuk K. (2021) pixy: Unbiased estimation of nucleotide diversity within and between populations. Molecular Ecology Resources. 21(3):878-889.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
