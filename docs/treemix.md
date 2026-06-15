# TreeMix群体历史 | TreeMix Population History & Gene Flow

**基于TreeMix推断群体历史与基因流, 含OptM最优m值估计和可视化 | Infer population history and gene flow with OptM optimization**

## 功能概述 | Overview

treemix 模块封装了 TreeMix 工具, 用于推断群体分化历史和基因流(gene flow)事件。TreeMix基于等位基因频率协方差构建群体迁移图(migration graph)。本模块提供三个子命令: prepare (VCF转TreeMix格式)、run (运行m值扫描+OptM+绘图)、all (一键完整流程), 自动完成LD过滤、频率计算、TreeMix格式转换、迁移边扫描、最优m值确定及树图/残差图绘制。

## 快速开始 | Quick Start

```bash
# 完整流程
biopytools treemix all -i input.vcf.gz -o treemix_output --cluster cluster.txt

# 仅运行TreeMix (已有.frq.gz格式输入)
biopytools treemix run -i input.frq.gz -o output --root outgroup --m-max 5
```

## 子命令 | Subcommands

### prepare - 输入准备 | Input Preparation

| 参数 | 描述 |
|------|------|
| `-i, --input` | VCF输入文件 |
| `-o, --output` | 输出目录 |
| `--cluster` | 样本分组文件(sample_id population) |
| `--pop-delimiter` | 自动推断群体时的分隔符(默认`_`) |
| `--ld-window` | LD窗口SNP数(默认50) |
| `--ld-step` | LD步长(默认10) |
| `--ld-r2` | LD r2阈值(默认0.2) |

### run - 运行TreeMix | Run TreeMix

| 参数 | 描述 |
|------|------|
| `-i, --input` | TreeMix输入文件(.frq.gz) |
| `-o, --output` | 输出目录 |
| `--m-max` | 测试m=0..m_max(默认10) |
| `-m` | 指定m值(单独运行或绘图) |
| `--root` | 外群群体名 |
| `--bootstrap` | bootstrap次数(默认1000) |
| `--replicates` | 每个m值的重复次数(默认10) |
| `-k` | SNP block大小(默认500) |
| `--noss` | 关闭样本量校正 |
| `--se` | 计算迁移权重标准误 |
| `-t, --threads` | 线程数(默认12) |

### all - 完整流程 | Full Pipeline

参数为 prepare + run 的并集, 同时支持 `--treemix-path`、`--plink-path`、`--bcftools-path`、`--r-path`、`--plotting-funcs`。

(运行 `biopytools treemix <subcommand> -h` 查看完整参数列表)

## 输出 | Output

- `*.frq.gz`: TreeMix输入格式(由plink生成)
- `m{N}/`: 每个m值的TreeMix输出(树图、残差、协方差)
- `OptM_result.csv`: OptM统计量(用于确定最优m)
- `treemix_tree_m{optimal}.png`: 最优m值的树图
- `residual_m{optimal}.png`: 残差图

## 依赖 | Dependencies

- **TreeMix**: 群体历史推断 (https://bitbucket.org/nygcresearch/treemix)
- **PLINK**: 频率计算和LD过滤
- **BCFtools**: VCF处理
- **R + OptM包**: 最优m值估计 (https://github.com/oliverfrainey/OptM)
- **plotting_funcs.R**: TreeMix树绘图脚本

## 引用 | Citation

- Pickrell J.K., Pritchard J.K. (2012) Inference of population splits and mixtures from genome-wide allele frequency data. PLoS Genetics. 8(11):e1002967.
- Fitak R.R. (2021) OptM: estimating the optimal number of migration edges on population trees using Treemix. Biology Methods and Protocols. 6(1):bpab017.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
