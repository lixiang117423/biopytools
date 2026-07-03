# DeepTMHMM 跨膜螺旋/信号肽预测模块

**预测蛋白质跨膜螺旋（TMH）与信号肽 | Predict transmembrane helices and signal peptides in proteins**

## 功能概述 | Overview

DeepTMHMM模块封装了 [DeepTMHMM 1.0](https://dtu.dk/bioinformatics/DeepTMHMM)（DTU，基于深度语言模型），预测蛋白质序列中的跨膜螺旋（TMH）、β-桶（β-barrel）和信号肽（Signal Peptide）。模块自动调用 `predict.py`，把原始输出整理成结构化 TSV（带英文列名，便于 R/Excel 处理），同时保留 `3line` 拓扑文件和 `GFF3` 供下游分析。支持断点续传：主输出齐全则自动跳过。

> 与 `tmhmm` 模块的区别：`tmhmm` 封装的是经典的 TMHMM 2.0c（仅 HMM 跨膜螺旋）；`deeptmhmm` 封装的是更新的 DeepTMHMM 1.0（深度学习，同时预测信号肽，精度更高）。

## 快速开始 | Quick Start

```bash
# 基本预测
biopytools deeptmhmm -i proteins.fa -o output_dir/

# 指定输出前缀
biopytools deeptmhmm -i proteins.fa -o output_dir/ --prefix sample1

# 指定conda环境与安装目录
biopytools deeptmhmm -i proteins.fa -o output_dir/ \
    --conda-env deeptmhmm_v.1.0 \
    --deeptmhmm-dir ~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入蛋白质FASTA文件 |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--prefix` | 输入文件名（去扩展名） | 输出文件前缀 |
| `--conda-env` | `deeptmhmm_v.1.0` | DeepTMHMM所在conda环境名 |
| `--deeptmhmm-dir` | `~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0` | DeepTMHMM安装目录（含 `predict.py`） |

（运行 `biopytools deeptmhmm -h` 查看完整参数列表）

### 工具路径配置

- `--deeptmhmm-dir` 按以下优先级查找：命令行参数 → 环境变量 `DEEPTMHMM_DIR` → `~/.config/biopytools/config.yml` → 默认路径
- `--conda-env` 可用环境变量 `DEEPTMHMM_ENV` 覆盖
- conda调用自动加 `--no-capture-output`，predict.py进度实时流式输出到作业 `.out`（§13.2.0）

## 输出 | Output

```
output_dir/
├── {prefix}_deeptmhmm_summary.tsv      # 整理后的结构化TSV ⭐
├── {prefix}_deeptmhmm_topologies.3line # 三行拓扑（>id | TYPE / 序列 / 标注）
├── {prefix}_deeptmhmm_tmr.gff3         # TMH/信号肽区段GFF3
├── {prefix}_deeptmhmm_results.md       # DeepTMHMM原始markdown报告
└── deeptmhmm.log                       # 运行日志
```

`{prefix}_deeptmhmm_summary.tsv` 列说明：

| 列名 | 描述 |
|------|------|
| `ID` | 蛋白质ID |
| `Length` | 蛋白长度（氨基酸数） |
| `Protein_Type` | DeepTMHMM判定的蛋白类型（如 `TMhelix`、`SignalPeptide`、无则空） |
| `Pred_TMHs` | 预测的跨膜螺旋数量 |
| `Signal_Peptide` | 信号肽（`no` 或 `yes (start-end)`） |
| `TM_Regions` | 跨膜螺旋区段（如 `10-32;45-67`，无则 `-`） |

运行结束后日志会输出统计：蛋白总数、`0个TMH` / `>=1个TMH` 数量、含信号肽数量。

## 依赖 | Dependencies

- **DeepTMHMM 1.0**（需单独安装，学术用户从 DTU 申请 Academic License）
- conda环境：`deeptmhmm_v.1.0`（环境内含 `predict.py` 运行所需的 Python/PyTorch/biopython）
- 安装目录需含 `predict.py`

## 引用 | Citation

- Hallgren, J., Tsirigos, K.D., Pedersen, M.D. et al. DeepTMHMM predicts alpha and beta transmembrane proteins using deep language models. *Nature Machine Intelligence*, 2022, 4: 598–609.

## 相关链接 | References

- [DeepTMHMM (DTU)](https://dtu.dk/bioinformatics/DeepTMHMM)
- [项目主页](https://github.com/lixiang117423/biopytools)
