# LAI 组装质量指数计算 | LTR Assembly Index (LAI) Calculator

**基于 LTR 逆转座子完整度评估基因组组装连续性 | Evaluate assembly continuity via LTR retrotransposon integrity**

## 功能概述 | Overview

LAI（LTR Assembly Index）是通过计算长末端重复逆转座子（LTR-RT）的完整度来衡量基因组组装质量的指标。完整 LTR-RT 的比例越高，说明组装越连续。LAI 已成为植物基因组组装质量评估的事实标准之一，与 N50、BUSCO 等互补。

本模块串联了 LTRharvest 与 LTR_finder 并行版进行 LTR 候选识别，再用 LTR_retriever 进行严格筛选和 LAI 计算。流程支持完整运行，也可按阶段拆分（仅识别 / 仅筛选 / 仅计算 LAI），并自动跳过已完成的步骤以方便断点续跑。三个工具分别运行在独立的 conda 环境中。

## 快速开始 | Quick Start

```bash
# 基本用法：完整流程
biopytools lai -i genome.fa -o output_dir

# 仅计算 LAI（已有 LTR_retriever 输出）
biopytools lai -i genome.fa -o output_dir -m calculate

# 指定线程和自定义 conda 环境
biopytools lai -i asm.fa -o lai_results -t 32 \
    --conda-harvest ~/miniforge3/envs/my_ltr_harvest
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 基因组 FASTA 文件 |
| `-o, --output` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `-m, --mode` | `full` | 运行模式（full/harvest/retrieve/calculate）|
| `--skip-completed / --no-skip-completed` | 跳过 | 是否跳过已完成的步骤 |
| `--conda-harvest` | `~/miniforge3/envs/ltr_harvest_parallel_v.1.2` | LTR_harvest 并行版 conda 环境 |
| `--conda-finder` | `~/miniforge3/envs/ltr_finder_parallel_v.1.3` | LTR_finder 并行版 conda 环境 |
| `--conda-retriever` | `~/miniforge3/envs/ltr_retriever_v.3.0.1` | LTR_retriever conda 环境 |

运行模式说明：
- `full`：完整流程（识别 → 筛选 → 计算）
- `harvest`：仅运行 LTR 候选识别
- `retrieve`：仅运行 LTR_retriever 筛选
- `calculate`：仅运行 LAI 计算

（运行 `biopytools lai -h` 查看完整参数列表）

## 输出 | Output

```
output_dir/
├── *.out.harvest.scatter      # LTRharvest 候选
├── *.out.finder.scatter       # LTR_finder 候选
├── *.pass.list                # LTR_retriever 通过列表
├── *.LTRlib                   # LTR 库
├── *.out.LAI                  # LAI 指数结果
├── *.LAI_summary              # LAI 汇总
└── lai.log                    # 运行日志
```

## 依赖 | Dependencies

- LTR_harvest 并行版（`--conda-harvest`）
- LTR_finder 并行版（`--conda-finder`）
- LTR_retriever（`--conda-retriever`，含 LAI 计算）
- CD-HIT、BLAST+（LTR_retriever 依赖）

## 引用 | Citation

- Ou, S. et al. Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. *Genome Biology* 20, 275 (2019).
- Ou, S. & Jiang, N. LTR_retriever: A highly accurate and sensitive program for identification of long terminal repeat retrotransposons. *Plant Physiology* 176, 1410-1422 (2018).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [LTR_retriever 官方仓库](https://github.com/oushujun/LTR_retriever)
