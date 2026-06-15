# Phytophthora效应子鉴定 | Phytophthora Effector Identification

**鉴定卵菌（Phytophthora等）蛋白中的 RxLR / CRN / NLP 等效应子候选 | Identify RxLR, CRN, NLP and other candidate effectors from oomycete proteins**

## 功能概述 | Overview

针对预测好的蛋白质集合，先用 SignalP 预测分泌信号肽，再通过 HMMER 搜索各效应子家族 HMM 模型（内置 RxLR、CRN、NLP、Protease、SCP、elicitin、YxSL，以及可选的 WY 结构域），输出每类效应子的候选蛋白与注释 TSV。支持单文件或多样本目录批量运行，并提供子命令合并多样本结果。

- 一键流程：SignalP 信号肽预测 -> HMM 搜索 + 基序注释
- 内置七类效应子 HMM 模型，无需用户额外准备
- 支持 SignalP 3.0、6.0 或同时运行两个版本
- 支持单样本 / 多样本目录自动识别，并提供 `merge` 子命令汇总

## 快速开始 | Quick Start

```bash
# 鉴定效应子（单文件或目录均可）
biopytools phyto-effector -i proteins.fa -o ./effector_out

# 多样本：直接传入目录，结果将分样品输出
biopytools phyto-effector -i ./proteins_dir/ -o ./effector_out

# 合并多样本运行结果
biopytools phyto-effector merge -i ./effector_out -o ./effector_merged
```

## 参数说明 | Parameters

### run 子命令（默认）| run subcommand

#### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入蛋白质 FASTA 文件或目录（目录中自动扫描 `.fa` / `.fasta` / `.faa` 等） |

#### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./phyto_effector_output` | 输出目录 |
| `--skip-signalp` | False | 跳过 SignalP 预测（已有结果时使用） |
| `--signalp-version` | `both` | SignalP 版本：`3` / `6` / `both` |
| `--signalp-mode` | `slow-sequential` | SignalP 运行模式：`fast` / `slow` / `slow-sequential` |
| `--organism` | `eukarya` | 生物类型：`eukarya` / `other` |
| `--signalp3-sprob-threshold` | `0.9` | SignalP 3.0 的 HMM Sprob 阈值 |
| `--use-wy-domain` | False | 同时搜索 WY 结构域（PF18634） |
| `--score-threshold` | `0.0` | HMM score 阈值 |
| `-e, --evalue` | `1e-5` | E-value 阈值（已弃用，保留兼容） |
| `-t, --threads` | `12` | 线程数 |

#### 程序路径可选覆盖 | Tool path overrides

| 参数 | 默认值 |
|------|--------|
| `--signalp-path` | `~/miniforge3/envs/signalp6/bin/signalp6` |
| `--signalp3-path` | `~/miniforge3/envs/signalp_v.3.0b/bin/signalp` |
| `--hmmsearch-path` | `~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch` |

每类效应子还可通过 `--rxlr-hmm` / `--crn-hmm` / `--nlp-hmm` / `--protease-hmm` / `--scp-hmm` / `--elicitin-hmm` / `--yxsl-hmm` / `--rxlr-wy-hmm` 替换默认内置 HMM。

（运行 `biopytools phyto-effector -h` 或 `biopytools phyto-effector run -h` 查看完整参数列表）

### merge 子命令 | merge subcommand

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | 必需 | 含多样本输出目录的父目录 |
| `-o, --output-dir` | `./phyto_effector_merged` | 合并结果输出目录 |

## 输出 | Output

输出目录典型结构（每类效应子独立子目录）：

- `01_signalp/`：SignalP 信号肽预测结果
- `02_rxlr/`、`06_candidates/rxlr_candidates.tsv`：RxLR 候选（含 WY 域注释）
- `03_candidates/crn_candidates.tsv`：CRN 候选
- `04_*` ~ `08_*`：NLP / Protease / SCP / elicitin / YxSL 候选
- `00_pipeline_info/`：软件版本与流程信息
- `99_logs/`：运行日志

`merge` 子命令会生成 `<type>_all.tsv`（如 `rxlr_all.tsv`），按样品汇总所有候选。

## 依赖 | Dependencies

- `signalp`（6.0）和/或 `signalp` 3.0b（信号肽预测）
- `hmmsearch`（HMMER 套件，效应子家族搜索）
- 可选 `tmhmm`（跨膜区预测）
- Python 包：`pandas`（merge 子命令）

## 引用 | Citation

效应子家族 HMM 来源随包内置。常用引用：
- Petersen TN, et al. SignalP 4.0: discriminating signal peptides from transmembrane regions. *Nature Methods*. 2011, 8(10):785-786.
- Eddy SR. Accelerated Profile HMM Searches. *PLoS Computational Biology*. 2011, 7(10):e1002195.

## 相关链接 | References

- [HMMER 项目](http://hmmer.org/)
- [SignalP 服务](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
- [项目主页](https://github.com/lixiang117423/biopytools)
