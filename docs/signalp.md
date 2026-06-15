# SignalP 6.0 信号肽预测模块

**基于深度学习预测蛋白质信号肽及其类型 | Predict signal peptides and their types using deep learning**

## 功能概述 | Overview

SignalP模块封装了 [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP-6.0)，使用BERT蛋白质语言模型编码器和条件随机场（CRF）解码器预测蛋白质信号肽。模块自动检测SignalP所在的Conda环境并包装命令运行，并在原始输出基础上额外生成中文可读格式和R友好的纯TSV汇总文件，便于下游统计分析和Excel导入。

## 快速开始 | Quick Start

```bash
# 真核生物信号肽预测（默认）
biopytools signalp -i proteins.faa -o output_dir/

# 细菌/古菌预测（预测所有信号肽类型）
biopytools signalp -i proteins.faa -o output_dir/ --organism other

# 高精度slow模式 + 生成图表
biopytools signalp -i proteins.faa -o output_dir/ --mode slow --format all
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入FASTA文件（氨基酸序列，支持 `.fa/.fasta/.faa/.ffn/.fna`） |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--organism, -org` | `eukarya` | 生物类型：`eukarya`/`euk`（真核，仅Sec/SPI）、`other`（细菌/古菌，全部类型） |
| `--format, -fmt` | `txt` | 输出格式：`txt`、`png`、`eps`、`all`、`none` |
| `--mode, -m` | `fast` | 预测模式：`fast`（小模型）、`slow`（完整模型并行，需>14GB内存）、`slow-sequential`（串行，内存低但慢约6倍） |
| `--bsize, -bs` | `12` | 批处理大小 |
| `--write-procs, -wp` | `12` | 写入进程数 |
| `--torch-num-threads, -tt` | `12` | PyTorch线程数 |
| `--model-dir, -md` | `None` | 模型权重目录 |
| `--skip-resolve` | `False` | 跳过冲突解析步骤 |
| `--keep-plots` | `False` | 保留plot文件（覆盖默认清理行为） |
| `--signalp-path` | `~/miniforge3/envs/signalp6/bin/signalp6` | SignalP可执行文件路径 |

（运行 `biopytools signalp -h` 查看完整参数列表）

### 预测类型 | Prediction Types

| 代码 | 含义 |
|------|------|
| `SP` / `SPI` | Sec/SPI型经典信号肽（分泌蛋白） |
| `LIPO` / `SPII` | Sec/SPII型脂蛋白信号肽 |
| `TAT` | Tat/SPI型双精氨酸信号肽 |
| `TATLIPO` | Tat/SPII型脂蛋白 |
| `PILIN` / `SPIII` | Sec/SPIII型菌毛蛋白 |
| `OTHER` | 无信号肽 |

## 输出 | Output

```
output_dir/
├── prediction_results.txt           # 原始预测结果
├── processed_entries.fasta          # 切除信号肽后的成熟蛋白序列
├── output.gff3                      # 信号肽GFF3注释
├── region_output.gff3               # 区域GFF3注释
├── output.json                      # JSON格式结果
├── prediction_results_readable.txt  # 中文可读格式（含统计和列说明）⭐
├── signalp_summary.tsv              # R友好的纯TSV汇总 ⭐
└── signalp_prediction.log           # 运行日志
```

`signalp_summary.tsv` 列：`蛋白质ID  预测代码  Has_SP  SP概率  SP起始  SP终止  切割位点起始  切割位点终止  切割置信度`，可直接用 `read.delim()` 读入R。

## 依赖 | Dependencies

- SignalP 6.0（需单独安装，参考官网）
- PyTorch（SignalP 6.0自带依赖）
- 默认程序路径：`~/miniforge3/envs/signalp6/bin/signalp6`
- slow模式需要完整模型文件（>14GB内存）

## 引用 | Citation

- Teufel, F. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nature Biotechnology*, 2022, 40:1023-1025.

## 相关链接 | References

- [SignalP 6.0 官网](https://services.healthtech.dtu.dk/service.php?SignalP-6.0)
- [SignalP 6.0 GitHub](https://github.com/fteufel/signalp-6.0)
- [项目主页](https://github.com/lixiang117423/biopytools)
