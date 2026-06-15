# SignalP 6.0信号肽预测工具|SignalP 6.0 Signal Peptide Prediction Tool

## 功能简介|Introduction

SignalP 6.0是基于深度学习的信号肽预测工具，使用BERT蛋白质语言模型编码器和条件随机场(CRF)解码器来预测蛋白质的信号肽及其类型。

SignalP 6.0 is a deep learning-based signal peptide prediction tool that uses a BERT protein language model encoder and conditional random field (CRF) decoder to predict signal peptides and their types.

## 安装|Installation

SignalP 6.0需要单独安装，请参考官方文档：
https://services.healthtech.dtu.dk/service.php?SignalP-6.0

本模块默认安装路径：
`/share/org/YZWL/yzwl_lixg/miniforge3/envs/signalp6/`

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 使用biopytools命令|Use biopytools command
biopytools signalp -i proteins.faa -o output_dir

# 使用Python模块|Use Python module
python -m biopytools.signalp -i proteins.faa -o output_dir
```

### 参数说明|Parameters

#### 必需参数|Required Parameters

| 参数|说明|
|------|------|
| `-i, --input` | 输入FASTA文件(氨基酸序列)|Input FASTA file (amino acid sequences) |
| `-o, --output-dir` | 输出目录|Output directory |

#### 可选参数|Optional Parameters

| 参数|默认值|说明|
|------|-------|------|
| `--organism, -org` | `eukarya` | 生物类型: `eukarya`, `other`, `euk`|Organism type |
| `--format, -fmt` | `txt` | 输出格式: `txt`, `png`, `eps`, `all`, `none`|Output format |
| `--mode, -m` | `fast` | 预测模式: `fast`, `slow`, `slow-sequential`|Prediction mode |
| `--bsize, -bs` | `12` | 批处理大小|Batch size |
| `--write-procs, -wp` | `12` | 写入进程数|Number of write processes |
| `--torch-num-threads, -tt` | `12` | PyTorch线程数|PyTorch threads |
| `--signalp-path` | 见下方|SignalP程序路径|SignalP program path |
| `--model-dir, -md` | `None` | 模型权重目录|Model weights directory |
| `--skip-resolve` | `False` | 跳过冲突解析|Skip conflict resolution |

### 使用示例|Examples

```bash
# 1. 基本预测|Basic prediction
biopytools signalp -i proteins.faa -o output

# 2. 细菌/古菌预测|Bacteria/Archaea prediction
biopytools signalp -i proteins.faa -o output --organism other

# 3. 使用slow模式|Use slow mode
biopytools signalp -i proteins.faa -o output --mode slow

# 4. 生成图表|Generate plots
biopytools signalp -i proteins.faa -o output --format all

# 5. 自定义线程数|Customize threads
biopytools signalp -i proteins.faa -o output --torch-num-threads 16 --write-procs 4
```

## 输出文件|Output Files

### 必需输出文件|Required Output Files

| 文件名|说明|
|--------|------|
| `prediction_results.txt` | 预测结果汇总|Prediction summary (ID, 预测类型, 概率, 剪切位点)|
| `processed_entries.fasta` | 剪切信号肽后的成熟蛋白序列|Mature proteins with SP removed |
| `output.gff3` | 信号肽GFF3格式注释|Signal peptide positions in GFF3 format |
| `region_output.gff3` | 区域GFF3格式|Signal peptide regions in GFF3 format |
| `output.json` | JSON格式结果|Results in JSON format |

### 可选输出文件|Optional Output Files

| 文件名|说明|
|--------|------|
| `SEQUENCE_plot.txt` | 每个序列的位置预测表|Per-position predictions |
| `SEQUENCE_plot.png` | 预测概率图表|Probability plot |
| `SEQUENCE_plot.eps` | EPS格式图表|EPS format plot |

## 预测类型|Prediction Types

SignalP 6.0可以预测以下类型的信号肽：

| 类型|说明|
|------|------|
| `SP` | Sec/SPI型信号肽|Sec/SPI signal peptide |
| `LIPO` | Sec/SPII型脂蛋白信号肽|Sec/SPII lipoprotein signal peptide |
| `TAT` | Tat/SPI型双精氨酸信号肽|Tat/SPI twin-arginine signal peptide |
| `TATLIPO` | Tat/SPII型脂蛋白|Tat/SPII lipoprotein |
| `PILIN` | Sec/SPIII型菌毛蛋白|Sec/SPIII pilin |
| `OTHER` | 无信号肽|No signal peptide |

## 预测模式说明|Mode Description

| 模式|说明|内存需求|速度|
|------|------|---------|-----|
| `fast` | 小模型，快速预测|低|快|
| `slow` | 完整模型，并行预测|>14GB|中等|
| `slow-sequential` | 完整模型，串行预测|低|慢(约为fast的6倍)|

## 注意事项|Notes

1. **输入文件格式**：必须是氨基酸序列FASTA文件（.faa, .fa, .fasta等）
2. **生物类型**：
   - `other`：细菌、古菌等（预测所有类型）|Bacteria, archaea (predict all types)
   - `eukarya`/`euk`：真核生物（仅预测Sec/SPI）|Eukaryotes (Sec/SPI only)
3. **模式选择**：slow模式需要提前安装完整模型文件
4. **GPU支持**：如有GPU可转换模型以加速预测（参考SignalP官方文档）

## 参考资料|References

- SignalP 6.0官网: https://services.healthtech.dtu.dk/service.php?SignalP-6.0
- GitHub仓库: https://github.com/fteufel/signalp-6.0
- 论文: [待补充]

## 版本历史|Version History

| 版本|日期|说明|
|------|------|------|
| 1.0.0 | 2026-02-04 | 初始版本|Initial release |
