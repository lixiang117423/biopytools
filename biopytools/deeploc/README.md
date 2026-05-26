# DeepLoc 2.1蛋白质亚细胞定位预测工具|DeepLoc 2.1 Protein Subcellular Localization Prediction Tool

## 功能简介|Introduction

DeepLoc 2.1是基于蛋白质语言模型的多标签亚细胞定位预测工具，可以：
- 预测10种亚细胞定位（细胞核、细胞质、细胞外、线粒体、细胞膜、内质网、叶绿体、高尔基体、溶酶体/液泡、过氧化物酶体）
- 预测4种膜蛋白类型（外周膜、跨膜、脂质锚定、可溶性）
- 预测排序信号

DeepLoc 2.1 is a multi-label subcellular localization predictor based on protein language models. It can:
- Predict 10 subcellular localizations (Nucleus, Cytoplasm, Extracellular, Mitochondrion, Cell membrane, ER, Chloroplast, Golgi, Lysosome/Vacuole, Peroxisome)
- Predict 4 membrane protein types (Peripheral, Transmembrane, Lipid-anchored, Soluble)
- Predict sorting signals

## 安装|Installation

DeepLoc 2.1通过Singularity容器运行，需要单独安装：
- Singularity版本: 3.8.7+
- 镜像位置: `/share/org/YZWL/yzwl_lixg/software/singularity/deeploc2.1_latest.sif`
- 数据库位置: `/share/org/YZWL/yzwl_lixg/software/deeploc/database`

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 使用biopytools命令|Use biopytools command
biopytools deeploc -i proteins.faa -o output_dir

# 使用Python模块|Use Python module
python -m biopytools.deeploc -i proteins.faa -o output_dir
```

### 参数说明|Parameters

#### 必需参数|Required Parameters

| 参数|说明|
|------|------|
| `-i, --input` | 输入FASTA文件(蛋白质序列)|Input FASTA file (protein sequences) |
| `-o, --output-dir` | 输出目录|Output directory |

#### 可选参数|Optional Parameters

| 参数|默认值|说明|
|------|-------|------|
| `--model` | `Fast` | 预测模型: `Fast`(ESM1b), `Accurate`(ProtT5)|Prediction model |
| `--device` | `cpu` | 计算设备: `cpu`, `cuda`, `mps`|Compute device |
| `--singularity-image` | 见下方|Singularity镜像路径|Singularity image path |
| `--database-dir` | 见下方|数据库目录|Database directory |
| `--singularity-exec` | 见下方|Singularity可执行文件路径|Singularity executable path |
| `--plot` | `False` | 绘制attention图|Plot attention values |

### 使用示例|Examples

```bash
# 1. 基本预测|Basic prediction
biopytools deeploc -i proteins.faa -o output

# 2. 使用Accurate模型(高精度)|Use Accurate model
biopytools deeploc -i proteins.faa -o output --model Accurate

# 3. 使用GPU加速|Use GPU acceleration
biopytools deeploc -i proteins.faa -o output --device cuda

# 4. 绘制attention图|Plot attention
biopytools deeploc -i proteins.faa -o output --plot
```

## 输出文件|Output Files

### 主要输出文件|Main Output File

| 文件名|说明|
|--------|------|
| `{input}_deeploc2.tsv` | 预测结果表格|Prediction results table |

### 表格格式|Table Format

TSV文件包含以下列：
1. **蛋白质ID** - Protein ID
2. **预测定位** - Predicted localization(s)
3. **排序信号** - Predicted sorting signal(s)
4. **各定位概率** - Probability for each localization (10列)

### 10个亚细胞定位|10 Subcellular Localizations

| 编号|定位|英文|
|----|------|------|
| 1 | 细胞核|Nucleus |
| 2 | 细胞质|Cytoplasm |
| 3 | 细胞外|Extracellular |
| 4 | 线粒体|Mitochondrion |
| 5 | 细胞膜|Cell membrane |
| 6 | 内质网|Endoplasmic reticulum |
| 7 | 叶绿体|Chloroplast |
| 8 | 高尔基体|Golgi apparatus |
| 9 | 溶酶体/液泡|Lysosome/Vacuole |
| 10 | 过氧化物酶体|Peroxisome |

### 4种膜蛋白类型|4 Membrane Protein Types

| 类型|英文|
|------|------|
| 外周膜蛋白|Peripheral membrane protein |
| 跨膜蛋白|Transmembrane protein |
| 脂质锚定蛋白|Lipid anchored protein |
| 可溶性蛋白|Soluble protein |

### 可选输出文件|Optional Output Files

| 文件名|说明|
|--------|------|
| `{input}_deeploc2_plot.png` | Attention可视化图|Attention visualization |
| `{input}_deeploc2_plot.txt` | Attention值文本文件|Attention values text file |

## 预测模式说明|Model Description

| 模型|说明|参数量|速度|精度|内存需求|
|------|------|--------|-----|-----|---------|
| `Fast` | ESM1b模型|650M|快|稍低|低|
| `Accurate` | ProtT5-XL-Uniref50模型|3B|慢|高|~32GB|

**建议**：
- `Fast`: 大规模预测、快速筛选
- `Accurate`: 高精度预测、少量蛋白质

## 注意事项|Notes

1. **输入文件格式**：必须是蛋白质序列FASTA文件（.faa, .fa, .fasta等）
2. **多标签预测**：一个蛋白质可以预测多个定位
3. **Accurate模型**：需要较多内存（~32GB），个人电脑不推荐
4. **GPU支持**：需要CUDA支持的NVIDIA GPU
5. **容器运行**：通过Singularity容器运行，无需安装Python依赖

## 参考资料|References

- DeepLoc 2.1官网: https://services.healthtech.dtu.dk/services/DeepLoc-2.1/
- 论文: https://doi.org/10.1093/nar/gkae237
- DeepLoc 2.0论文: https://doi.org/10.1093/nar/gkac278

## 版本历史|Version History

| 版本|日期|说明|
|------|------|------|
| 1.0.0 | 2026-02-04 | 初始版本|Initial release |
