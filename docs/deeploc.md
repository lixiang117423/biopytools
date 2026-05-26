# DeepLoc蛋白质亚细胞定位预测模块

**蛋白质亚细胞定位预测工具，基于DeepLoc 2.1 | Protein Subcellular Localization Prediction Tool Based on DeepLoc 2.1**

## 功能概述 | Overview

DeepLoc模块是对DeepLoc 2.1工具的封装和增强，通过Singularity容器运行DeepLoc 2.1，预测蛋白质的亚细胞定位。该模块在原始DeepLoc输出基础上，自动生成中文翻译的TSV格式结果和带样式的Excel文件，大幅提升结果的可读性和易用性。

## 主要特性 | Key Features

- **两种预测模型**: Fast（ESM1b，快速）和Accurate（ProtT5，高精度）
- **多设备支持**: CPU、CUDA（GPU）、MPS（Apple Silicon）
- **中文结果翻译**: 10种亚细胞定位、7种信号类型、4种膜蛋白类型完整翻译
- **多格式输出**: 原始TSV、可读TSV（中文表头）、Excel（带样式）
- **Attention可视化**: 可选生成注意力权重图
- **Singularity容器化**: 隔离运行环境，避免依赖冲突
- **自动Conda检测**: 智能检测Singularity所在Conda环境并自动包装命令

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用Fast模型（默认）进行预测
biopytools deeploc -i proteins.faa -o output_dir

# 使用Accurate模型（高精度）
biopytools deeploc -i proteins.faa -o output_dir --model Accurate

# 使用GPU加速
biopytools deeploc -i proteins.faa -o output_dir --device cuda

# 绘制attention图
biopytools deeploc -i proteins.faa -o output_dir --plot
```

### 完整示例 | Full Example

```bash
# 使用Accurate模型 + GPU + attention图
biopytools deeploc \
    -i /path/to/proteins.faa \
    -o /path/to/output \
    --model Accurate \
    --device cuda \
    --plot
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTA文件（蛋白质序列） | `-i proteins.faa` |
| `-o, --output-dir` | 输出目录 | `-o output_dir` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--model` | `Fast` | 预测模型：`Fast`（ESM1b，快速）或 `Accurate`（ProtT5，高精度） |
| `--device` | `cpu` | 计算设备：`cpu`、`cuda`（GPU）、`mps`（Apple Silicon） |
| `--plot` | `False` | 是否绘制attention图 |
| `--singularity-image` | `~/software/singularity/deeploc2.1_latest.sif` | Singularity镜像路径 |
| `--database-dir` | `~/software/deeploc/database` | 数据库目录 |
| `--singularity-exec` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` | Singularity可执行文件路径 |

### 支持的输入文件格式 | Supported Input Formats

| 格式 | 扩展名 | 说明 |
|------|--------|------|
| FASTA | `.fa`, `.faa`, `.fasta`, `.ffn`, `.fna` | 标准蛋白质序列文件 |

## 预测模型对比 | Model Comparison

| 特性 | Fast (ESM1b) | Accurate (ProtT5) |
|------|-------------|-------------------|
| 底层模型 | ESM1b | ProtT5 |
| 预测速度 | 快 | 慢 |
| 预测精度 | 标准 | 更高 |
| 适用场景 | 大规模批量预测 | 少量高精度预测 |

## 输出文件结构 | Output Structure

```
output_dir/
├── deeploc_prediction.log              # 运行日志
├── {input}_deeploc2.tsv                # DeepLoc原始预测结果
├── {input}_deeploc2_plot.png           # Attention可视化图（--plot时生成）
├── results_*.csv                       # DeepLoc原始CSV结果
├── deeploc_results_readable.tsv        # 可读TSV格式（中文表头）⭐
└── deeploc_results.xlsx                # Excel格式（带样式）⭐
```

⭐ = 本模块新增的格式化输出

## 结果格式说明 | Result Format Details

### 可读TSV文件 (deeploc_results_readable.tsv)

包含以下列（中英文对照表头）：

| 列名 | 说明 |
|------|------|
| 蛋白质ID | 输入FASTA中的序列ID |
| 亚细胞定位 | 预测的亚细胞定位（中文翻译） |
| 信号序列 | 预测的信号序列类型（中文翻译） |
| 膜蛋白类型 | 膜蛋白分类（中文翻译） |
| 细胞质 ~ 过氧化物酶体 | 10种定位的概率（百分比格式） |
| 外周膜 ~ 可溶性 | 4种膜蛋白类型的概率（百分比格式） |

### 亚细胞定位类型 | Subcellular Localization Types

| 英文 | 中文 |
|------|------|
| Cytoplasm | 细胞质 |
| Nucleus | 细胞核 |
| Extracellular | 细胞外 |
| Cell membrane | 细胞膜 |
| Mitochondrion | 线粒体 |
| Plastid | 质体 |
| Endoplasmic reticulum | 内质网 |
| Lysosome/Vacuole | 溶酶体/液泡 |
| Golgi apparatus | 高尔基体 |
| Peroxisome | 过氧化物酶体 |

### 信号序列类型 | Signal Types

| 英文 | 中文 |
|------|------|
| Signal peptide | 信号肽 |
| Transmembrane domain | 跨膜结构域 |
| Chloroplast transit peptide | 叶绿体转运肽 |
| Nuclear localization signal | 核定位信号 |
| Nuclear export signal | 核输出信号 |
| Mitochondrial transit peptide | 线粒体转运肽 |
| Peroxisomal targeting signal | 过氧化物酶体靶向信号 |

### 膜蛋白类型 | Membrane Protein Types

| 英文 | 中文 |
|------|------|
| Peripheral | 外周膜 |
| Transmembrane | 跨膜 |
| Lipid anchor | 脂质锚定 |
| Soluble | 可溶性 |

## 工作流程 | Workflow

```
1. 参数解析与验证
   ↓
2. 配置初始化
   - 验证输入文件格式
   - 检查Singularity镜像、数据库、可执行文件
   - 创建输出目录
   ↓
3. 日志初始化
   - 同时输出到文件和终端
   ↓
4. 构建执行命令
   - 自动检测Singularity的Conda环境
   - 构建singularity exec命令
   - 挂载输入目录、输出目录、数据库目录
   ↓
5. 运行DeepLoc 2.1
   - 在Singularity容器中执行deeploc2
   ↓
6. 检查输出文件
   - 验证预测结果文件
   - 验证attention图（如果启用）
   ↓
7. 格式化结果 ⭐
   - 生成可读TSV（中文表头，概率转百分比）
   - 生成Excel文件（带样式和列宽设置）
   ↓
8. 输出最终摘要
```

⭐ = 本模块新增的功能

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Singularity** (>= 3.8.7)
  - 推荐通过Conda安装: `conda install -c conda-forge singularity`
- **DeepLoc 2.1 Singularity镜像**
  - 默认路径: `~/software/singularity/deeploc2.1_latest.sif`
- **DeepLoc数据库**
  - 默认路径: `~/software/deeploc/database`
- **openpyxl**（可选，用于生成Excel）
  - 安装: `pip install openpyxl`

### Python依赖 | Python Dependencies

- Python >= 3.7
- 标准库: `argparse`, `logging`, `os`, `pathlib`, `subprocess`, `sys`, `re`, `shutil`, `csv`, `glob`, `dataclasses`
- 可选: `openpyxl`（Excel输出）

## 使用示例 | Usage Examples

### 示例1: 标准预测

```bash
# 使用默认Fast模型预测蛋白质亚细胞定位
biopytools deeploc \
    -i /path/to/proteins.faa \
    -o /path/to/output

# 查看结果
cat /path/to/output/deeploc_results_readable.tsv
```

### 示例2: 高精度预测 + GPU

```bash
# 使用Accurate模型和GPU加速（适合小规模高精度需求）
biopytools deeploc \
    -i proteins.faa \
    -o output_accurate \
    --model Accurate \
    --device cuda

# 查看Excel结果
libreoffice output_accurate/deeploc_results.xlsx
```

### 示例3: 带Attention可视化

```bash
# 生成attention权重图，用于分析模型关注的序列区域
biopytools deeploc \
    -i proteins.faa \
    -o output_with_plot \
    --plot

# 查看attention图
# 输出: output_with_plot/{input}_deeploc2_plot.png
```

### 示例4: 自定义Singularity路径

```bash
# 使用自定义的Singularity安装路径和镜像
biopytools deeploc \
    -i proteins.faa \
    -o output \
    --singularity-image /custom/path/deeploc2.1.sif \
    --database-dir /custom/path/database \
    --singularity-exec /custom/path/bin/singularity
```

## 故障排查 | Troubleshooting

### 1. Singularity镜像不存在

**症状**: `Singularity镜像不存在|Singularity image not found: ~/software/singularity/deeploc2.1_latest.sif`

**解决**:
- 检查镜像路径是否正确
- 下载DeepLoc 2.1 Singularity镜像并放置到正确路径
- 或使用`--singularity-image`指定正确路径

### 2. 数据库目录不存在

**症状**: `数据库目录不存在|Database directory not found: ~/software/deeploc/database`

**解决**:
- 检查数据库目录路径
- 使用`--database-dir`指定正确路径

### 3. Singularity可执行文件不存在

**症状**: `Singularity可执行文件不存在|Singularity executable not found`

**解决**:
```bash
# 安装Singularity
conda install -c conda-forge singularity

# 或指定正确的Singularity路径
biopytools deeploc -i proteins.faa -o output \
    --singularity-exec /path/to/singularity
```

### 4. Excel文件未生成

**症状**: 日志显示`未安装openpyxl，跳过Excel生成`

**解决**:
```bash
pip install openpyxl
```

### 5. Conda环境未自动检测

**症状**: Singularity命令直接调用失败

**解决**:
- 确保Singularity安装在Conda环境中（路径包含`/envs/`）
- 或手动使用`conda run`包装命令

## 性能参考 | Performance Reference

| 模型 | 设备 | 适用规模 |
|------|------|----------|
| Fast (ESM1b) | CPU | 大规模批量预测（>1000条序列） |
| Fast (ESM1b) | CUDA | 中大规模批量预测 |
| Accurate (ProtT5) | CUDA | 少量高精度预测（<500条序列） |
| Accurate (ProtT5) | CPU | 少量序列，不追求速度 |

## 与原始DeepLoc 2.1的对比 | Comparison with Original DeepLoc 2.1

| 特性 | 原始DeepLoc 2.1 | biopytools deeploc |
|------|-----------------|-------------------|
| 运行方式 | 需要自行配置Python环境 | Singularity容器化，开箱即用 |
| Conda环境 | 需要手动激活 | 自动检测并包装 |
| 结果格式 | 仅英文TSV/CSV | 中英文对照TSV + Excel |
| 概率显示 | 小数格式 | 百分比格式（可读TSV） |
| Excel输出 | 无 | 带样式的Excel文件 |
| 日志记录 | 无 | 完整的文件+终端日志 |

## 版本历史 | Version History

- **v1.0.0**: 初始版本
  - Fast/Accurate双模型支持
  - Singularity容器化运行
  - 中文结果翻译（定位、信号、膜蛋白）
  - TSV和Excel格式化输出
  - Attention可视化支持
  - 自动Conda环境检测

## 许可证 | License

本模块遵循biopytools项目的许可证。DeepLoc 2.1请参考其原始许可证。
