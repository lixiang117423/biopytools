# Merge-DeepBSA 合并模块

**Windows版DeepBSA结果合并工具 | Windows DeepBSA Results Merger**

## 功能概述 | Overview

Merge-DeepBSA 模块用于合并 Windows 版 DeepBSA 的运行结果，将分散在各方法目录中的 QTL 结果和绘图数据整合为统一格式，便于后续分析和可视化。

## 主要特性 | Key Features

- **自动方法检测**: 扫描输入目录，自动识别 CSV、values.txt、npy 文件对应的分析方法
- **QTL结果合并**: 将所有方法的 QTL 结果合并为 `merged_results.xlsx`，自动标注来源方法
- **绘图数据提取**: 生成 `plot_data_for_R.csv`，包含平滑值和阈值，可直接用于 R 绘图
- **智能数据源**: 优先从 npy 文件读取（与原始软件完全一致），npy 不存在时 fallback 到 values.txt 重新计算
- **方法筛选**: 支持指定合并特定方法，默认合并所有检测到的方法

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 合并所有方法的结果
biopytools merge-deepbsa -i ./jicaiBSA_Visualize_Results -o ./merged_results

# 只合并指定方法
biopytools merge-deepbsa -i ./results -o ./merged -m DL,K,ED4,SNP

# 自定义平滑参数（仅 fallback 模式生效）
biopytools merge-deepbsa -i ./results -o ./merged --smooth-func Tri-kernel --smooth-frac 0.15
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-dir` | Windows版DeepBSA输出目录 | `-i ./jicaiBSA_Visualize_Results` |
| `-o, --output-dir` | 合并结果输出目录 | `-o ./merged_results` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --methods` | 全部 | 要合并的方法，逗号分隔（DL,K,ED4,SNP,SmoothG,SmoothLOD,Ridit） |
| `--smooth-func` | `LOWESS` | 平滑函数，仅 fallback 模式使用（LOWESS/Tri-kernel） |
| `--smooth-frac` | `0.1` | 平滑窗口比例，仅 fallback 模式使用 |

## 支持的方法 | Supported Methods

| 方法 | 描述 |
|------|------|
| **DL** | Deep Learning 基于深度学习 |
| **K** | K-Method 基于K值 |
| **ED4** | Euclidean Distance 4 基于欧氏距离 |
| **SNP** | SNP-Index SNP指数法（输入/输出中显示为 ΔSNP） |
| **SmoothG** | Smooth G-value 平滑G值法 |
| **SmoothLOD** | Smooth LOD 平滑LOD值法 |
| **Ridit** | Ridit Analysis Ridit分析法 |

## 输入文件格式 | Input File Formats

输入目录应为 Windows 版 DeepBSA 的输出目录，包含以下文件：

### CSV 结果文件

文件名格式：`0-{method}-Tri-kernel-smooth-0.1-{value}.csv`

```
0-DL-Tri-kernel-smooth-0.1-auto.csv
0-ED4-Tri-kernel-smooth-0.1-auto.csv
0-ΔSNP-Tri-kernel-smooth-0.1-auto.csv
```

### values.txt 原始数据文件

文件名格式：`{method} values.txt`

```
DL values.txt
ED4 values.txt
ΔSNP values.txt
```

文件内容为 TSV 格式（染色体、位置、值）：
```
Chr01	100000	0.0234
Chr01	100500	0.0156
...
```

### npy 绘图数据文件（可选，优先使用）

文件名格式：`all_data_for_plot_{method}.npy`

```
all_data_for_plot_DL.npy
all_data_for_plot_ED4.npy
all_data_for_plot_ΔSNP.npy
```

## 输出文件 | Output Files

```
output_dir/
├── merged_results.xlsx          # 各方法QTL合并表
├── plot_data_for_R.csv         # 绘图数据表（含平滑值和阈值）
└── merge_deepbsa.log           # 运行日志
```

### merged_results.xlsx

合并所有方法的 QTL 结果，包含以下列：

| 列名 | 描述 |
|------|------|
| **Method** | 来源方法（SNP方法显示为 ΔSNP） |
| **QTL** | QTL编号 |
| **Chr** | 染色体 |
| **Left** | 左侧位置 |
| **Peak** | 峰值位置 |
| **Right** | 右侧位置 |
| **Value** | QTL值 |
| **Source_File** | 原始CSV文件名 |

注意：仅保留 Value 列不为 "-" 的有效 QTL 结果。

### plot_data_for_R.csv

包含各方法的绘图数据，可直接用于 R 语言绑图：

| 列名 | 描述 |
|------|------|
| **Method** | 方法名 |
| **Chromosome** | 染色体 |
| **Position** | 位置 |
| **Value** | 原始值 |
| **Smooth_Value** | 平滑值 |
| **Threshold** | 阈值 |

## 工作流程 | Workflow

```
1. 扫描输入目录
   - 检测 CSV、values.txt、npy 文件
   - 识别可用的分析方法
   ↓
2. 合并QTL结果
   - 遍历所有CSV文件
   - 提取方法名，过滤无效行
   - 添加 Method 和 Source_File 列
   - 合并保存为 merged_results.xlsx
   ↓
3. 提取绘图数据
   - 优先从 npy 文件读取（与原始软件一致）
   - npy 不存在时 fallback 到 values.txt 重新计算
   - 保存为 plot_data_for_R.csv
   ↓
4. 输出统计信息
   - 各方法QTL数量
   - 总计统计
```

## 数据源优先级 | Data Source Priority

绘图数据的提取按以下优先级：

1. **npy 文件**（优先）：从 `all_data_for_plot_{method}.npy` 读取，与 DeepBSA 原始软件的绘图数据完全一致
2. **values.txt**（fallback）：从 `{method} values.txt` 读取原始数据，使用 `PlotDataCalculator` 重新计算平滑值和阈值

当 npy 文件存在时直接使用，避免重新计算带来的微小数值差异。

## 故障排查 | Troubleshooting

### 输入目录中未找到数据文件

**错误信息**: `输入目录中未找到CSV、values.txt或npy文件`

**解决**: 确保输入目录是 Windows 版 DeepBSA 的输出目录，包含 `*.csv`、`* values.txt` 或 `all_data_for_plot_*.npy` 文件。

### CSV文件被跳过

**日志信息**: `跳过非标准CSV文件`

**解决**: 确保 CSV 文件名符合 `0-{method}-Tri-kernel-smooth-{frac}-{value}.csv` 格式。

### 平滑值与原始软件不一致

**解决**: 确保输入目录中包含 npy 文件（`all_data_for_plot_*.npy`）。npy 文件存储了原始软件计算的精确绘图数据，使用 values.txt 重新计算可能存在微小数值差异。

## 版本历史 | Version History

- **v1.0**: 初始版本
  - QTL结果合并（merged_results.xlsx）
  - 绘图数据提取（plot_data_for_R.csv）
  - npy 文件优先读取
  - values.txt fallback 计算
