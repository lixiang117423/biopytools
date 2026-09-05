# PredGPI - GPI锚定蛋白预测 | PredGPI GPI-anchored Protein Prediction

一句话理解：**判断每条蛋白是不是「GPI 锚定蛋白」——即会不会用一根脂质尾巴把自己挂在细胞膜外表面**。输入一个蛋白质 FASTA，输出一张带分类和置信度的表格，并标出预测的切割位点。

## 功能概述 | Overview

- 基于 PredGPI（HMM + SVM 联合模型）预测 GPI 锚定蛋白
- 输出每条蛋白是否 GPI 锚定、预测切割位点、假阳性率(FPR)、HMM/SVM 打分与分类标签
- 提供保守模型（`--conservative`）可选，换用更严格的高特异性模型
- 纯 Python 计算，无外部命令行调用；自动处理非标准氨基酸（U→C、Z/B/X→A）与过短序列（≤40 aa 直接判非 GPI）
- 断点续传：主输出 TSV 已存在则跳过

## 快速开始 | Quick Start

```bash
biopytools predgpi -i proteins.fa -o output_dir/
```

最小输入：一个蛋白质 FASTA 文件 + 输出目录。输出前缀默认取输入文件名。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GPI 锚定蛋白 | 一种「用脂质尾巴拴在细胞膜外表面」的蛋白：蛋白主体在膜外，靠一条糖脂尾巴固定 |
| 切割位点 | GPI 锚合成时在蛋白 C 端切开、接上脂质尾巴的位置 |
| FPR(假阳性率) | 该预测「误报为 GPI 蛋白」的估计概率，越小越可信 |
| HMM / SVM | 两种机器学习模型：HMM 看序列模式，SVM 做分类，两者打分结合判定 |

## 输入 | Input

蛋白质 FASTA 文件（支持 .fa / .faa / .fasta / .fna / .ffn，含 .gz 压缩）：

```text
>protein_01
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLS...
```

- 序列必须是氨基酸（蛋白质）；非标准氨基酸会被自动替换为相近的标准氨基酸
- 长度 ≤40 的序列过短，会被直接判为非 GPI（FPR 记为 1）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：输入蛋白文件、输出目录。没有默认值。

### 安装目录与模型 | Home & model

**通俗理解|In plain words:** `--predgpi-home` 指向 PredGPI 安装目录（含模型文件 GPIDAT），部署时已配好，**普通用户不用动**。`--conservative` 换用更严格的保守模型，**牺牲灵敏度换特异性**——只在需要「尽量少误报」时开启。

### 输出前缀 | Prefix

**通俗理解|In plain words:** 决定输出文件名，默认取输入文件名，**一般不用动**。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 加载 HMM/SVM 模型后逐条预测，把原始打分翻译成人能看懂的分类，写成 TSV。

```text
蛋白质 FASTA
    │
    ▼
(断点续传) {prefix}_predgpi.tsv 已存在? → 是则跳过
    │
    ▼
加载 HMM + SVM 模型(GPIDAT)
    │
    ▼
逐条预测: 替换非标准AA → 短序列直接判非GPI → 其余跑 predGpipe
    │
    ▼
按 FPR 分类 → 写 {prefix}_predgpi.tsv + software_versions.yml
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml     # 软件/模型版本与参数记录
├── {prefix}_predgpi.tsv          # 主结果:预测汇总表
└── 99_logs/
    └── predgpi.log               # 运行日志
```

- 主结果 `{prefix}_predgpi.tsv` 列：`ID / Length / GPI_Anchored / Cleavage_Site / FPR / HMM_LogProb / SVM_Score / Probability / Classification`

## 结果解读 | Interpreting Results

### 1. 汇总表（`{prefix}_predgpi.tsv`）

**通俗理解|In plain words:** 每行一个蛋白，先看 `GPI_Anchored`（是不是），再看 `Classification`（多可信）。

- `GPI_Anchored`：True/False，判定标准是 FPR ≤ 0.01
- `Cleavage_Site`：预测的 GPI 锚切割位点（蛋白序列坐标）；非 GPI 为 `-`
- `FPR`：假阳性率估值，越小越可信；`FPR = 1` 表示短序列或预测失败（被标记非 GPI）
- `Probability` 与 `Classification`：由 FPR 映射的分类标签

### 2. 分类阈值 | Classification thresholds

**通俗理解|In plain words:** FPR 越低，越敢拍胸脯说它是 GPI 蛋白。

| 分类 | FPR 范围 | 含义 |
|------|----------|------|
| Highly Probable | ≤ 0.0015 | 高度可信的 GPI 锚定蛋白 |
| Probable | ≤ 0.005 | 较可信 |
| Weakly Probable | ≤ 0.01 | 弱可信（仍判为 GPI） |
| Improbable | > 0.01 | 不太可能是 GPI（判为非 GPI） |

### 3. 好坏判据 | Judgment

- **Highly Probable / Probable**：可作为候选 GPI 锚定蛋白；Highly Probable 几乎不会误报
- **Weakly Probable**：谨慎看待，建议用其他方法（如实验验证或其它预测工具）交叉确认
- 若整批蛋白大部分被判非 GPI，属正常——GPI 锚定蛋白在蛋白组里本来就是少数

## 参数选择建议 | Parameter Guidance

- **`--conservative`**：默认关（更灵敏，可能多一点弱阳性）；做「精选候选」要低误报时开启
- **`--predgpi-home`**：仅部署环境不同时才动
- 短序列（≤40 aa）无需特别处理，程序会自动判为非 GPI

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入蛋白质FASTA｜Input protein FASTA |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--predgpi-home` | — |  | predgpi安装目录｜predgpi install directory |
| `--conservative` | — |  | 使用保守omega模型｜Use conservative omega model |
| `--prefix` | — |  | 输出前缀(默认输入文件名)｜Output prefix (default: input filename) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入蛋白质FASTA｜Input protein FASTA |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `--predgpi-home` | `~/software/predgpi` |  | predgpi安装目录｜predgpi install directory |
| `--conservative` | — | store_true | 使用保守omega模型｜Use conservative omega model |
| `--prefix` | — |  | [STR] 输出前缀(默认输入文件名)｜Output prefix (default: input filename) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- PredGPI 安装目录（默认 `~/software/predgpi`，含 GPIDAT 模型目录：PHMM.TOT.ss.mod / PHMM.TOT.ss.mod_CSDGN / MOD）
- Python 环境可导入 predgpi 包（predgpilib、predGpipe、readFasta）
- Python 3 + PyYAML（记录版本）

## 常见问题 | FAQ

**Q1：换 `--conservative` 重跑，为什么结果没变？**
断点续传按 `{prefix}_predgpi.tsv` 是否存在判断。换模型重跑前先删掉旧 TSV（或换输出目录），否则直接复用旧结果。

**Q2：报「模型文件不存在 / GPIDAT 目录不存在」？**
说明 PredGPI 没装好或路径不对。确认 `--predgpi-home` 指向的目录下有 `GPIDAT/` 及其中的模型文件；`--conservative` 还需 `PHMM.TOT.ss.mod_CSDGN`。

**Q3：为什么有些蛋白 FPR=1 被判非 GPI？**
两种情况：序列过短（≤40 aa）无法可靠预测，或 predgpi 内部对个别序列抛异常，此时程序按约定标记为非 GPI（FPR=1），并在日志里告警。

**Q4：PredGPI 的版本信息？**
PredGPI 无 `--version` 命令，版本以发表文献标注（Pierleoni et al., BMC Bioinformatics 2008, 9:392），写入 `00_pipeline_info/software_versions.yml`。