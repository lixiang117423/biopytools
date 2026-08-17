# DeepLoc - 蛋白质亚细胞定位预测 | DeepLoc Protein Subcellular Localization Prediction

一句话理解：**给每条蛋白质序列猜出它在细胞里的「住址」**——是待在细胞核里、挂在细胞膜上、还是被分泌到细胞外，顺便告诉你它是跨膜蛋白还是可溶性蛋白。输入一个蛋白质 FASTA 文件，输出一张带定位结果和置信概率的表格。

## 功能概述 | Overview

- 基于 DeepLoc 2.1 深度学习模型，预测蛋白质的亚细胞定位与膜关联性
- 覆盖 10 类定位：细胞质、细胞核、细胞外、细胞膜、线粒体、质体、内质网、溶酶体/液泡、高尔基体、过氧化物酶体
- 同时预测信号类型（信号肽、跨膜结构域、叶绿体转运肽、核定位信号等）与膜蛋白类型（外周膜/跨膜/脂质锚定/可溶性）
- 两种模型可选：Fast（ESM1b，快）与 Accurate（ProtT5，高精度）
- 通过 Singularity 容器运行，免去手动装依赖；输出原始 TSV 之外，自动转出中文易读 TSV 与 Excel
- 断点续传：主输出已存在时跳过昂贵的容器预测，只重做轻量格式化

## 快速开始 | Quick Start

```bash
biopytools deeploc -i proteins.faa -o output_dir
```

最小输入：一个蛋白质 FASTA 文件（.fa / .faa / .fasta / .ffn / .fna）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 亚细胞定位 | 蛋白质在细胞里的「工作地址」：像一座城市里，快递要送到小区、写字楼还是工厂 |
| 信号肽 | 蛋白质 N 端的一小段「地址标签」，告诉细胞该把它送到哪里，送到后常被剪掉 |
| 跨膜结构域 | 蛋白质「横穿」细胞膜的一段，像一根钉子钉穿木板，把蛋白固定在膜上 |
| 膜蛋白类型 | 蛋白和细胞膜的「关系」：可溶性=自由漂浮；跨膜=钉穿膜；外周膜=贴在膜表面；脂质锚定=用一根「油绳」拴在膜上 |
| ESM1b / ProtT5 | 两个不同的深度学习模型：前者轻量快速，后者更大更准（对应 Fast/Accurate） |
| Singularity 容器 | 把软件连同所有依赖打包成的「便携实验室」，一条命令即可在隔离环境里运行 |
| device | 用什么硬件算：cpu（普通 CPU）、cuda（NVIDIA GPU）、mps（苹果 Apple Silicon GPU） |

## 输入 | Input

蛋白质 FASTA 文件，每条序列一个条目：

```text
>protein_001
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
```

- 标准 FASTA 格式；序列必须是氨基酸（蛋白质），不接受核酸
- 扩展名需为 .fa / .faa / .fasta / .ffn / .fna 之一（否则校验报错）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这两个必须给：要预测哪些蛋白（输入文件），结果写到哪（输出目录）。没有默认值，缺一个就跑不起来。

### 模型与设备 | Model & device

**通俗理解|In plain words:** 决定「用哪个脑子、用哪种硬件」。Fast 用轻量模型，速度快、适合几万条序列的大蛋白组初筛；Accurate 用更大的模型，更准但更慢，适合重点关注的少数蛋白。device 默认 cpu，没有 GPU 就不用动；有 NVIDIA 卡改 cuda、苹果 M 系列芯片可改 mps 能明显加速。**一般先用默认 Fast+cpu。**

### 运行环境路径 | Runtime paths

**通俗理解|In plain words:** 这三个是「工具装在哪」的路径，指向 Singularity 容器镜像、DeepLoc 数据库目录和 singularity 可执行文件。部署时管理员已配好默认值，**普通用户一般不用动**；只有当你的镜像/数据库装在不同位置时才需要指定。

### 绘图选项 | Plotting

**通俗理解|In plain words:** `--plot` 让程序额外画注意力图（attention），用于看模型「重点关注序列的哪一段」。绝大多数场景不需要，**默认关着即可**。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先检查结果是否已算过（有就直接复用），没有就用 Singularity 容器跑 DeepLoc 2.1，再把原始 CSV 转成中文易读的 TSV 和 Excel。

```text
蛋白质 FASTA
    │
    ▼
(断点续传检查) 主输出 {输入名}_deeploc2.tsv 是否已存在?
    ├─ 是 → 跳过容器预测,直接走格式化
    └─ 否 → singularity exec deeploc2 (挂载输入/输出/数据库目录)
              │
              ▼
         生成 {输入名}_deeploc2.tsv + results_*.csv
    │
    ▼
格式化: results_*.csv → deeploc_results_readable.tsv + deeploc_results.xlsx
    │
    ▼
记录软件版本到 00_pipeline_info/
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml             # 软件/镜像/数据库版本与参数记录
├── {输入名}_deeploc2.tsv                 # DeepLoc 原始预测结果(主输出,续传判据)
├── results_*.csv                         # DeepLoc 容器输出的原始 CSV
├── deeploc_results_readable.tsv          # 中文易读 TSV(概率转百分比)
├── deeploc_results.xlsx                  # 中文 Excel(带样式)
├── {输入名}_deeploc2_plot.png            # 注意力图(仅 --plot 时)
└── 99_logs/
    └── deeploc.log                       # 运行日志
```

- `{输入名}` 为输入 FASTA 去掉扩展名后的文件名（如 proteins.faa → proteins）
- `deeploc_results_readable.tsv` 是给用户看的主结果：前 4 列是蛋白 ID、定位、信号、膜类型（已翻译成中文），后 14 列是各定位/膜类型的概率（已转成百分比）

## 结果解读 | Interpreting Results

### 1. 可读 TSV（`deeploc_results_readable.tsv`）

**通俗理解|In plain words:** 每行一个蛋白，前几列是「结论」，后 14 列是「各选项的得票率」。看结论列即可；想较真再看概率列。

- `蛋白质ID`：序列名；`亚细胞定位`：预测最可能的定位（可为多个，用 | 分隔，如「细胞核|细胞质」）
- `信号序列`：预测到的信号类型（信号肽、跨膜结构域、线粒体转运肽等），没有则为 `-`
- `膜蛋白类型`：外周膜 / 跨膜 / 脂质锚定 / 可溶性
- 概率列（细胞质、细胞核、细胞外、细胞膜……共 10 个定位 + 4 个膜类型）：每个取值 0%–100%，**最高的一类即预测定位**

### 2. 好坏判据 | Judgment

- **概率越高越可信**：定位概率 >80% 的结果基本可靠；50%–80% 谨慎看待；多个定位概率接近时说明该蛋白可能是「双定位」或模型不确定
- 模型给出「信号肽 + 细胞外」的蛋白，通常是分泌蛋白的典型特征
- 预测「跨膜 + 细胞膜」的蛋白，是膜蛋白的典型特征；两者一致可互相印证

## 参数选择建议 | Parameter Guidance

- **`--model`**：蛋白组级大扫描用 Fast；对少数关键蛋白要精准结论时用 Accurate（如验证某个候选效应蛋白的定位）
- **`--device`**：默认 cpu；有 NVIDIA GPU 用 cuda、Apple Silicon 用 mps 可显著加速，尤其 Accurate 模型
- **`--plot`**：仅在需要看模型关注区域时开启；日常批量预测关闭以省时省空间
- **其余路径参数**：部署好后就别动，换环境时再指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件(蛋白质序列)｜Input FASTA file (protein sequences) |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--model` | `Fast` | Fast/Accurate | 预测模型｜Prediction model |
| `--device` | `cpu` | cpu/cuda/mps | 计算设备｜Compute device |
| `--singularity-image` | `~/software/singularity/deeploc2.1_latest.sif` |  | Singularity镜像路径｜Singularity image path |
| `--database-dir` | `~/software/deeploc/database` | Path | 数据库目录｜Database directory |
| `--singularity-exec` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity可执行文件路径｜Singularity executable path |
| `--plot` | — |  | 绘制attention图｜Plot attention values |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件(蛋白质序列)｜Input FASTA file (protein sequences) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--model` | `Fast` | Fast/Accurate | 预测模型｜Prediction model (default: Fast) Fast: ESM1b模型(快速)｜Fast: ESM1b model (fast) Accurate: ProtT5模型(高精度)｜Accurate: ProtT5 model (accurate) |
| `--device` | `cpu` | cpu/cuda/mps | 计算设备｜Compute device (default: cpu) |
| `--singularity-image` | `~/software/singularity/deeploc2.1_latest.sif` |  | Singularity镜像路径｜Singularity image path |
| `--database-dir` | `~/software/deeploc/database` |  | 数据库目录｜Database directory |
| `--singularity-exec` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity可执行文件路径｜Singularity executable path |
| `--plot` | — | store_true | 绘制attention图｜Plot attention values |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Singularity（可执行文件，默认 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`，通过 conda 环境自动包装调用）
- DeepLoc 2.1 Singularity 镜像（默认 `~/software/singularity/deeploc2.1_latest.sif`）
- DeepLoc 数据库目录（默认 `~/software/deeploc/database`）
- Python 3（openpyxl 用于生成 Excel，缺失时自动跳过 Excel）

## 常见问题 | FAQ

**Q1：换模型重跑，为什么结果没变？**
断点续传按主输出 `{输入名}_deeploc2.tsv` 是否存在判断。换 `--model`/`--device` 重跑前，先删掉输出目录里的 `{输入名}_deeploc2.tsv`（或换一个输出目录），否则会复用旧结果只重做格式化。

**Q2：报「Singularity 镜像不存在 / 数据库目录不存在」怎么办？**
这是部署问题：默认路径里没有镜像或数据库。用 `--singularity-image`、`--database-dir` 指定实际位置，或联系管理员安装。

**Q3：能预测核酸序列吗？**
不能。本工具只接受蛋白质（氨基酸）FASTA；扩展名不对会被校验拦下。

**Q4：Excel 没生成？**
生成 Excel 需要 openpyxl 包。缺失时程序会跳过 Excel 但保留 TSV；安装 `pip install openpyxl` 后重跑即可（此时已有主输出，会走续传，直接重新格式化）。