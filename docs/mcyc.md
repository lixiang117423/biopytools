# MCyc - 甲烷循环基因丰度分析 | Methane Cycle Gene Abundance Analysis

一句话理解：**用 MCycDB 数据库给宏基因组/宏转录组 reads 做比对，算出「产甲烷、甲烷氧化等甲烷循环相关基因」在各样本里的丰度**，产出原始计数、TPM、CLR 三张矩阵。

## 功能概述 | Overview { #overview }

- 输入一个样本清单文件（样本名 + FASTQ 路径），自动准备数据并跑 MCycDB Perl 流程
- 用 Diamond 比对 MCycDB 蛋白库，统计甲烷循环相关基因的 reads 数
- 输出三张矩阵：RawCounts（原始计数）、TPM（每百万归一化）、CLR（对数比变换）
- Diamond 索引缺失时自动构建；已有中间结果可跳过比对（断点续传）
- 运行结束自动清理临时文件与数据库软链接（可 --keep-temp 保留）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools mcyc samples.txt
```

最小输入：一个两列样本清单文件。输出默认写到当前目录（可用 `-o` 指定）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 甲烷循环 | 产甲烷菌「造」甲烷、甲烷氧化菌「吃」甲烷的碳循环过程 |
| MCycDB | 专门收录甲烷循环相关基因的参考数据库 |
| 丰度矩阵 | 「基因 × 样本」的数值表，行列交叉处是该基因在某样本里的量 |
| RawCounts | 比对上去的原始 reads 数，直接反映「测了多少」 |
| TPM(此处实为CPM) | 每百万 reads 归一化，把不同测序深度的样本拉到同一尺度再比较 |
| CLR | 中心对数比变换，把相对丰度变成便于统计的对数尺度 |
| Diamond | 一种快速的蛋白序列比对工具，把 reads 比对到蛋白库 |

## 输入 | Input { #input }

- **样本清单文件**（位置参数）：每行一个样本，`样本名 + 文件路径` 两列（制表符或空格分隔），以 `#` 开头的行被忽略。
- **双端数据**：R1、R2 用分号 `;` 连接，程序会把两端合并（cat）成一个文件再比对。
- **单端数据**：只写一个路径即可。

示例 `samples.txt`：

```text
# 样本名  文件路径(双端用分号连接)
S1  /data/S1_R1.fastq.gz
S2  /data/S2_R1.fastq.gz;/data/S2_R2.fastq.gz
```

## 参数说明 | Parameters { #parameters }

### 输入与输出 | Input & output

**通俗理解|In plain words:** 位置参数是样本清单文件；`-o` 输出目录默认是**当前目录**（不是某个子目录），三张矩阵直接写在那里。`-t` 线程数默认 12（CLI）用于 seqkit 统计等步骤。

### 数据库 | Database

**通俗理解|In plain words:** `--mcyc-base` 指定 MCycDB 目录（默认 `~/software/MCycDB/MCycDB-main`），程序会从中找 `MCycDB_2021.faa`（蛋白库）、`id2gene.map`（基因映射）、`MCycDB_2021.dmnd`（Diamond 索引）。**数据库装好后一般不用动**；Diamond 索引不存在时会自动构建。

### 跳过与保留 | Skip & keep

**通俗理解|In plain words:** `--skip-diamond` 跳过比对，直接复用已有的中间结果文件 `MCyc_Raw_Temp.txt`（适合中断后只重算矩阵）；`--keep-temp` 保留临时文件与数据库软链接以便排查。**日常用默认即可。**

## 分析流程 | Pipeline { #pipeline }

```text
检查依赖(diamond/seqkit/perl)
    │
    ▼
准备数据库(必要时 diamond makedb + 建软链接)
    │
    ▼
准备数据(清单解析 → staging 目录软链接/合并双端)
    │
    ▼
统计 reads 数(seqkit stats → sample_info.txt)
    │
    ▼
Diamond 比对(perl MCycDB_FunctionProfiler.PL，已有结果则跳过)
    │
    ▼
计算三张矩阵(RawCounts / TPM / CLR)
    │
    ▼
清理临时文件与软链接
```

## 输出 | Output { #output }

```text
<输出目录，默认当前目录>/
├── Matrix_1_RawCounts.txt      # 原始计数矩阵
├── Matrix_2_TPM.txt            # TPM(实为每百万 CPM)矩阵
├── Matrix_3_CLR.txt            # CLR 变换矩阵
└── mcyc_analysis.log           # 运行日志
```

中间临时产物（`mcyc_staging_fastq/`、`sample_info.txt`、`MCyc_Raw_Temp.txt`、数据库软链接 `MCycDB_2021.*`）默认在运行结束清理，加 `--keep-temp` 保留。

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. Matrix_1_RawCounts.txt（原始计数）

**通俗理解|In plain words:** 行是甲烷循环相关基因，列是样本，数值是比对到的原始 reads 数。它受测序深度影响大，样本间直接比较不公平。

### 2. Matrix_2_TPM.txt（归一化丰度）

**通俗理解|In plain words:** 按「每 100 万 reads」折算后的丰度（实现为 CPM：reads / 该样本总 reads × 1e6），消除了测序深度差异。**跨样本比较丰度时优先用它。** 注意：这里的 TPM 不是转录组意义上的长度归一化 TPM，而是每百万计数（CPM）。

### 3. Matrix_3_CLR.txt（CLR 变换）

**通俗理解|In plain words:** 对每个样本做中心对数比变换（先加 1 伪计数，再取 log 并减几何均值对数），把组成型数据变成适合做 PCA、聚类、相关性分析的对数尺度。

### 4. 好坏判据

- **比对率低/矩阵多为 0**：检查 reads 是核酸还是蛋白（脚本按 nucl 模式跑），或样本本就不含甲烷循环相关类群。
- **样本间总数差异悬殊**：RawCounts 会明显倾斜，应改用 TPM/CLR。
- **CLR 出现极端值**：通常是某样本 reads 极少导致的伪计数主导，考虑剔除极低深度样本。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **常规流程**：`biopytools mcyc samples.txt -o results/` 即可。
- **跨样本比较**：用 TPM 或 CLR 矩阵，不要用 RawCounts。
- **中断后只想重算矩阵**：加 `--skip-diamond`（复用已有 MCyc_Raw_Temp.txt）。
- **排查问题**：加 `--keep-temp` 保留中间文件。
- **数据库不在默认位置**：用 `--mcyc-base /path/to/MCycDB-main`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--output, -o` | — | Path | 输出目录路径｜Output directory path |
| `--threads, -t` | `12` | int | 线程数｜Thread count |
| `--mcyc-base` | — | Path | MCycDB基础路径｜MCycDB base path |
| `--skip-diamond` | — |  | 跳过Diamond比对｜Skip Diamond alignment |
| `--keep-temp` | — |  | 保留临时文件｜Keep temporary files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-o, --output` | — |  | 输出目录路径｜Output directory path (default: current directory) |
| `-t, --threads` | `4` | int | 线程数｜Thread count (default: 4) |
| `--mcyc-base` | — |  | MCycDB基础路径｜MCycDB base path |
| `--skip-diamond` | — | store_true | 跳过Diamond比对｜Skip Diamond alignment |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep temporary files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- diamond（蛋白比对，需在 PATH 中）
- seqkit（reads 统计，需在 PATH 中）
- perl（运行 MCycDB_FunctionProfiler.PL，需在 PATH 中）
- MCycDB 数据库（默认 `~/software/MCycDB/MCycDB-main`，含 MCycDB_2021.faa / id2gene.map / MCycDB_2021.dmnd）
- Python 包：pandas、numpy（矩阵计算）

> 说明：本模块的 diamond/seqkit/perl 通过 `shutil.which` 在 PATH 中查找并直接调用（shell 执行），**不走 conda run**。

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
部分支持。若 `MCyc_Raw_Temp.txt` 已存在且非空，会跳过 Diamond 比对直接进入矩阵计算；样本的 staging 文件已存在也会跳过重新准备。中断后用 `--skip-diamond` 可显式复用已有比对结果。

**Q2：为什么输出矩阵出现在当前目录而不是子目录？**
`-o` 默认是当前目录，且工作目录（临时文件、软链接）也是当前目录。建议在专用目录下运行，或显式 `-o` 到目标目录，避免污染。

**Q3：线程数到底默认多少？**
CLI（`biopytools mcyc`）默认 12；直接 `python -m biopytools.mcyc` 时 argparse 默认 4。两者不一致，如在意请显式 `-t`。

**Q4：双端数据怎么传？**
清单文件里 R1、R2 用分号连接（`R1.fq.gz;R2.fq.gz`），程序会把两端 cat 合并成一个文件再比对。

**Q5：报「MCycDB 文件不存在」？**
检查 `--mcyc-base` 指向的目录是否含 `MCycDB_2021.faa` 与 `id2gene.map`；缺失说明数据库未放对位置。

**Q6：TPM 和 CPM 有什么区别？**
本模块的 TPM 实际按 CPM（每百万计数）计算，未做基因长度归一化；名称沿用 MCycDB 习惯。如需真 TPM 请自行换算。