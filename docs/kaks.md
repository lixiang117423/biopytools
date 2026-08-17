# Ka/Ks 计算 | Ka/Ks Calculation

一句话理解：**给两个物种的「同源基因对」算一个比值 Ka/Ks，用来判断基因在进化中是「被严格看管」还是「在变新花样」**。

输入两个物种的 CDS 序列文件 + 一份「谁和谁是一对」的配对表，输出每对基因的 Ka、Ks、Ka/Ks 比值，以及选择压力分类（负选择/中性/正选择）。

## 功能概述 | Overview { #overview }

- 调用 KaKs_Calculator 2.0 批量计算同源基因对的 Ka（非同义替换率）与 Ks（同义替换率）
- 支持 17 种计算方法（默认 GMYN，即 Gamma 修正的 Yang-Nielsen 方法）
- 输入侧做严格质控：字符、N 比例、序列长度、密码子完整性逐一校验，不合格的配对自动剔除
- 输出逐对详细结果（TSV/CSV）+ Excel 汇总 + JSON 统计，并自动给每个基因打上选择压力标签
- 断点续传：无（每次运行都从头计算，临时文件默认自动清理，可用 `--keep-temp` 保留）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools kaks -1 species1.cds.fasta -2 species2.cds.fasta -p pairs.txt -o results/
```

最小输入：两个物种的 CDS FASTA + 一个两列表格（TSV/CSV，列1=物种1序列ID，列2=物种2序列ID）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| CDS | 编码蛋白的 DNA 序列，长度一定是 3 的倍数（3 个碱基编 1 个氨基酸） |
| 同源基因对 | 两个物种里「同一个基因」，由配对表指定谁和谁比 |
| 同义替换(Ks) | 碱基变了但氨基酸没变的替换，像「同义词替换」，一般认为不受选择影响 |
| 非同义替换(Ka) | 碱基变了氨基酸也变的替换，像「改意思」，会被自然选择盯上 |
| Ka/Ks(omega) | 两者比值：<1 负选择，约=1 中性，>1 正选择，是核心结论 |
| 负选择/净化选择 | 「有害改动被清除」，基因很保守，功能重要 |
| 正选择 | 「有益改动被保留」，基因在快速适应新环境 |

## 输入 | Input { #input }

### 两个 FASTA 文件

物种1、物种2 各自的 CDS 序列（核酸，A/T/C/G/N）。程序会对序列做质控：只接受 A/T/C/G/N 字符；N 比例不超过 10%；长度 50-50000 bp；长度必须是 3 的倍数（否则整条被剔除）。

### 配对文件

TSV（制表符）或 CSV（逗号）格式，2 到 3 列：第1列=物种1序列ID，第2列=物种2序列ID，第3列（可选）=配对名（不写则自动编号 pair_0001）。ID 必须与 FASTA 文件里的序列名完全一致。

```text
gene001    gene001
gene002    gene002
gene003    gene003
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 两个 FASTA、一个配对表、一个输出目录，四个缺一不可。

### 计算方法 | Method

**通俗理解|In plain words:** 管「用哪套公式估计 Ka/Ks」。默认 GMYN 是公认的推荐方法（Gamma 修正、更贴近真实替换过程），**一般不用动**。其余 16 种（YN/MYN/NG/LWL/LPB/GY/MS/MA 及各自 Gamma 版）是不同年代、不同假设的方法，只在你要和某篇文献的口径对齐时才换。

### 运行选项 | Runtime

**通俗理解|In plain words:** `--kaks-path` 指定 KaKs_Calculator 程序路径（默认自动找，一般不用动）；`--verbose` 开详细日志（排查问题用）；`--temp-dir` 自定义临时目录；`--keep-temp` 保留中间 AXT 文件（调试用，正常不需要）。

## 分析流程 | Pipeline { #pipeline }

```text
步骤1: 解析配对表 -> 得到需要的序列 ID 集合
   |
   v
步骤2: 用 seqkit 从两个 FASTA 里抽出这些 ID 的子集
   |
   v
步骤3: 序列质控（字符/N比例/长度/密码子完整性），剔除不合格配对
   |
   v
步骤4: 生成 AXT 格式输入 -> KaKs_Calculator 批量计算
   |
   v
步骤5: 解析结果 -> 添加选择压力分类/质量标记 -> 输出 TSV/CSV/XLSX/JSON
   |
   v
步骤6: 清理临时文件（除非 --keep-temp）
```

## 输出 | Output { #output }

```text
output_dir/
|-- kaks_detailed.tsv        # 逐对详细结果（制表符，含 seq1_id/seq2_id/Ka/Ks/Ka/Ks/分类）
|-- kaks_detailed.csv        # 逐对详细结果（逗号分隔，方便 Excel 打开）
|-- kaks_summary.xlsx        # Excel 汇总（多 sheet：结果/汇总/选择分布/质量分布/解释）
|-- summary_stats.json       # 统计汇总（均值/中位数/分位数/选择分布等）
|-- 99_logs/
|   +-- kaks_analysis.log    # 运行日志
+-- 00_pipeline_info/        # 流程元数据目录（预留）
```

## 结果解读 | Interpreting Results { #interpreting }

- **`kaks_detailed.tsv` / `.csv`**：核心结果。关键列：
  - `Ka`、`Ks`、`Ka/Ks`（即 omega）：比值是结论核心
  - `Selection_Type`：自动分类（Strong_Negative / Moderate_Negative / Neutral / Weak_Positive / Strong_Positive）
  - `Significance`：结合 Fisher 检验 P 值判断正/负选择是否显著
  - `Quality_Flag`：Good/Warning/Poor（异常高 omega、异常低 Ks、缺失值会被标记）
- **`kaks_summary.xlsx`**：适合快速浏览，Results sheet 是完整表，其余 sheet 是各类分布汇总。
- **`summary_stats.json`**：程序可读的统计结果，含平均 omega、选择压力分布、生物学解释文本。
- **好坏判据**：绝大多数基因 omega 应明显 <1（受负选择约束）；出现一批 omega >1 且 P<0.05 的基因，才提示正选择。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规全基因组 Ka/Ks 分析**：默认 GMYN，什么都不用改。
- **要和某篇老文献口径一致**：按对方方法选 `--method`（如经典 NG、YN）。
- **序列质量差、想先看看到底哪些被剔**：加 `--verbose`，日志里会逐条列出剔除原因（字符/N比例/长度/密码子）。
- **调试中间 AXT 文件**：加 `--keep-temp`，会在 `<output>/tmp` 下保留。
- **配对特别多（>1000 对）**：程序只发警告不会停，但耗时明显变长，建议分批。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--fasta1, -1` | 必填 | Path | 第一个FASTA文件｜First FASTA file |
| `--fasta2, -2` | 必填 | Path | 第二个FASTA文件｜Second FASTA file |
| `--pairs, -p` | 必填 | Path | 序列配对文件｜Sequence pair file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--method, -m` | `GMYN` | GMYN/MYN/YN/NG/LWL/LPB/MLWL/MLPB/GY/MS/MA/GNG/GLWL/GLPB/GMLWL/GMLPB/GYN | 计算方法｜Calculation method |
| `--kaks-path` | `KaKs_Calculator` | str | KaKs_Calculator软件路径｜Path to KaKs_Calculator executable |
| `--verbose, -v` | — |  | 启用详细日志｜Enable verbose logging |
| `--temp-dir` | — | Path | 自定义临时目录｜Custom temporary directory |
| `--keep-temp` | — |  | 保留临时文件｜Keep temporary files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-1, --fasta1, --species1` | 必填 | str | 第一个FASTA文件 (物种1 CDS序列)｜First FASTA file (species 1 CDS sequences) |
| `-2, --fasta2, --species2` | 必填 | str | 第二个FASTA文件 (物种2 CDS序列)｜Second FASTA file (species 2 CDS sequences) |
| `-p, --pairs, --pair-file` | 必填 | str | 序列配对文件 (TSV/CSV格式)｜Sequence pair file (TSV/CSV format) |
| `-o, --output, --output-dir` | 必填 | str | 输出目录｜Output directory for results |
| `-m, --method, --calc-method` | `GMYN` | GMYN/MYN/YN/NG/LWL/LPB/MLWL/MLPB/GY/MS/MA/GNG/GLWL/GLPB/GMLWL/GMLPB/GYN | 计算方法 (默认: GMYN)｜Calculation method (default: GMYN) |
| `-t, --threads` | `12` | int | 线程数 (默认: 12)｜Thread count (default: 12) |
| `--kaks-path, --calculator-path` | `KaKs_Calculator` | str | KaKs_Calculator可执行文件路径｜Path to KaKs_Calculator executable |
| `-v, --verbose, --debug` | — | store_true | 启用详细日志记录｜Enable verbose logging |
| `--temp-dir, --tmp-dir` | — | str | 自定义临时目录｜Custom temporary directory |
| `--keep-temp, --no-cleanup` | — | store_true | 保留临时文件 (用于调试)｜Keep temporary files (for debugging) |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- KaKs_Calculator 2.0（默认路径 `~/miniforge3/envs/phylo/bin/KaKs_Calculator`，conda 环境 `phylo`，可用 `KAKS_PATH` 环境变量或 `--kaks-path` 覆盖）
- seqkit（默认路径 `~/miniforge3/envs/misc/bin/seqkit`，conda 环境 `misc`，用于抽序列子集）
- Python 包：BioPython、pandas、openpyxl（写 Excel）

## 常见问题 | FAQ { #faq }

**Q1：为什么有的配对结果里没有出现？**
输入质控阶段被剔除的配对不会进入结果：序列含非法字符、N 比例 >10%、长度超出 50-50000、长度不是 3 的倍数，都会导致该配对被跳过。加 `--verbose` 看日志里的「过滤 X 个配对」与逐条失败原因。

**Q2：序列长度不是 3 的倍数会怎样？**
质控直接判该序列无效（密码子不完整），相关配对被剔除。请确认输入是完整 CDS（去掉 UTR、内含子）。注意：程序内部写 AXT 时对「配对里两条长度不一致」会**自动截断到较短者**（对齐到 3 的倍数），这只针对长度一致性问题，不是密码子问题。

**Q3：KaKs_Calculator 要 AXT 格式，我要自己准备吗？**
不需要。程序自动把两个 FASTA + 配对表转成 AXT 格式喂给 KaKs_Calculator，你只需给 FASTA 和配对表。

**Q4：换参数重跑会不会复用旧结果？**
不会，本模块**没有断点续传**，每次从头计算。但临时目录会落到 `<output>/tmp` 下并在结束后清理，重跑是干净的重来。

**Q5：omega 都是 NaN 或 0 怎么回事？**
通常意味着 Ks=0（两序列几乎相同，无法估计替换率）或 KaKs_Calculator 对该对计算失败。这类对会被标成 Poor（Missing_Values），属正常现象，重点看有效计算的那部分。
