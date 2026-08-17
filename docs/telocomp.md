# TeloComp 端粒鉴定 | TeloComp Telomere Identification

一句话理解：**用 ONT/HiFi 长读长测序数据，找每条染色体两端的端粒重复序列，判断组装是否「完整到头」**。
输入一条基因组 FASTA 和 ONT/HiFi reads，输出含端粒的 reads（BAM）、左右端端粒序列、以及每条染色体两端端粒 reads 数的报告。

## 功能概述 | Overview

- Filter_1：从 ONT/HiFi reads 中筛出含端粒重复序列的 reads（比对回基因组，输出 BAM）
- Filter_2：同时有 ONT 和 HiFi 时才运行，提取左右两端端粒序列（trim_L / trim_R）
- 端粒位置分析：统计每条染色体左右两端的端粒 reads 数，生成报告和汇总
- 端粒重复序列可自定义：植物默认 `CCCTAAA`，动物可改 `TTAGGG`
- 自动创建基因组索引（`.fai`），无需手动准备

## 快速开始 | Quick Start

```bash
biopytools telocomp -g genome.fa -o output --ont ont.fastq.gz --hifi hifi.fastq.gz
```

最小输入：一条基因组 FASTA + ONT 或 HiFi reads 至少一种（Filter_2 需两种都有，见 FAQ）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 端粒(telomere) | 染色体两端的「帽子」，由重复序列组成，保护染色体 |
| 端粒重复序列 | 一段反复出现的小序列，如植物 CCCTAAA、动物 TTAGGG |
| ONT | 牛津纳米孔测序，超长读长，适合看染色体两端 |
| HiFi | PacBio 高精度长读长，读长中等、精度高 |
| soft-clip | 比对时 read 两端「对不上基因组」被切掉的部分，端粒 reads 常见此特征 |
| BAM | 比对结果的标准二进制格式 |

## 输入 | Input

### 基因组 FASTA

待检测的组装。程序会自动用 samtools 生成 `.fai` 索引。

### ONT / HiFi reads

至少提供一种（`--ont` 或 `--hifi`）。**Filter_2 步骤需要 ONT 和 HiFi 同时提供**，只有一种时会跳过 Filter_2（见 FAQ）。

## 参数说明 | Parameters

### 数据输入 | Input data

**通俗理解|In plain words:** `--ont` 和 `--hifi` 是两套可选的 reads，可以只给一种，也可以都给。都给时能跑完整的 Filter_1 + Filter_2 + 端粒位置分析；只给一种则只跑 Filter_1 和位置分析。

### 端粒重复序列 | Telomere motif

**通俗理解|In plain words:** `-m/--motif` 是端粒的重复单元，**植物默认 `CCCTAAA`，动物要改成 `TTAGGG`**，一定要按物种改对，否则一个端粒都找不到。`-M/--motif-num` 是重复单元的碱基数（默认 7，对应 CCCTAAA 的 7 个碱基），改 motif 时通常也要同步改它。

### 覆盖度与流程控制 | Coverage & control

**通俗理解|In plain words:** `-c/--coverage`（0-100）是 Filter_2 的覆盖度参数，默认 100，**一般不用动**。`--skip-filter` 跳过 Filter 步骤（只做位置分析，适合已有 BAM 的场景）；`--no-visualization` 跳过可视化/位置分析汇总。

## 分析流程 | Pipeline

```text
基因组 FASTA + ONT/HiFi reads
    │
    ▼
步骤1: 检查/创建基因组索引(samtools faidx)
    │
    ▼
步骤2: Filter_1 - 筛出含端粒的reads → filter1_ont.bam / filter1_hifi.bam
    │
    ▼
步骤3: Filter_2 - 提取左右端端粒序列(需ONT+HiFi都有) → filter2_output/trim_L, trim_R
    │
    ▼
步骤4: 端粒位置分析 - 统计每条染色体两端端粒reads数 → telomere_positions.txt
    │
    ▼
步骤5: 汇总 → telomere_summary.txt
```

## 输出 | Output

```text
output/
├── filter1_ont.bam         # ONT 含端粒 reads 的比对结果
├── filter1_hifi.bam        # HiFi 含端粒 reads 的比对结果
├── filter2_output/
│   ├── trim_L/             # 左端端粒序列(.fa)
│   └── trim_R/             # 右端端粒序列(.fa)
├── telomere_positions.txt  # 端粒位置检测报告(每条染色体两端reads数)
├── telomere_summary.txt    # 结果汇总
├── telocomp.log            # 运行日志
└── genome.fa.fai           # 基因组索引(自动生成)
```

- `telomere_positions.txt`：核心结果，逐条染色体列出左/右端各有多少 ONT、HiFi 端粒 reads
- `filter2_output/trim_L` / `trim_R`：实际提取出的左右端端粒序列，可用于后续端粒长度分析
- `filter1_ont.bam` / `filter1_hifi.bam`：含端粒 reads 的比对，做位置分析的基础

## 结果解读 | Interpreting Results

- **染色体两端都有端粒 reads = 组装「完整到头」**，是最理想的结果
- **只有一端有、另一端没有**：可能该端组装不完整，或端粒重复序列在该物种不匹配（先核对 `-m`）
- **某条染色体完全没有端粒 reads**：可能是中间 contig（本就没有端粒），或 motif 不对
- **ONT 和 HiFi 的 reads 数可以互相印证**：两者都在同一端检出，可信度更高
- 报告中「未检测到」不代表一定缺失，端粒区在测序/组装里本来就容易被截断

## 参数选择建议 | Parameter Guidance

- **植物**：默认 `CCCTAAA` 直接用
- **动物**：`-m TTAGGG -M 6`（TTAGGG 是 6 个碱基）
- **只有一种 reads**：只跑 Filter_1 + 位置分析即可，Filter_2 会自动跳过
- **已有 Filter_1 的 BAM、只想重跑位置分析**：`--skip-filter`
- **只要 BAM 不要报告**：`--no-visualization`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `--ont` | — |  | ONT数据文件｜ONT data file |
| `--hifi` | — |  | HiFi数据文件｜HiFi data file |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-c, --coverage` | `100` | int | 覆盖度参数(0-100)｜Coverage parameter (0-100) |
| `-m, --motif` | `CCCTAAA` |  | 端粒重复序列(植物默认CCCTAAA, 动物TTAGGG)｜Telomeric repeat sequence (plant: CCCTAAA, animal: TTAGGG) |
| `-M, --motif-num` | `7` | int | 端粒重复序列碱基数｜Number of bases in telomere motif |
| `--skip-filter` | — |  | 跳过Filter步骤｜Skip filter steps |
| `--no-visualization` | — |  | 跳过可视化｜Skip visualization |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--ont` | — |  | ONT数据文件｜ONT data file |
| `--hifi` | — |  | HiFi数据文件｜HiFi data file |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-c, --coverage` | `100` | int | 覆盖度参数(0-100)｜Coverage parameter (0-100) |
| `-m, --motif` | `CCCTAAA` |  | 端粒重复序列｜Telomeric repeat sequence |
| `-M, --motif-num` | `7` | int | 端粒重复序列碱基数｜Number of bases in telomere motif |
| `--skip-filter` | — | store_true | 跳过Filter步骤｜Skip filter steps |
| `--no-visualization` | — | store_true | 跳过可视化｜Skip visualization |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- conda 环境 `telocomp`（默认 `~/miniforge3/envs/telocomp`），提供 python 运行环境和 samtools
- TeloComp 工具集（默认 `~/software/telocomp/TeloComp-1.0.0/bin`），含 `telocomp_Filter_1`、`telocomp_Filter_2`
- GenomeSyn（默认 `~/software/GenomeSyn/GenomeSyn-main/GenomeSyn-1.2.7/bin`）
- Python 包 `pysam`（位置分析/统计 BAM 时使用）

## 常见问题 | FAQ

**Q1：动物物种一个端粒都找不到？**
端粒重复序列默认是植物的 `CCCTAAA`。动物请改用 `-m TTAGGG -M 6`。

**Q2：为什么 Filter_2 没跑、filter2_output 是空的？**
Filter_2 需要 **ONT 和 HiFi 同时提供**。只有一种 reads 时会自动跳过并创建空目录，这是预期行为——你仍可拿到 Filter_1 的 BAM 和端粒位置报告。

**Q3：能看到端粒图吗？**
当前模块只实现了 Filter 步骤（端粒 reads 识别/过滤）和位置分析，**未实现完整的 Complement 可视化绘图**（需要额外 WGS reads、NextPolish 和 Assembly 步骤）。报告 `telomere_positions.txt` 是主要产出。

**Q4：有断点续传吗？**
没有。重跑会重新执行各步骤。索引已存在时不会重复创建，但 Filter 和位置分析会重做。

**Q5：报错找不到 TeloComp bin / GenomeSyn bin 怎么办？**
这两个是默认路径 `~/software/telocomp/TeloComp-1.0.0/bin` 和 `~/software/GenomeSyn/.../GenomeSyn-1.2.7/bin`，且必须真实存在（代码会校验）。若你的安装路径不同，需先确认对应目录存在。
