# CPhasing 基因组分相与挂载 | CPhasing Genome Phasing and Scaffolding

一句话理解：**用 Hi-C / Pore-C 数据把多倍体基因组里「来自不同亲本」的两套（或多套）染色体分开，并挂载到染色体级别**，解决「多倍体组装分不清哪条染色体属于哪个亚基因组」的问题。

## 功能概述 | Overview

- 封装 CPhasing 工具集，支持 pipeline、mapper、alleles、hypergraph、hyperpartition、plot、chimeric、collapse 等全部子命令
- 默认 pipeline 模式：一命令完成 Hi-C 比对、分相（phasing）、染色体挂载全流程
- 内置亚基因组聚类（`--enable-haplotype-cluster`），默认对 pipeline/scaffolding 开启
- 支持多种分相模式（phasing / haploid / basal / basal_withprune）与预设（precision / sensitive / very-sensitive / nofilter）
- 对上游 CPhasing「退出码 0 但实际出错」的静默失败做了检测兜底

## 快速开始 | Quick Start

```bash
biopytools cphasing -i genome.fa --hic1 R1.fq.gz --hic2 R2.fq.gz -t 12 -n 16:2
```

注意：运行前必须先激活 CPhasing 环境（见「依赖」）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 多倍体 | 基因组里有多套同源染色体（如异源四倍体有 A、B 两套） |
| 分相（phasing） | 把混在一起的「两套染色体」区分开，各自还原 |
| 亚基因组 | 多倍体里来自某个祖先亲本的那套染色体（如 A 亚基因组） |
| 挂载（scaffolding） | 把碎片按正确顺序串成完整染色体 |
| Hi-C / Pore-C | 反映 DNA 空间邻近性的测序技术，用来判断哪些片段该连在一起 |
| 分组数（-n） | 期望分成几组，形如 `16:2` 表示「16 条染色体、2 个亚基因组」；`0` 为自动 |

## 输入 | Input

### 基因组

FASTA 格式的组装（pipeline 模式必需）。用 `-i` / `--input` 指定。

### Hi-C 读段

两个配对文件，用 `--hic1` / `--hic2` 指定（支持 .fq.gz 压缩）。

```text
biopytools cphasing -i genome.fa --hic1 R1.fq.gz --hic2 R2.fq.gz -t 12 -n 16:2
```

## 参数说明 | Parameters

### 子命令 | Subcommand

**通俗理解|In plain words:** CPhasing 是一套工具集，第一个位置参数决定跑哪个子命令，不写默认 `pipeline`（最常用的完整流程）。想只跑某一步（如只要比对用 mapper），写对应子命令名即可。

### 必需参数（pipeline 模式）| Required

**通俗理解|In plain words:** pipeline 模式下，基因组、Hi-C R1、Hi-C R2 三者缺一不可。其他子命令不需要这些。

### 分组与模式 | Groups and mode

**通俗理解|In plain words:** `-n` 是「分成几组」的关键参数：`16:2` 意思是最终 16 条染色体、分属 2 个亚基因组，`0` 让程序自动判断。`--mode` 选分相策略（多倍体分相用默认 phasing；纯合/单倍体场景用 haploid；basal 用于祖先型组装）。`--preset` 是灵敏度预设，precision 最严格、nofilter 不过滤，一般用默认 precision。

### Hi-C 比对器 | Hi-C aligner

**通俗理解|In plain words:** `--hic-aligner` 选比对工具，默认 `_chromap`（内部 chromap）速度最快；可选 chromap / minimap2 / bwa-mem2。`--hic-mapper-k` / `--hic-mapper-w` 是比对器的 kmer 与窗口参数，一般不用动。`--mapping-quality` 是最低比对质量门槛。

### 步骤控制与高级选项 | Step control and advanced

**通俗理解|In plain words:** `--steps` / `--skip-steps` 用于只跑或跳过指定步骤。`--hcr` 启用高置信区域分析；`--pattern` 指定酶切位点模式（如 AAGCTT）；`--low-memory` 在内存紧张时开启。这些都属于进阶选项，正常流程一般不用碰。

## 输出 | Output

CPhasing 的输出由上游工具生成在输出目录（默认 `./cphasing_output`）下，内容随子命令而异。pipeline 模式通常包含：

```text
cphasing_output/
├── 00_data/        # 数据准备与比对中间产物
├── ...             # CPhasing 各步骤输出（分组/分相/挂载结果）
├── fasta/          # 挂载后的染色体级 FASTA（按亚基因组分）
└── 99_logs/
    └── cphasing.log  # biopytools 封装层的运行日志
```

关键结果：分相与挂载后的染色体级 FASTA（通常按亚基因组分别输出），以及各步骤的中间统计。具体文件名以 CPhasing 版本为准，以 `99_logs/cphasing.log` 里记录的完整命令与输出路径为准。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看分相挂载结果，核心是「每个亚基因组是否拿到了预期的染色体数、挂载是否连续」。

- **染色体数目**：最终每个亚基因组的染色体数应符合预期（如 16:2 应得到两套各 16 条上下）。
- **挂载连续性**：染色体级 FASTA 里 N 的数量越少、contig 越长，说明挂载越连续。
- **分组是否合理**：日志里会记录分组信息，若某亚基因组染色体数明显偏多/偏少，检查 -n 分组数是否设对。
- **静默失败**：本封装层会扫描 stderr 里的 Traceback 等异常模式，即使 CPhasing 退出码为 0 也会判失败，留意日志里的「命令静默失败」提示。

## 参数选择建议 | Parameter Guidance

- **-n 分组数**：异源四倍体常见 `染色体数:2`；不确定就用 `0` 自动。
- **--mode**：异源多倍体分相用 phasing；分相完成后要合并成单倍体用 haploid；basal 用于构建祖先型参考。
- **--hic-aligner**：大基因组默认 `_chromap` 最快；比对异常时可换 minimap2 或 bwa-mem2 排查。
- **--preset**：首次跑用 precision；数据质量差、比对上不去时可换 sensitive 或 very-sensitive。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 基因组FASTA文件｜Genome FASTA file |
| `--hic1` | — |  | Hi-C R1 reads文件｜Hi-C R1 reads file |
| `--hic2` | — |  | Hi-C R2 reads文件｜Hi-C R2 reads file |
| `-n, --groups` | `0` |  | 分组数｜Number of groups (例如: "8:4" ｜ e.g., "8:4", "0"=自动｜auto) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--mode` | `phasing` | phasing/haploid/basal/basal_withprune | 分相模式｜Phasing mode |
| `--preset` | `precision` | precision/sensitive/very-sensitive/nofilter | 分析预设｜Analysis preset |
| `-o, --output-dir` | `./cphasing_output` |  | 输出目录｜Output directory |
| `--steps` | — |  | 运行指定步骤｜Run specified steps (例如: "1,2,3" ｜ e.g., "1,2,3") |
| `--skip-steps` | — |  | 跳过步骤｜Skip steps (例如: "1,2" ｜ e.g., "1,2") |
| `--hic-aligner` | `_chromap` | _chromap/chromap/minimap2/bwa-mem2 | Hi-C比对器｜Hi-C aligner |
| `--hic-mapper-k` | — | int | Hi-C mapper kmer大小｜Hi-C mapper kmer size |
| `--hic-mapper-w` | — | int | Hi-C mapper窗口大小｜Hi-C mapper window size |
| `--mapping-quality` | `0` | int | 最小比对质量｜Minimum mapping quality |
| `--hcr` | — |  | 启用高置信区域｜Enable high confidence regions |
| `--pattern` | — |  | 酶切位点模式｜Restriction enzyme pattern (例如: AAGCTT) |
| `--low-memory` | — |  | 低内存模式｜Low memory mode |
| `--no-haplotype-cluster` | `False` |  | 禁用亚基因组聚类（默认开启，仅在pipeline/scaffolding子命令生效）｜Disable subgenome clustering (default ON; only applies to pipeline/scaffolding) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 基因组FASTA文件｜Genome FASTA file (pipeline模式必需｜required for pipeline) |
| `--hic1` | — |  | Hi-C R1 reads文件｜Hi-C R1 reads file (pipeline模式必需｜required for pipeline) |
| `--hic2` | — |  | Hi-C R2 reads文件｜Hi-C R2 reads file (pipeline模式必需｜required for pipeline) |
| `-n, --groups` | `0` |  | 分组数｜Number of groups (例如: "8:4" ｜ e.g., "8:4", "0"=自动｜auto) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (默认｜default: 12) |
| `--mode` | `phasing` | phasing/haploid/basal/basal_withprune | 分相模式｜Phasing mode (默认｜default: phasing) |
| `--preset` | `precision` | precision/sensitive/very-sensitive/nofilter | 分析预设｜Analysis preset (默认｜default: precision) |
| `-o, --output-dir` | `./cphasing_output` |  | 输出目录｜Output directory (默认｜default: ./cphasing_output) |
| `--steps` | — |  | 运行指定步骤｜Run specified steps (pipeline模式) |
| `--skip-steps` | — |  | 跳过步骤｜Skip steps (pipeline模式) |
| `--hic-aligner` | `_chromap` | _chromap/chromap/minimap2/bwa-mem2 | Hi-C比对器｜Hi-C aligner (默认｜default: _chromap) |
| `--hic-mapper-k` | — | int | Hi-C mapper kmer大小｜Hi-C mapper kmer size |
| `--hic-mapper-w` | — | int | Hi-C mapper窗口大小｜Hi-C mapper window size |
| `--mapping-quality` | `0` | int | 最小比对质量｜Minimum mapping quality (默认｜default: 0) |
| `--hcr` | — | store_true | 启用高置信区域｜Enable high confidence regions |
| `--pattern` | — |  | 酶切位点模式｜Restriction enzyme pattern (例如: AAGCTT) |
| `--low-memory` | — | store_true | 低内存模式｜Low memory mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- CPhasing（用 pixi 安装，**非普通 conda 环境**）
- 运行前必须先执行：`source ~/software/CPhasing_v0.3.0/bin/activate_cphasing`（路径按实际安装位置调整）
- cphasing-rs（Rust 二进制，CPhasing 内部 18 处调用，缺失会导致静默失败，通常随 activate_cphasing 一起设置到 PATH）
- 软件目录默认 `~/software/CPhasing/CPhasing-main`（可用环境变量 CPHASING_DIR 覆盖）

## 常见问题 | FAQ

**Q1：报「cphasing 不在 PATH 中」？**
CPhasing 用 pixi 安装，必须先 `source ~/software/CPhasing_v0.3.0/bin/activate_cphasing` 激活，普通 `conda activate` 无效。

**Q2：cphasing 能跑但结果不完整、没有明显报错？**
很可能是 cphasing-rs 不在 PATH。activate_cphasing 会一并设置它；缺失时 CPhasing 内部多处静默失败。用 `which cphasing-rs` 确认。

**Q3：断点续传支持吗？**
本封装层不额外做断点续传，直接委托给 CPhasing 自身（CPhasing 会跳过已完成步骤）。用 `--steps` / `--skip-steps` 可手动控制重跑范围。

**Q4：为什么默认加了 --enable-haplotype-cluster？**
亚基因组聚类（把染色体归到各亚基因组）对多倍体分相是关键一步，故对 pipeline/scaffolding 默认开启。想关闭用 `--no-haplotype-cluster`。

**Q5：想只跑比对或只跑某一步？**
用子命令直接指定，如 `biopytools cphasing mapper -i genome.fa --hic1 ... --hic2 ...`，其余参数经 `--` 之后透传给 CPhasing。

