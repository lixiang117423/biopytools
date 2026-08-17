# RepeatMask - 重复序列屏蔽 | Repeat Masking

一句话理解：**给一个基因组自动从头识别重复序列(RepeatModeler)，再用 RepeatMasker 把重复区「糊掉/标记」(屏蔽)**，产出屏蔽后的基因组和重复注释，供下游基因预测使用。

## 功能概述 | Overview { #overview }

- RepeatModeler 从头建重复库（可选 `-LTRStruct` 结构发现、`-quick` 快速模式）
- RepeatMasker 用该库做全基因组屏蔽，产出屏蔽基因组 + 注释
- 屏蔽模式可选 soft(小写) / hard(N) / x(X)
- 可选 Dfam/Repbase 物种库（`-s species`，默认使用 Dfam）做已知重复补充
- 断点续传：数据库、库文件、masked 文件存在即自动跳过（见 FAQ）

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools repeatmask -i genome.fa -o output
```

最小输入：一个基因组 FASTA。默认从头建库 + 软屏蔽。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 重复序列 | 基因组里反复出现的大量相似片段（主要是转座子），会干扰基因预测 |
| 从头建库(de novo) | 不看参考、直接从这份基因组「归纳」出重复序列代表库 |
| 屏蔽(masking) | 把重复区标记出来，让基因预测软件跳过；软屏蔽=转小写，硬屏蔽=换成 N |
| 软屏蔽(soft) | 重复区写成小写字母，序列还在、信息不丢，最常用 |
| 硬屏蔽(hard) | 重复区换成字母 N，彻底抹掉；`x` 模式换成 X（视觉更醒目） |
| Dfam / Repbase | 公共的重复序列/转座子数据库，按物种提供已知重复，补充从头库漏掉的 |
| LTR 结构发现 | 专门找「长末端重复」结构，提高 LTR 逆转座子的建库精度 |

## 输入 | Input { #input }

### 基因组 FASTA

必需，标准 FASTA 格式（可含多序列）：

```text
>chr1
ATGCGTACGTACGTAGCTA...
>chr2
TTGCAAGCTAGCATCGATC...
```

## 参数说明 | Parameters { #parameters }

### 流程控制 | Pipeline control

**通俗理解|In plain words:** `--skip-modeler` 跳过从头建库（你已有库、只想直接屏蔽时用）；`--no-dfam` 不用 Dfam 公共库（只靠从头库）；`--no-ltr` 关掉 RepeatModeler 的 LTR 结构发现（更快但 LTR 建库变差）；`--modeler-quick` 用 RepeatModeler 快速模式（牺牲一点灵敏度换速度）。**默认从头建库 + LTR 结构 + Dfam，一般不用动。**

本组参数：`--skip-modeler`、`--no-dfam`、`--no-ltr`、`--modeler-quick`（均为开关，默认关）。

### 屏蔽设置 | Masking options

**通俗理解|In plain words:** `--masking-mode` 决定「重复区怎么标」——soft(小写)最常用、序列不丢；hard(N)/x(X)彻底抹掉。`-s species` 指定物种名以用 Dfam/Repbase 的对应库（如 `-s maize`），不指定则只靠从头库。

本组参数：`--masking-mode`（默认 soft，可选 soft/hard/x）、`-s/--species`（默认无）。

### 运行与路径 | Runtime & paths

**通俗理解|In plain words:** `-t` 线程数；三个软件路径默认指向 conda 环境 `repeat`，**装好后一般不用动**，装在别处时才改。

本组参数：`-t/--threads`（默认 12）、`--repeatmodeler-path`（默认 `~/miniforge3/envs/repeat/bin/RepeatModeler`）、`--repeatmasker-path`（默认 `~/miniforge3/envs/repeat/bin/RepeatMasker`）、`--builddatabase-path`（默认 BuildDatabase）。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先把基因组做成 RepeatModeler 能读的数据库，再从头建库，最后用库屏蔽。

```text
基因组 FASTA
    |
    v
1. 基因组统计(序列数/长度/N50/GC)
    |
    v
2. BuildDatabase 建库(可 --skip-modeler 跳过) -> {base}_db*
    |
    v
3. RepeatModeler 从头建库(可选 -LTRStruct/-quick)
    -> {base}_db-families.fa
    |
    v
4. RepeatMasker 屏蔽(软/硬, 可选 Dfam 物种库)
    -> repeatmasker_output/<基因组>.masked 等
```

其中 `{base}` 为 `<基因组名>_repeat`（如 `genome.fa` → `genome_repeat`）。

## 输出 | Output { #output }

```text
<输出目录>/
├── {base}_db*                        # BuildDatabase 数据库文件(二进制)
├── {base}_db-families.fa             # RepeatModeler 从头重复库
├── builddatabase.log                 # 建库日志
├── repeatmodeler.log                 # RepeatModeler 日志
├── repeatmask.log                    # 主日志
└── repeatmasker_output/
    ├── <基因组>.out                  # RepeatMasker 注释(核心)
    ├── <基因组>.out.gff              # 注释 GFF
    ├── <基因组>.tbl                  # 统计表(重复分类占比)
    ├── <基因组>.masked               # 屏蔽后基因组(原始名)
    ├── <基因组>.masked.fa            # 屏蔽后基因组(标准名, 核心)
    └── repeatmasker_run.log          # RepeatMasker 日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 屏蔽后基因组（`<基因组>.masked.fa`）

**通俗理解|In plain words:** 最核心产物——重复区已按所选模式标记的基因组，直接喂给基因预测软件（BRAKER/AUGUSTUS 等）。

- soft 模式下重复区是小写字母；hard 下是 N；x 下是 X

### 2. 注释与统计（`<基因组>.out` / `.tbl`）

- `.out`：每个重复的位置、类型、分歧度明细；`.tbl`：按类型汇总的占比表
- **总重复占比**是核心指标（多数植物 30%-80%）；异常低提示建库漏检，异常高提示假阳性（如把基因区也标成重复）

### 3. 从头库（`{base}_db-families.fa`）

RepeatModeler 归纳出的重复序列代表库，可复用于其它相关物种的屏蔽。

## 参数选择建议 | Parameter Guidance { #guidance }

- **`--masking-mode`**：基因预测用默认 `soft`（信息不丢）；做比对/序列分析需彻底去除重复时用 `hard`
- **`-s species`**：有对应 Dfam 物种库时加上可补充已知重复；不确定物种名可先不加，纯从头库也能跑
- **`--modeler-quick`**：大基因组想先快速出结果时开；追求灵敏度保持默认
- **`--skip-modeler`**：已用 repeat-analyzer 等拿到库，只想屏蔽时开启，省时明显
- **`--threads`**：RepeatModeler/RepeatMasker 吃核，建议 12-32

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./repeatmask_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-s, --species` | — |  | 物种名称｜Species name for Dfam/Repbase database |
| `--repeatmodeler-path` | `~/miniforge3/envs/repeat/bin/RepeatModeler` |  | RepeatModeler路径｜RepeatModeler path |
| `--repeatmasker-path` | `~/miniforge3/envs/repeat/bin/RepeatMasker` |  | RepeatMasker路径｜RepeatMasker path |
| `--builddatabase-path` | `BuildDatabase` |  | BuildDatabase路径｜BuildDatabase path |
| `--skip-modeler` | — |  | 跳过RepeatModeler｜Skip RepeatModeler step |
| `--no-dfam` | — |  | 不使用Dfam数据库｜Do not use Dfam database |
| `--masking-mode` | `soft` | soft/hard/x | 屏蔽模式｜Masking mode |
| `--no-ltr` | — |  | 不运行LTR结构发现｜Do not run LTR structural discovery |
| `--modeler-quick` | — |  | RepeatModeler快速模式｜RepeatModeler quick mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./repeatmask_output` |  | 输出目录｜Output directory (default: ./repeatmask_output) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-s, --species` | — |  | 物种名称(用于Dfam/Repbase数据库)｜Species name for Dfam/Repbase database |
| `--repeatmodeler-path` | `~/miniforge3/envs/repeat/bin/RepeatModeler` |  | RepeatModeler路径｜RepeatModeler path |
| `--repeatmasker-path` | `~/miniforge3/envs/repeat/bin/RepeatMasker` |  | RepeatMasker路径｜RepeatMasker path |
| `--builddatabase-path` | `BuildDatabase` |  | BuildDatabase路径｜BuildDatabase path |
| `--skip-modeler` | — | store_true | 跳过RepeatModeler步骤｜Skip RepeatModeler step |
| `--no-dfam` | — | store_true | 不使用Dfam数据库｜Do not use Dfam database |
| `--masking-mode` | `soft` | soft/hard/x | 屏蔽模式｜Masking mode: soft(小写｜lowercase), hard(N), x(X) (default: soft) |
| `--use-ltr` | `True` | store_true | 运行LTR结构发现｜Run LTR structural discovery (default: enabled) |
| `--no-ltr` | — | store_true | 不运行LTR结构发现｜Do not run LTR structural discovery |
| `--modeler-quick` | — | store_true | RepeatModeler快速模式｜RepeatModeler quick mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- RepeatModeler + BuildDatabase（默认 `~/miniforge3/envs/repeat/bin/`）
- RepeatMasker（默认 `~/miniforge3/envs/repeat/bin/`）
- conda 环境名 `repeat`（模块会自动检测命令路径并 `conda run -n repeat` 包装）
- Dfam/Repbase 数据库（可选，`-s species` 时使用）

## 常见问题 | FAQ { #faq }

**Q1：重跑会跳过已完成的步骤吗？**
会。断点续传按产物存在性判断：数据库文件、`-families.fa` 库文件、masked+out 文件已存在即跳过对应步骤。换参数重跑前需删除对应旧产物（或换输出目录）。

**Q2：`--skip-modeler` 后没库文件怎么办？**
跳过 RepeatModeler 时 RepeatMasker 只会用 Dfam/Repbase 物种库（若给了 `-s`）。若两者都没有，屏蔽会缺库。建议要么不跳过，要么提供 `-s species`。

**Q3：报 RepeatModeler/RepeatMasker/BuildDatabase 未找到？**
检查 PATH，或用 `--repeatmodeler-path` / `--repeatmasker-path` / `--builddatabase-path` 指定完整路径。

**Q4：`--no-ltr` 有什么影响？**
关掉 RepeatModeler 的 `-LTRStruct` 步骤，LTR 逆转座子的建库质量会下降，但速度变快。植物基因组建议保持默认（开）。

**Q5：soft 和 hard 屏蔽怎么选？**
soft 保留下游信息（推荐）；hard 彻底抹掉重复，适合纯比对、防止重复区干扰的场合。
