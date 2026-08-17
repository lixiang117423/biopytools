# Phytophthora 效应子鉴定 | Phytophthora Effector Identification

一句话理解：**自动从一批病原(卵菌)蛋白序列里，把「效应子」(帮助侵染宿主的分泌蛋白)挑出来，覆盖 RxLR、CRN、NLP、Protease、SCP、elicitin、YxSL 七类，并给出每条候选的证据(信号肽、HMM 命中、基序位置、跨膜域)。**

## 功能概述 | Overview { #overview }

- 鉴定七类效应子：RxLR、CRN、NLP、Protease、SCP、elicitin、YxSL
- 支持单文件或多样本目录输入，多样本按样本独立运行并自动合并汇总
- SignalP 3.0 / 6.0 信号肽预测，可二选一或 both 取并集(更全)
- 各类型用 hmmsearch 做 HMM 搜索；RxLR 额外有 BLASTP 同源补充和正则基序兜底
- RxLR 类型额外做 TMHMM 跨膜域过滤(效应子应是分泌蛋白，不应有跨膜域)
- 所有 HMM 模型内置(bundled)，无需自备
- 断点续传：已完成步骤自动跳过

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools phyto-effector -i proteins.fa -o effector_out
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 效应子(effector) | 病原菌分泌到宿主里「搞破坏 / 躲免疫」的蛋白 |
| 信号肽(signal peptide) | 蛋白 N 端一段「出关凭证」，指导蛋白分泌到细胞外 |
| RxLR / EER | 卵菌效应子标志性的短序列基序，像「身份标签」 |
| CRN | 另一类效应子，标志是 LFLAK 基序 |
| HMM | 一种「序列模式模板」，能模糊匹配同一家族的蛋白 |
| hmmsearch | 用 HMM 去蛋白库里搜同源序列的工具 |
| E-value / score | 匹配可信度打分，E-value 越小 / score 越大越可信 |
| 跨膜域(TM) | 蛋白嵌进细胞膜的部分；效应子是分泌蛋白，正常不该有 |
| 基序(motif) | 序列上固定的小段特征 |
| 信号肽切割位点 | 信号肽被剪掉的位置，跨膜域过滤时用来判断 TM 是否在信号肽区外 |

## 输入 | Input { #input }

一个蛋白质 FASTA 文件，或一个装 FASTA 文件的目录(自动识别 `.fa` `.fasta` `.faa` `.pep` `.protein` `.prot`)：

- 单文件：所有序列视为一个样本，结果直接写在输出目录下
- 目录：每个 FASTA 文件视为一个样本，结果写到 `<输出目录>/<样本名>/`，并在顶层合并各样本候选

```text
>protein1
MSKTRKVL...
>protein2
MAAQ...
```

## 参数说明 | Parameters { #parameters }

### 输入输出 | Input & output { #parameters-io }

**通俗理解|In plain words:** `-i` 给 FASTA(文件或目录)，`-o` 给输出目录。目录里有多个 FASTA 就是多样本模式，会自动合并汇总。

### SignalP 信号肽预测 | SignalP { #parameters-signalp }

**通俗理解|In plain words:** 信号肽是「这蛋白会分泌」的证据，几乎所有效应子都要求有信号肽。`--signalp-version` 选 3、6 或 both(默认，两者都跑取并集，最全也最慢)；`--organism` 一般固定 eukarya；`--signalp-mode` 决定 SignalP 6 的精度/速度取舍，默认 slow-sequential 最准。`--skip-signalp` 在只想纯跑 HMM、不关心信号肽时打开(会放宽候选标准)。这些路径参数(SignalP / hmmsearch 的安装位置)一般不用动，走默认 conda 环境路径即可。

### HMM 搜索 | HMM search { #parameters-hmm }

**通俗理解|In plain words:** 每个效应子类型用对应的 HMM 模板去库里搜同源蛋白。各 `--xxx-hmm` 默认都用内置模型，一般不用动；`--score-threshold` 是 HMM 打分下限，默认 0(不过滤)，调大能减少误报但也可能漏掉真效应子。`-e/--evalue` 已弃用，改用 `--score-threshold`。

### RxLR 专属 | RxLR-specific { #parameters-rxlr }

**通俗理解|In plain words:** `--use-wy-domain` 让 RxLR 额外搜索 WY 结构域(某些 RxLR 效应子带这个结构)，默认关。RxLR 的 BLASTP 同源补充用内置参考序列，一般不用动。

### 运行资源 | Compute { #parameters-compute }

**通俗理解|In plain words:** `--threads` 控制 hmmsearch / BLASTP / SignalP 的并行度，越大越快，默认 12。

## 分析流程 | Pipeline { #pipeline }

```text
输入 FASTA(单文件或目录)
    │
    ▼
SignalP 信号肽预测(3.0 / 6.0 / both，各类型共享一次)
    │
    ├── RxLR：hmmsearch + BLASTP 同源 + 正则基序 → TMHMM 跨膜域过滤 → 基序注释
    ├── CRN：hmmsearch → LFLAK 基序注释
    └── NLP / Protease / SCP / elicitin / YxSL：hmmsearch → (可选)基序注释
    │
    ▼
各类型候选 TSV(多样本时顶层合并 + software_versions.yml)
```

## 输出 | Output { #output }

单样本模式(输出目录下，各类型一个子目录)：

```text
effector_out/
├── rxlr/
│   ├── 01_signalp/                 # 信号肽结果(全类型共享，只跑一次)
│   ├── 02_hmmsearch/               # rxlr.domtblout
│   ├── 03_blastp/                  # BLASTP 同源结果
│   ├── 04_wy_domain/               # WY 结构域(仅 --use-wy-domain)
│   ├── 05_tmhmm/                   # TMHMM 跨膜域结果 + 过滤清单
│   └── 06_candidates/rxlr_candidates.tsv   # RxLR 候选(核心结果)
├── crn/
│   ├── 02_hmmsearch/               # crn.domtblout
│   └── 03_candidates/crn_candidates.tsv     # CRN 候选
├── nlp/04_nlp/nlp_candidates.tsv
├── protease/05_protease/protease_candidates.tsv
├── scp/06_scp/scp_candidates.tsv
├── elicitin/07_elicitin/elicitin_candidates.tsv
├── yxsl/08_yxsl/yxsl_candidates.tsv
└── 99_logs/pipeline.log            # 运行日志
```

多样本模式：顶层为每个样本建子目录(内部结构同上)，并在顶层 `rxlr/rxlr_candidates_all_samples.tsv` 等位置放合并后的候选，`00_pipeline_info/software_versions.yml` 记录软件版本。

核心结果文件是各类型的 `*_candidates.tsv`，关键列：`Effector_Type`、`Sample`、`Sequence_ID`、`Sequence`、`SignalP`、`SP_Cleavage_Site`、`Source`(HMM/BLASTP/Regex)、`HMM_Score`、`RxLR`/`EER`(Yes/No + 位置)、`TMHMM_TMs` 等。

## 结果解读 | Interpreting Results { #interpreting-results }

- 候选 TSV 的 `SignalP=Yes` 是「有分泌证据」；`HMM=Yes` 是「命中了该类型模板」；RxLR 类型的 `RxLR=Yes` / `EER=Yes` 是「有标志基序」
- RxLR 候选经过 TMHMM 过滤后，不应有信号肽区外的跨膜域(`TMHMM_TMs` 为 0 或在信号肽区内)
- 三个证据(信号肽 + HMM + 基序)凑得越齐，越像真效应子；只有单一线索的要谨慎
- Protease 类型要求必须有信号肽(`require_sp`)来控制误报，这是刻意的
- 某类型 candidates 为空：可能是这个物种确实没有该类效应子，也可能是 HMM 阈值 / 数据质量问题

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规鉴定：全默认即可(both 信号肽 + 全类型 HMM)，直接看各类型 candidates TSV
- 只想快速预览：`--signalp-version 6`(只跑 SignalP 6，跳过 3.0)或 `--skip-signalp` 先看 HMM 命中面
- 误报偏多想收紧：把 `--score-threshold` 从 0 往上调(如 10~20)，牺牲一点召回换精度
- 做 RxLR 专项且关心 WY 结构域：加 `--use-wy-domain`
- 有大量序列(数万条)：SignalP 3.0 会自动分块处理，无需手动干预

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件或目录｜Input FASTA file or directory |
| `-o, --output-dir` | `./phyto_effector_output` |  | 输出目录｜Output directory |
| `--skip-signalp` | — |  | 跳过SignalP预测｜Skip SignalP prediction |
| `--signalp-path` | `~/miniforge3/envs/protein/bin/signalp6` |  | SignalP程序路径｜SignalP program path |
| `--organism` | `eukarya` | eukarya/other | 生物类型｜Organism type |
| `--signalp-mode` | `slow-sequential` | fast/slow/slow-sequential | SignalP运行模式｜SignalP run mode |
| `--signalp-version` | `both` | 3/6/both | SignalP版本｜SignalP version (3/6/both) |
| `--signalp3-path` | `~/miniforge3/envs/signalp_v.3.0b/bin/signalp` |  | SignalP 3.0程序路径｜SignalP 3.0 program path |
| `--signalp3-sprob-threshold` | `0.9` | float | SignalP 3.0 HMM Sprob阈值｜SignalP 3.0 HMM Sprob threshold |
| `--hmmsearch-path` | `~/miniforge3/envs/protein/bin/hmmsearch` |  | hmmsearch程序路径｜hmmsearch program path |
| `--rxlr-hmm` | — |  | RxLR HMM文件(默认内置)｜RxLR HMM file (default: bundled) |
| `--use-wy-domain` | — |  | 同时搜索WY结构域｜Also search WY domain |
| `--rxlr-wy-hmm` | — |  | WY HMM文件(默认内置)｜WY HMM file (default: bundled) |
| `--crn-hmm` | — |  | CRN HMM文件(默认内置)｜CRN HMM file (default: bundled) |
| `--nlp-hmm` | — |  | NLP HMM文件(默认内置)｜NLP HMM file (default: bundled) |
| `--protease-hmm` | — |  | Protease HMM文件(默认内置)｜Protease HMM file (default: bundled) |
| `--scp-hmm` | — |  | SCP HMM文件(默认内置)｜SCP HMM file (default: bundled) |
| `--elicitin-hmm` | — |  | elicitin HMM文件(默认内置)｜elicitin HMM file (default: bundled) |
| `--yxsl-hmm` | — |  | YxSL HMM文件(默认内置)｜YxSL HMM file (default: bundled) |
| `-e, --evalue` | `1e-05` | float | E-value阈值(已弃用)｜E-value threshold (deprecated) |
| `--score-threshold` | `0.0` | float | HMM score阈值｜HMM score threshold |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件或目录｜Input FASTA file or directory |
| `-o, --output-dir` | `./phyto_effector_output` |  | 输出目录｜Output directory (default: ./phyto_effector_output) |
| `--skip-signalp` | — | store_true | 跳过SignalP预测｜Skip SignalP prediction |
| `--signalp-path` | `~/miniforge3/envs/protein/bin/signalp6` |  | SignalP程序路径｜SignalP program path |
| `--organism` | `eukarya` | eukarya/other | 生物类型｜Organism type (default: eukarya) |
| `--signalp-mode` | `slow-sequential` | fast/slow/slow-sequential | SignalP运行模式｜SignalP run mode (default: slow-sequential) |
| `--signalp-version` | `both` | 3/6/both | SignalP版本｜SignalP version: 3, 6, or both (default: both) |
| `--signalp3-path` | `~/miniforge3/envs/signalp_v.3.0b/bin/signalp` |  | SignalP 3.0程序路径｜SignalP 3.0 program path |
| `--signalp3-sprob-threshold` | `0.9` | float | SignalP 3.0 HMM Sprob阈值｜SignalP 3.0 HMM Sprob threshold (default: 0.9) |
| `--hmmsearch-path` | `~/miniforge3/envs/protein/bin/hmmsearch` |  | hmmsearch程序路径｜hmmsearch program path |
| `--blastp-path` | `~/miniforge3/envs/annot/bin/blastp` |  | blastp程序路径｜blastp program path |
| `--rxlr-blastp-queries` | — |  | RxLR BLASTP查询序列FASTA(默认内置)｜RxLR BLASTP query FASTA (default: bundled) |
| `--tmhmm-path` | `~/miniforge3/envs/protein/bin/tmhmm` |  | tmhmm程序路径｜tmhmm program path |
| `--rxlr-hmm` | — |  | RxLR HMM文件路径(默认内置)｜RxLR HMM file path (default: bundled) |
| `--use-wy-domain` | — | store_true | 同时搜索WY结构域｜Also search WY domain |
| `--rxlr-wy-hmm` | — |  | WY HMM文件路径(默认内置)｜WY HMM file path (default: bundled) |
| `-e, --evalue` | `1e-05` | float | E-value阈值(已弃用，使用--score-threshold)｜E-value threshold (deprecated) |
| `--score-threshold` | `0.0` | float | HMM score阈值｜HMM score threshold (default: 0.0) |
| `--crn-hmm` | — |  | CRN HMM文件路径(默认内置)｜CRN HMM file path (default: bundled) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--…-hmm` | — |  | … HMM文件路径(默认内置)｜… HMM file path (default: bundled) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- signalp6：`~/miniforge3/envs/protein/bin/signalp6`(conda 环境 protein)
- signalp(3.0)：`~/miniforge3/envs/signalp_v.3.0b/bin/signalp`(conda 环境 signalp_v.3.0b)
- hmmsearch：`~/miniforge3/envs/protein/bin/hmmsearch`(conda 环境 protein)
- blastp / makeblastdb：`~/miniforge3/envs/annot/bin/blastp`(conda 环境 annot)
- tmhmm：`~/miniforge3/envs/protein/bin/tmhmm`(conda 环境 protein)
- Python 3 + pandas + yaml

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持。SignalP 结果、各类型 domtblout、BLASTP 结果、TMHMM 过滤清单等按「输出文件存在性」跳过。换参数重跑前先删对应旧产物，否则会复用旧结果。

**Q2：SignalP 3.0 对超长蛋白库会失败吗？**
不会。超过 3800 条序列会自动按 300 条/块切分，单块失败还会递归拆半重试，运行结束自动清理临时文件。

**Q3：`-e/--evalue` 还有用吗？**
已弃用，现在用 `--score-threshold` 控制 HMM 过滤。`-e` 只影响 RxLR 的 BLASTP 步骤。

**Q4：为什么 Protease 候选里有些没写进结果？**
Protease 家族扩充到 203 个 Pfam 家族后命中面变广，特意加了「必须有信号肽」的要求来控误报，无信号肽的命中会被跳过。

**Q5：多样本目录里文件重名会怎样？**
程序会检查样本名(stem)是否冲突，重名直接报错，改文件名再跑。

**Q6：只想重新跑某一种效应子类型？**
目前是全类型一起跑；但配合断点续传，已完成类型的产物不会被重算，成本可控。