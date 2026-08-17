# PanEDTA 泛基因组转座子注释 | PanEDTA Pan-genome TE Annotation

一句话理解：把多个基因组的转座子(跳跃基因)放到一起联合鉴定，生成一个跨基因组非冗余的"泛转座子库"，再注释回各个基因组。

## 功能概述 | Overview

- 对多个基因组做转座子(TE)的联合鉴定与注释，而非逐个单跑
- 通过上游 panEDTA.sh(EDTA 套件的一部分)执行
- 可选 CDS 文件用于排除"长得像转座子的基因"，可选人工筛选的 TE 库做引导
- 可控是否执行全基因组注释、是否覆盖旧结果
- 输出泛转座子库及各基因组的 TE 注释

## 快速开始 | Quick Start

```bash
biopytools panedta -i genomes.txt -t 24
```

`-i` 是基因组列表文件(每行一个基因组 FASTA 路径)。输出默认写到 `./panedta_output`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 转座子(TE) | "跳跃基因"，能在基因组里复制自己、搬来搬去的 DNA 片段 |
| 泛转座子库(panTE library) | 把多个基因组里鉴定出的 TE 去重合并成的一份"代表库" |
| 全长拷贝(fl-copy) | 结构完整、没被截断的 TE 拷贝；用它筛掉"残缺/半截"的假阳性 |
| CDS | 真正编码蛋白质的序列，用来把"长得像 TE 的基因"从 TE 库里剔除 |
| 筛选库(curated lib) | 人工校对过的 TE 序列库，作为鉴定的"已知种子" |

## 输入 | Input

- **基因组列表文件**(必需)：每行一个基因组 FASTA 路径(空格/制表符分隔的第一列即路径)

```text
/path/to/genomeA.fa
/path/to/genomeB.fa
/path/to/genomeC.fa
```

- 可选 `-c` CDS 序列文件
- 可选 `-l` 人工筛选的 TE 库

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 给基因组列表，`-o` 给输出目录(默认 `./panedta_output`)。

### 参考库 | Reference libraries

**通俗理解|In plain words:** `-c`(CDS)用来剔除"误把基因当转座子"的情况；`-l`(筛选库)提供已知的可靠 TE 种子，提高鉴定精度。**两者都是可选的，有就传，没有也能跑。**

### 行为控制 | Behavior

**通俗理解|In plain words:** `--fl-copy` 是"全长拷贝数达到多少才算真的 TE"的阈值，调高=更严格(要更多全长证据)、调低=更宽松；`--anno` 控制是否做全基因组注释(默认 1 做)；`--overwrite` 是否覆盖旧结果(默认 0 不覆盖，即已完成的基因组会跳过)。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-list` | 必填 |  | 基因组列表文件｜Genome list file |
| `-c, --cds` | — |  | CDS序列文件｜CDS sequences file |
| `-l, --curatedlib` | — |  | 筛选TE库｜Curated TE library |
| `-f, --fl-copy` | `3` | int | 全长拷贝数阈值｜Full-length copy number cutoff |
| `-a, --anno` | `1` | IntRange | 执行全基因组注释｜Perform whole-genome annotation |
| `--overwrite` | `0` | IntRange | 覆盖已有结果｜Overwrite existing results |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-o, --output-dir` | `./panedta_output` | Path | 输出目录｜Output directory |
| `--edta-path` | — | Path | EDTA安装路径｜EDTA installation path |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 交给上游 panEDTA.sh，由它完成"逐个基因组找 TE → 合并成泛 TE 库 → (可选)注释回各基因组"。

```text
基因组列表文件
    │
    ▼
验证列表中的基因组文件
    │
    ▼
构建并执行: bash panEDTA.sh -g genome_list [-c cds] [-l curatedlib] -f fl_copy -a anno -o overwrite -t threads
    │
    ▼
检查输出目录并结束
```

## 输出 | Output

```text
panedta_output/
└── (panEDTA.sh 的产物,由上游工具决定具体文件结构)
```

PanEDTA 的输出由上游 `panEDTA.sh` 生成，典型包含泛转座子库(panTE library)与各基因组的 TE 注释文件，具体文件名以运行日志中的输出目录内容为准。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 重点看两样东西：合并出来的"泛转座子库"(共有多少类 TE、每类多少拷贝)，以及各基因组里 TE 的分布注释。

- 泛 TE 库的大小反映该群体的 TE 多样性；库越精简冗余越低
- 各基因组 TE 注释可用于比较"谁的转座子多、哪类转座子扩张了"
- 若提供了 CDS，注意核对"被剔除的类基因序列"是否合理

## 参数选择建议 | Parameter Guidance

- **`--fl-copy`**：默认 3(至少 3 个全长拷贝才认定为真 TE)；基因组少或组装质量一般时可降到 1~2 以保留更多候选
- **`-c`**：有可靠的 CDS 集合(如同物种注释)就传，能明显减少"基因被误判为 TE"
- **`-l`**：有历史积累的人工筛选库就传，鉴定精度更高
- **`--overwrite 1`**：换了参数想重算时使用；默认 0 会跳过已有结果

## 依赖 | Dependencies

- EDTA 套件(含 panEDTA.sh)，conda 环境 `EDTA_v.2.2.2`（默认路径 `~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA/panEDTA.sh`）
- bash 与 EDTA 依赖的系列工具(RepeatModeler、RepeatMasker 等，随 EDTA 环境提供)

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
本包装器本身不做跳过逻辑，但 `--overwrite`(默认 0)会交给上游 panEDTA.sh 判断：已完成的基因组不覆盖。想强制重算用 `--overwrite 1`。

**Q2：找不到 panEDTA.sh？**
用 `--edta-path` 指向 EDTA 安装目录(该目录下应有 panEDTA.sh)，或确认 conda 环境 `EDTA_v.2.2.2` 已安装。

**Q3：基因组列表里文件不存在会怎样？**
验证阶段会对不存在的文件打警告，但只要有至少一个有效基因组文件就继续运行。

**Q4：CDS 是必需的吗？**
不是。CDS 和筛选库都是可选项，不传也能跑，只是鉴定精度可能略低。
