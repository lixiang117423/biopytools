# needle_identity - 序列两两一致性计算 | Pairwise Sequence Identity (EMBOSS needle)

一句话理解：**把输入 FASTA 里的序列两两比对一遍，算出每对序列「长得有多像」的百分比表，用来判断哪些序列是重复/近亲、亲缘关系有多近。**

## 功能概述 | Overview { #overview }

- 用 EMBOSS needle 对输入 FASTA 中所有序列做两两全局比对（共 n(n-1)/2 对）
- 输出一张 TSV 表：每行一对序列的 identity、similarity、比对长度、空位、得分
- 多线程并行（ThreadPoolExecutor），线程数可调（默认 12）
- 自动检测并 conda run 调用 needle，无需手动指定 conda 环境
- 生成 software_versions.yml 与运行日志

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools needle-identity -i sequences.fa -o output_dir
```

最小输入：一个含至少 2 条序列的 FASTA 文件（序列 ID 不能重复）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| FASTA | 用纯文本存序列的通用格式，> 开头一行是名字，下面是序列字母 |
| 比对(alignment) | 把两条序列并排对齐（允许插入空位）看哪些位置字母相同，像把两段文字对齐找相同字 |
| identity(一致性) | 对齐后「相同位置字母相同」的比例(%)，本工具的核心结果 |
| similarity(相似性) | 相似氨基酸的比例(%)，蛋白质比对时比 identity 更宽松、更有参考意义 |
| gap(空位) | 为对齐而在某条序列里插入的空位，对应字母的缺失/插入 |
| gapopen / gapextend | 开一个空位 / 延长已有空位的扣分，数值越大比对越「舍不得」开空位 |
| 罚分(penalty) | 比对软件给「不匹配/空位」打的分，越低越好 |

## 输入 | Input { #input }

- FASTA 文件(-i)：至少 2 条序列；序列 ID（> 后第一个词）必须唯一，否则直接报错。

示例：

```text
>seqA
ATGGCTAGCTAGCTAGCTA
>seqB
ATGGCTAGCGAGCTAGCTA
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** -i 是输入 FASTA，-o 是输出目录（默认 ./output）。输出文件名为「输入文件名（去扩展名）.needle_identity.tsv」，日志在 99_logs/ 下。这两个必填，其余都有合理默认。

### 比对罚分 | Gap penalties

**通俗理解|In plain words:** --gapopen（默认 10.0）与 --gapextend（默认 0.5）直接透传给 needle，控制开空位/延长空位的扣分。调大 gapopen 会让比对更「怕」开空位（更紧凑、少空位），调小则更容忍空位。绝大多数场景用默认即可，不需要动；除非序列空位特别多、结果不理想，才考虑微调。

### 运行参数 | Runtime

**通俗理解|In plain words:** --threads（默认 12）控制并行比对的线程数，序列对越多提速越明显；--needle-path 指定 needle 可执行文件路径，不指定时按「环境变量 NEEDLE_PATH 到 配置文件 到 默认 ~/miniforge3/envs/protein/bin/needle」的优先级自动查找。一般不用动。

## 分析流程 | Pipeline { #pipeline }

```text
解析 FASTA（校验至少 2 条、ID 不重复）
    |
    v
每条序列写入临时文件 tmp/
    |
    v
两两组合，并行调用 needle 比对
    |
    v
解析 needle 报告(Length/Identity/Similarity/Gaps/Score)
    |
    v
写结果表 {stem}.needle_identity.tsv 与版本信息、日志
    |
    v
清理 tmp/
```

## 输出 | Output { #output }

```text
output_dir/
|-- 00_pipeline_info/
|   -- software_versions.yml               # EMBOSS 版本与参数存档
|-- {输入文件名去扩展名}.needle_identity.tsv   # 两两一致性结果表（核心）
-- 99_logs/
    -- {输入文件名去扩展名}.needle_identity.log
```

结果表（TSV，制表符分隔）列：

```text
seq1  seq2  identity_percent  matches  aligned_length  gaps  similarity_percent  score
```

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 打开 TSV，每一行回答「这两条序列有多像」。核心看 identity_percent：100 表示完全相同（可能是重复序列），越接近 100 越像。

- seq1 / seq2：参与比对的序列 ID
- identity_percent：一致性百分比，最常用列，数值越高两条序列越接近
- matches：对齐后完全相同的字母数
- aligned_length：比对总长度（含空位）
- gaps：空位数量；空位很多说明两条序列长度/结构差异大
- similarity_percent：相似性百分比（氨基酸的保守替换也算相似），蛋白质比对时比 identity 更有意义
- score：needle 原始比对得分，负数越接近 0 表示越好

好坏判据：核酸序列 identity 达到 95% 左右通常视为「几乎相同/近源」；蛋白质相似性达到 60-70% 往往提示同源。阈值因物种和目的而异，需结合业务判断。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规两两比对：全部默认即可（threads=12、gapopen=10.0、gapextend=0.5）
- 序列条数很多（大于 50）：两两组合数 = n(n-1)/2，会快速增长；适当调大 --threads 提速，注意 needle 单对比对也吃内存
- 空位很多导致比对难看：可尝试调小 --gapextend（让延伸空位更便宜、更容忍插入缺失）
- 自定义 needle：用 --needle-path 或环境变量 NEEDLE_PATH 指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA序列文件｜Input FASTA file |
| `--output-dir, -o` | `./output` |  | 输出目录｜Output directory (default: ./output) |
| `--needle-path` | — |  | needle可执行文件路径｜needle executable path |
| `--threads` | `12` | int | 并行线程数｜Threads (default: 12) |
| `--gapopen` | `10.0` | float | gap开放罚分｜Gap open penalty (default: 10.0) |
| `--gapextend` | `0.5` | float | gap延伸罚分｜Gap extend penalty (default: 0.5) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA序列文件｜Input FASTA file |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory (default: ./output) |
| `--needle-path` | — |  | needle可执行文件路径｜needle executable path (default: ~/miniforge3/envs/protein/bin/needle) |
| `--threads` | `12` | int | 并行线程数｜Threads (default: 12) |
| `--gapopen` | `10.0` | float | gap开放罚分｜Gap open penalty (default: 10.0) |
| `--gapextend` | `0.5` | float | gap延伸罚分｜Gap extend penalty (default: 0.5) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- EMBOSS 套件（命令行 needle），默认走 conda 环境 protein（~/miniforge3/envs/protein/bin/needle），通过 conda run 自动检测并调用
- Biopython（仅用于解析输入 FASTA）
- Python 3 与 PyYAML（写版本信息）

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持。每次重跑都会重新计算所有序列对（临时目录 tmp/ 运行结束即清理）。序列很多时建议先拿几条试跑确认参数再全量。

Q2：为什么报「至少需要 2 条序列」或「序列 ID 重复」？
工具要做「两两」比对，至少 2 条；且每条序列的 ID（> 后第一个词）必须唯一，重复 ID 无法区分结果，会直接报错。

Q3：needle 报不可用怎么办？
默认按「NEEDLE_PATH 环境变量、~/.config/biopytools/config.yml、~/miniforge3/envs/protein/bin/needle」顺序查找。找不到时用 --needle-path 显式指定你机器上的 needle 路径。

Q4：identity 和 similarity 有什么区别？
identity 只算「字母完全相同」；similarity 额外把「性质相近的氨基酸」也计为相似（仅蛋白质有意义）。两条蛋白质可能 identity 低但 similarity 高，说明同源但变异较大。

Q5：gapopen / gapextend 改大会怎样？
gapopen 变大，比对更不愿意开新空位（空位变少、比对更紧凑）；gapextend 变大，更不愿意把空位拉长。通常只在默认结果明显不合理时才调。