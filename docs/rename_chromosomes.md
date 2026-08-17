# 染色体重命名 | Rename Chromosomes

一句话理解：**把 FASTA 里五花八门的序列名统一改成 Chr01、Chr02 这种标准染色体命名**，让下游工具（挂载、绘图、注释）能正确识别染色体。
输入一个 FASTA，指定染色体条数 n，输出前 n 条已改名的序列（可选保留其余序列）。

## 功能概述 | Overview

- 把 FASTA 前 n 条序列统一重命名为 `Chr01`、`Chr02`…`ChrNN`（固定两位零填充）
- 默认只输出前 n 条；加 `--keep-all` 时，第 n 条之后改名为 `HiC_scaffold_01` 起并一并输出
- 纯 awk 流式处理，不占内存，适合超大 FASTA
- 不带参数映射文件，重命名后原序列描述信息不保留（见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools rename-chromosomes -i genome.fa -o renamed.fa -n 20
```

最小输入：一个 FASTA + 染色体条数 `-n`，输出前 20 条改名为 Chr01~Chr20 的序列。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 染色体(chromosome) | 基因组的「大件」，组装最终想得到的顶层单元 |
| 序列 ID | FASTA 里 `>` 号后面那段名字，是本工具要改的东西 |
| 零填充 | 数字前面补 0，让 Chr1 写成 Chr01，排序时不会乱 |
| scaffold | 尚未定到染色体级别的大片段，用 HiC_scaffold 命名标记 |

## 输入 | Input

标准 FASTA 文件。序列顺序即编号顺序：第 1 条改名为 Chr01、第 2 条 Chr02，依此类推，所以**输入前请先按染色体从大到小排好序**，否则最大的染色体可能被排到小编号。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-n/--number` 是「你有几条染色体」，程序只处理前 n 条。定错了要么丢染色体、要么把 scaffold 当染色体编号。

### 保留其余序列 | Keep the rest

**通俗理解|In plain words:** 默认只输出前 n 条（把 scaffold 丢掉）。若你想「染色体 + 未定位 scaffold 一起保留」，加 `--keep-all`，其余序列会被改成 `HiC_scaffold_01` 起并接在后面输出。

## 输出 | Output

- 输出 FASTA：`-o` 指定的文件，前 n 条命名为 Chr01~ChrNN（`--keep-all` 时其余为 HiC_scaffold_01 起）
- 日志文件：输出目录下 `rename_chromosomes_YYYYMMDD_HHMMSS.log`

## 结果解读 | Interpreting Results

- 用 `grep "^>" renamed.fa` 看序列头，应看到整齐的 Chr01、Chr02…
- 序列头里**只有新 ID，不含原描述**（awk 只保留新名字），这是预期行为
- 若输出的序列条数比预期少，检查 `-n` 是否给对，或是否忘了 `--keep-all`

## 参数选择建议 | Parameter Guidance

- **只要染色体**：`-n 染色体数`，不加其它开关
- **染色体 + scaffold 都要**：加 `--keep-all`
- **输入顺序很关键**：先把最大的序列排前面再跑，保证 Chr01 对应最大的染色体

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `--output, -o` | 必填 | Path | 输出FASTA文件路径｜Output FASTA file path |
| `--number, -n` | 必填 | int | 染色体数量｜Number of chromosomes |
| `--keep-all` | `False` |  | 输出所有序列（默认只输出前n个序列）｜Output all sequences (default: only first n) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `-o, --output` | 必填 |  | 输出FASTA文件路径｜Output FASTA file path |
| `-n, --number` | 必填 | int | 染色体数量｜Number of chromosomes |
| `--keep-all` | — | store_true | 输出所有序列（默认只输出前n个序列）｜Output all sequences (default: only first n) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- 系统自带的 `awk`（无需 conda 环境、无需额外软件）

## 常见问题 | FAQ

**Q1：原序列的描述信息去哪了？**
重命名后序列头只保留 `ChrNN` 这种新 ID，原描述（如 `chromosome 1`）不保留。若需要保留原始 ID 对照，请改用 `biopytools rename-genome-id`（会生成新旧 ID 映射文件）。

**Q2：为什么编号和大小对不上？**
编号按**输入顺序**来，不是按大小自动排序。请先按序列长度从大到小排序再输入。

**Q3：想保留未挂载的 scaffold 怎么办？**
加 `--keep-all`，第 n 条之后会改成 HiC_scaffold_01 起并输出。

**Q4：有断点续传吗？**
没有。本工具是单步流式处理，重跑即重新生成，很快无需续传。
