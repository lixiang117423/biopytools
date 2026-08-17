# 基因组 ID 顺序重命名 | Rename Genome ID

一句话理解：**把 FASTA 里的所有序列按出现顺序重新编号命名（默认 Chr01、Chr02…）**，可顺便只保留前 N 条并生成一张「旧 ID → 新 ID」对照表。
输入一个 FASTA，输出重命名后的 FASTA 和一份 ID 映射文件。

## 功能概述 | Overview

- 按顺序把每条序列重命名为 `Chr01`、`Chr02`…（前缀、零填充宽度均可自定义）
- 可用 `-n/--chr-count` 只输出前 N 条作为染色体，其余丢弃
- 默认生成 ID 映射文件（旧 ID → 新 ID → 原描述），可用 `--no-mapping` 关闭
- 用 Biopython 解析，兼容各种 FASTA 变体；纯 Python 实现，无需外部命令行软件

## 快速开始 | Quick Start

```bash
biopytools rename-genome-id -i input.fa -o output.fa -n 20
```

最小输入：一个 FASTA + 输出文件。`-n 20` 表示只保留前 20 条并重命名为 Chr01~Chr20。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 序列 ID | FASTA 里 `>` 号后面的名字 |
| 前缀(prefix) | 新名字开头的固定部分，默认 `Chr` |
| 零填充 | 数字前补 0，如 Chr01 而不是 Chr1，保证字典序不乱 |
| ID 映射文件 | 一张「改名前 ↔ 改名后」对照表，方便回查原来的名字 |

## 输入 | Input

标准 FASTA 文件。重命名按**文件里的出现顺序**进行：第 1 条 → Chr01、第 2 条 → Chr02…。因此输入前请先按染色体从大到小排序，让 Chr01 对应最大的染色体。

## 参数说明 | Parameters

### 重命名规则 | Renaming rules

**通俗理解|In plain words:** `-p/--prefix` 是新名字的前缀（默认 `Chr`）；`-w/--padding-width` 是数字补几位 0（默认 2，即 01、02…，超过 99 条会自动变 3 位）；`--no-zero-padding` 则不做补零（Chr1、Chr2…）。这些一般用默认即可。

### 染色体提取 | Chromosome extraction

**通俗理解|In plain words:** `-n/--chr-count` 表示「只保留前 N 条当染色体」，`0` 表示全部保留并重命名。当你只想把前 20 条变成染色体、丢掉后面小 scaffold 时用它。

### ID 映射 | ID mapping

**通俗理解|In plain words:** 默认会生成映射文件（`输出文件名_stem_id_mapping.txt`），记录每条序列的旧 ID、新 ID、原描述。`--no-mapping` 可关闭；`--mapping-file` 可自定义映射文件路径。

## 输出 | Output

```text
输出目录/
├── output.fa                    # 重命名后的 FASTA(-o 指定)
├── output_id_mapping.txt        # 旧ID → 新ID → 描述 对照表(默认生成)
└── fasta_id_renamer.log         # 运行日志
```

- `output.fa`：重命名后的序列，序列头为新 ID（原描述不保留）
- `output_id_mapping.txt`：三列 TSV（Original_ID / New_ID / Description），改名前后的对应关系都在这里

## 结果解读 | Interpreting Results

- 用 `grep "^>" output.fa` 检查序列头是否为整齐的 Chr01、Chr02…
- 序列头只保留新 ID，**原描述信息被丢弃**，但完整保存在映射文件里，需要回查时看映射文件
- 若只输出了前 N 条，其余序列被丢弃（不写入任何文件），确认 `-n` 是否符合预期

## 参数选择建议 | Parameter Guidance

- **只要染色体、丢弃 scaffold**：`-n 染色体数`
- **全部重命名、一条不丢**：不写 `-n`（默认 0）
- **不要 Chr 前缀**：`-p Scaffold` 或 `-p contig`
- **想要 Chr1 而非 Chr01**：`--no-zero-padding`
- **输入顺序最重要**：先按长度排序再输入，让编号和大小对应

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件｜Input FASTA file |
| `--output, -o` | 必填 |  | 输出FASTA文件｜Output FASTA file |
| `--prefix, -p` | `Chr` |  | 序列前缀｜Sequence prefix |
| `--no-zero-padding` | — |  | 不使用零填充(如Chr1而非Chr01)｜Do not use zero padding |
| `--padding-width, -w` | `2` | int | 填充宽度｜Padding width |
| `--chr-count, -n` | `0` | int | 染色体数量，提取前N条作为染色体｜Chromosome count, extract first N as chromosomes |
| `--no-mapping` | — |  | 不保存ID映射文件｜Do not save ID mapping file |
| `--mapping-file` | — |  | ID映射文件路径｜ID mapping file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 包 `biopython`（`Bio.SeqIO`，用于解析 FASTA）

## 常见问题 | FAQ

**Q1：和 `biopytools rename-chromosomes` 有什么区别？**
两者都重命名序列，但侧重点不同：本工具**逐条顺序编号**（前缀/补零宽度可定制、能生成 ID 映射文件、可只提取前 N 条）；`rename-chromosomes` 固定输出 Chr01 格式、额外支持把其余序列改名为 HiC_scaffold。需要保留旧 ID 对照就用本工具。

**Q2：原序列描述会丢吗？**
输出 FASTA 里只保留新 ID（原描述被新 ID 覆盖），但原描述完整保存在 `id_mapping.txt` 里。

**Q3：有断点续传吗？**
没有。本工具单步流式处理，重跑即重新生成，无需续传。

**Q4：映射文件不想要能关吗？**
能，加 `--no-mapping`；或想换位置用 `--mapping-file 路径`。
