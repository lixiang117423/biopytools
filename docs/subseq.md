# subseq 序列子集提取 | FASTA Subsequence Extraction

一句话理解：**从一个大 FASTA 文件里「按需挑出」你关心的若干条序列（按 ID 清单、按名字模式、或按长度范围），导出成一个小 FASTA**。

## 功能概述 | Overview

- 三种互斥的提取方式：按 ID 列表、按名称模式匹配、按长度范围筛选
- 按 ID 列表提取时，保持列表顺序输出，并报告哪些 ID 未找到
- 模式匹配支持 contains / startswith / endswith / regex 四种，可忽略大小写
- 基于 Biopython，支持标准 FASTA（含 `.gz` 压缩）

## 快速开始 | Quick Start

```bash
biopytools subseq -i input.fasta -l id_list.txt -o output.fasta
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| FASTA ID | 名称行 `>` 后面第一段文字（到第一个空白字符为止），匹配时用的就是它 |
| ID 列表 | 一个纯文本文件，每行一个 ID，告诉你「要挑哪些」 |
| 正则(regex) | 一种用符号描述「名字长什么样」的写法，适合复杂匹配规则 |

## 输入 | Input

- 输入 FASTA：标准 FASTA 格式（Biopython 解析，`.fa` / `.fasta` / `.fa.gz` 等均可）。
- ID 列表：每行一个 ID，空行自动跳过。

ID 列表示例：

```text
seq1
seq3
seq5
```

> 注意：ID 匹配用的是 FASTA 名称行 `>` 后的**第一段**（Biopython 的 `record.id`），即到第一个空格/制表符为止。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 告诉工具「从哪个文件取、写到哪个文件」。

- `-i, --input`：输入 FASTA 文件路径。
- `-o, --output`：输出 FASTA 文件路径。

### 提取方式（三者互斥）| Extraction method

**通俗理解|In plain words:** 三种挑法只能选一种。给 `-l` 按清单挑，给 `-p` 按名字挑，都不给就按长度挑。

- `-l, --id-list`：按 ID 列表提取。
- `-p, --pattern`：按名称模式匹配提取。
- 都不给时自动进入长度筛选模式。

### 模式匹配 | Pattern matching

**通俗理解|In plain words:** 只在使用 `-p` 时生效。`--pattern-type` 决定「怎么算匹配」：contains 是名字里含这段，startswith/endswith 是开头/结尾是这段，regex 是正则。`--ignore-case` 打开后不区分大小写。

- `--pattern-type`：`contains` / `startswith` / `endswith` / `regex`（默认 `contains`）。
- `--ignore-case`：忽略大小写。

### 长度筛选 | Length filtering

**通俗理解|In plain words:** 只在长度筛选模式（不给 `-l` / `-p`）时生效。给一个长度区间，只留长度落在区间内的序列。

- `--min-length`：最小序列长度（默认 0）。
- `--max-length`：最大序列长度（默认不限）。

### 其他 | Others

**通俗理解|In plain words:** 控制输出顺序和日志位置，**一般不用动**。

- `--no-order`：不保持 ID 列表顺序（默认保持）。
- `--log-dir`：日志输出目录（默认当前目录 `.`）。

## 输出 | Output

```text
output.fasta              # 提取出的序列子集
subseq_extraction.log     # 运行日志（位于 --log-dir，默认当前目录）
```

- 输出仍为标准 FASTA，序列内容原样保留。
- 按 ID 列表模式：按列表顺序输出，日志会打印匹配到 / 未找到的 ID 数量。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看日志里「成功匹配 ID 数」和「未找到的 ID」列表，确认挑出来的对不对、有没有漏。

- 匹配到多少、缺多少，日志里都有数字；未找到的 ID 会逐条列出（超过 10 个只显示前 10 个）。
- 若一个都没匹配到，工具会报错退出，需检查 ID 是否和名称行第一段一致。
- 模式匹配和长度筛选模式会打印匹配到的序列条数。

## 参数选择建议 | Parameter Guidance

- **按清单挑**：用 `-l`，最直接、可控，还能知道漏了哪些。
- **按名字挑**：想一次挑「同一类」序列（如所有 `chr1_*`）用 `-p chr1_` + `--pattern-type startswith`。
- **复杂规则**：用 `-p` + `--pattern-type regex` 写正则（如 `^scaffoldd+$`）。
- **只留长序列**：都不给 `-l` / `-p`，配合 `--min-length` 即可。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `--output, -o` | 必填 | Path | 输出FASTA文件路径｜Output FASTA file path |
| `--id-list, -l` | — | Path | ID列表文件路径｜ID list file path |
| `--pattern, -p` | — |  | 模式匹配字符串｜Pattern matching string |
| `--pattern-type` | `contains` | contains/startswith/endswith/regex | 模式类型｜Pattern type |
| `--ignore-case` | — |  | 忽略大小写｜Case insensitive |
| `--min-length` | `0` | int | 最小序列长度｜Minimum sequence length |
| `--max-length` | — | int | 最大序列长度｜Maximum sequence length |
| `--no-order` | — |  | 不保持ID列表顺序｜Do not keep ID list order |
| `--log-dir` | `.` | Path | 日志输出目录｜Log output directory |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `-o, --output` | 必填 |  | 输出FASTA文件路径｜Output FASTA file path |
| `-l, --id-list` | — |  | ID列表文件路径｜ID list file path |
| `-p, --pattern` | — |  | 模式匹配字符串｜Pattern matching string |
| `--length-only` | — | store_true | 仅使用长度筛选｜Use only length filtering |
| `--pattern-type` | `contains` | contains/startswith/endswith/regex | 模式类型｜Pattern type (default: contains) |
| `--ignore-case` | — | store_true | 忽略大小写｜Case insensitive |
| `--min-length` | `0` | int | 最小序列长度｜Minimum sequence length (default: 0) |
| `--max-length` | — | int | 最大序列长度｜Maximum sequence length (default: unlimited) |
| `--no-order` | — | store_true | 不保持ID列表顺序｜Do not keep ID list order |
| `--log-dir` | `.` |  | 日志输出目录｜Log output directory (default: current directory) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 包 `biopython`（Biopython，用于 FASTA 解析）。

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。每次运行重新解析并覆盖输出文件。

**Q2：为什么 ID 明明在文件里却匹配不到？**
匹配用的是名称行 `>` 后的第一段（到空格/制表符为止）。若你的 ID 含空格或依赖完整描述，需先确认第一段一致，或用 `--pattern-type contains` 匹配。

**Q3：三种方式能同时用吗？**
不能。`-l`、`-p`、长度筛选三者互斥，同时给会报错。

**Q4：`--no-order` 有什么用？**
默认按 ID 列表顺序输出（便于和列表对齐）；加 `--no-order` 则按原 FASTA 中的出现顺序输出。
