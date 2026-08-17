# split-fasta-id FASTA ID 分割 | Split FASTA Header IDs

一句话理解：**把 FASTA 里「又长又乱」的序列名称行（header）按分隔符拆开，只保留你想要的第 N 段，让序列 ID 变短变干净**，方便下游比对和结果展示。

## 功能概述 | Overview

- 从每条序列名称行中提取指定位置（0 起始）的一个元素作为新 ID
- 分隔符可自动检测，也支持空格、制表符、两者混合或任意自定义字符（如 `,`、`|`）
- 可选保留原文件备份、是否跳过空名称行
- 单遍流式处理，纯 Python 实现，无第三方软件依赖

## 快速开始 | Quick Start

```bash
biopytools split-fasta-id -i input.fasta -o output.fasta -p 0
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| FASTA header | 每条序列开头以 `>` 开始的那一行，通常包含 ID 和一堆描述信息 |
| 分隔符 | 把 header 里「一长串文字」切成几段的刀口，常见是空格或制表符 |
| 位置(position) | 切成几段后，取第几段；从 0 数起，0 是第一段 |

## 输入 | Input

- 标准 FASTA 文件（扩展名 `.fasta` / `.fa` / `.fas`）。
- 名称行格式示例：`>chr1 gene=xxx description here`。按空格切分会得到 `chr1`、`gene=xxx`、`description`、`here` 四段。

```text
>chr1 gene=xxx description here
ACGTACGT
>chr2 another example
TTTTGGGG
```

## 参数说明 | Parameters

### 分割参数 | Splitting

**通俗理解|In plain words:** 决定「按什么切、切完取哪一段」，是本工具的核心。位置越界时会自动回退到第一段，所以填错也通常不会报错、只会取到第一段。

- `-p, --position`：提取位置，0 起始（默认 0，取第一段）。
- `-d, --delimiter`：分隔符（默认 `auto`）。可选 `auto`（自动检测）、`space`（空格）、`tab`（制表符）、`both`（空格和制表符都算），或任意单个字符如 `,`、`|`。

### 处理选项 | Processing options

**通俗理解|In plain words:** 控制几个「边角情况」怎么处理，**一般不用动**。

- `--keep-original`：把原文件复制一份 `<输入文件>.backup` 再改。
- `--no-skip-empty`：默认跳过空名称行，加此选项后不跳过（保留空名称行）。
- `--preserve-comments`：声明为「保留名称行中的注释」，当前实现中尚未生效（见 FAQ）。

### 输出 | Output

**通俗理解|In plain words:** 指定写到哪里。

- `-o, --output`：输出 FASTA 路径（默认 `output.fasta`）。

## 输出 | Output

```text
output.fasta            # 新 ID 的 FASTA（序列内容不变）
fasta_split.log         # 运行日志（与输出文件同目录）
input.fasta.backup      # 仅在加 --keep-original 时生成
```

- 序列字母内容原样保留，只改名称行。
- 日志里会打印成功处理 / 跳过的序列条数。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开输出文件，看每条 `>` 后面是不是只剩你想要的那一段、有没有重复 ID。

- 分割正常时，`>chr1 gene=xxx` 取位置 0 会变成 `>chr1`。
- 若某条 header 无法切出指定位置，会回退到第一段；若整行都处理不了，会保留原名称行并在日志里记 warning。
- 多条序列可能被切成同一个 ID（本工具不去重），后续若用 ID 取序列要注意这一点。

## 参数选择建议 | Parameter Guidance

- **`-p`**：通常取 0（第一段，一般是 ID）。NCBI 风格的 `>gi|123|ref|...` 想取登录号可结合 `-d |` 与相应位置。
- **`-d`**：拿不准就用默认 `auto`（按前 10 条名称行自动判断空格还是制表符）。名称里既可能有空格又有制表符时，用 `both`。
- **`--keep-original`**：第一次跑、怕改坏时加上，确认无误后再删掉备份。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `--output, -o` | `output.fasta` | Path | 输出FASTA文件路径｜Output FASTA file path |
| `--position, -p` | `0` | int | 提取位置(0表示第一个元素)｜Extract position (0 means first element) |
| `--delimiter, -d` | `auto` | str | 分隔符类型｜Delimiter type: "auto"(auto detect), "space", "tab", "both"(space and tab), or any character like "," or "｜" |
| `--keep-original` | — |  | 保留原始文件作为备份｜Keep original file as backup |
| `--no-skip-empty` | — |  | 不跳过空的序列名称行｜Do not skip empty sequence name lines |
| `--preserve-comments` | — |  | 保留序列名称行中的注释｜Preserve comments in sequence name lines |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `-o, --output` | `output.fasta` |  | 输出FASTA文件路径｜Output FASTA file path |
| `-p, --position` | `0` | int | 提取位置(0表示第一个元素)｜Extract position (0 means first element) |
| `-d, --delimiter` | `auto` |  | 分隔符类型｜Delimiter type: "auto"(auto detect), "space", "tab", "both"(space and tab), or any character like "," or "｜" |
| `--keep-original` | — | store_true | 保留原始文件作为备份｜Keep original file as backup |
| `--no-skip-empty` | — | store_true | 不跳过空的序列名称行｜Do not skip empty sequence name lines |
| `--preserve-comments` | — | store_true | 保留序列名称行中的注释｜Preserve comments in sequence name lines |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3 标准库，无第三方软件、无 conda 环境依赖。

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。每次运行重新处理整个文件并覆盖输出。

**Q2：`--preserve-comments` 为什么没有效果？**
当前代码中该参数只被声明、未在分割逻辑中使用，属于已知未实现项。需要保留注释时，请先备份原文件或改用其他方式。

**Q3：位置填大了会怎样？**
不会报错，会自动回退到第一段。若担心取错，跑完检查一下输出 header。

**Q4：`auto` 是怎么判断分隔符的？**
取前 10 条名称行，统计空格和制表符出现次数，制表符多就用 tab，有空格就用 space，都没有则默认按空格。
