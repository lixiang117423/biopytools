# Assembly 转 AGP 格式 | Assembly to AGP Converter

一句话理解：**把 ALLHiC / JBAT 等挂载流程产出的 `.assembly` 文件，翻译成数据库通用的 AGP 格式，并顺便按长度排出「前 N 条染色体」清单**，方便提交和下游工具衔接。

## 功能概述 | Overview { #overview }

- 把 `.assembly` 文本（描述 contig 在 scaffold 上的排列与方向）转成标准 9 列 AGP 格式
- 每条 contig 之间插入缺口（gap），缺口大小由 `-g` 控制，默认 100 bp
- 按 scaffold 长度从大到小排序，取前 N 条生成 `chr.list`（染色体名 + 长度），可直接喂给 ALLHiC 等下游
- 输出前缀自动去掉 `.assembly` / `.agp` 后缀，避免文件名重复
- 已存在输出文件时不覆盖，需 `--force` 才会重写

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools assembly2agp -a corrected_asm.FINAL.assembly -p output_prefix -n 12
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| .assembly 文件 | 挂载流程输出的「摆放说明」，告诉你在每条染色体上 contig 怎么排队、方向朝哪 |
| AGP | NCBI 等数据库通用的「组装描述」格式，把「哪段是序列、哪段是缺口」逐行写清楚 |
| chr.list | 染色体清单文件，两列（名字 + 长度），很多下游工具靠它知道有哪些染色体 |
| gap（缺口） | 两条 contig 之间不确定的空档，默认填 100 bp |
| 方向（orientation） | contig 在 scaffold 上是正着（+）还是反着（-）放 |

## 输入 | Input { #input }

`.assembly` 文本文件，两种行：

- 以 `>` 开头：contig 定义行，格式为 `>contig名称 编号 长度 ...`（空格分隔，第 2 个字段是编号、第 3 个字段是长度）
- 普通行：空格分隔的一串 contig 编号，一个普通行 = 一条 scaffold；编号前带负号表示该 contig 反向放置

```text
>ptg000001l 1 190127
>ptg000002l 2 5326271
1 -2
2 1
```

上面示例表示：scaffold 1 由 contig 1（正向）和 contig 2（反向）拼成，scaffold 2 由 contig 2（正向）和 contig 1（正向）拼成。`#` 开头的行会被忽略。

## 参数说明 | Parameters { #parameters }

### 必需输入 | Required

**通俗理解|In plain words:** `-a` 是输入的 `.assembly` 文件；`-p` 是输出文件名的前缀；`-n` 告诉工具「最终染色体取前几条」，取值要大于 0，通常就是目标染色体数（如 12、21）。

相关参数：`-a, --assembly`、`-p, --prefix`、`-n, --num-chromosomes`。

### 输出与缺口 | Output and gap

**通俗理解|In plain words:** `-o` 指定输出目录，默认当前目录；`-g` 是 contig 之间缺口的长度，默认 100 bp，一般不用改；`-f, --force` 用于覆盖已存在的输出文件（默认检测到同名文件会报错退出，防止误覆盖）。

相关参数：`-o, --output-dir`、`-g, --gap`、`-f, --force`。

### 日志控制 | Logging

**通俗理解|In plain words:** 这些只影响屏幕输出的详略，不影响结果文件。`-v` 加一次显示 INFO、加两次显示 DEBUG；`--quiet` 只显示错误；`--log-file` 把完整日志写到指定文件。**一般不用动。**

相关参数：`-v, --verbose`、`--quiet`、`--log-file`。

## 输出 | Output { #output }

```text
输出目录/
├── {prefix}.agp       # 标准 9 列 AGP 文件
└── {prefix}.chr.list  # 染色体清单：scaffold名 <TAB> 长度（按长度降序，取前 N 条）
```

- AGP 文件 9 列：Chromosome、Start、End、Order、Tag（W=contig / U=gap）、Contig_ID、Contig_start、Contig_end、Orientation（+/-）。每条 contig（W）后跟一个缺口（U）行，缺口长度即 `-g` 的值
- chr.list 每行两列：染色体名和总长，按长度从大到小，是下游挂载/可视化工具的常见输入

## 结果解读 | Interpreting Results { #interpreting-results }

- `chr.list` 第一行就是最长的那条染色体，长度约等于该染色体的总长
- AGP 里 `U` 行越多，说明这条 scaffold 被拆得越碎（缺口多），连续程度越低
- `Orientation` 为 `-` 的 contig 表示它在组装时被反向放置，属正常现象，不代表错误

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- `-n` 按物种实际染色体数填写；不确定时可先给个大一点的数，再按 `chr.list` 的长度分布截取
- 重跑同一前缀且不想被「文件已存在」拦住时，加 `-f`
- 需要记录过程时加 `--log-file run.log`，平时靠 `-v` 即可

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--assembly, -a` | 必填 |  | 组装文件｜Assembly file path |
| `--prefix, -p` | 必填 |  | 输出前缀｜Output prefix |
| `--output-dir, -o` | `.` | Path | 输出目录｜Output directory path |
| `--gap, -g` | `100` | int | Scaffold间隙(bp)｜Scaffold gap size in bp |
| `--num-chromosomes, -n` | 必填 | int | 染色体数量｜Number of chromosomes |
| `--force, -f` | — |  | 强制覆盖｜Force overwrite existing files |
| `--verbose, -v` | — |  | 详细模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅ERROR)｜Quiet mode (ERROR only) |
| `--log-file` | — | Path | 日志文件｜Log file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-a, --assembly` | 必填 |  | Input assembly file path |
| `-p, --prefix` | 必填 |  | Output prefix (will be used for both AGP and chr.list files) |
| `-o, --output-dir` | `.` |  | Output directory path (default: current directory) |
| `-g, --gap` | `100` | int | Gap size between scaffolds in bp (default: 100) |
| `-n, --num-chromosomes` | 必填 | int | Number of chromosomes to include in chr.list file |
| `-v, --verbose` | `0` | count | Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | Quiet mode (only ERROR) |
| `--log-file` | — |  | Log file path |
| `-f, --force` | — | store_true | Force overwrite existing files |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- pandas（DataFrame 处理）
- 无 conda 环境、无外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：重跑时报「output file already exists」？**
这是防误覆盖保护。确认要覆盖时加 `-f, --force`，或换一个 `-p` 前缀 / `-o` 目录。

**Q2：`chr.list` 里的名字为什么都是 scaffold_N？**
程序会按长度排序后把前 N 条重命名为 `scaffold_1`、`scaffold_2`……，方便下游统一识别，不保留原 contig 名。

**Q3：会断点续传吗？**
不会。转换是一次性完成的，没有中间步骤可跳过。
