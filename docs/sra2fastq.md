# sra2fastq SRA 转 FASTQ | SRA to FASTQ Conversion

一句话理解：**把从 NCBI SRA 数据库下载的 `.sra` 归档解包成能直接用于下游分析的 FASTQ 测序文件**，支持多线程加速和批量处理。

## 功能概述 | Overview

- 优先用 `parallel-fastq-dump` 多线程高速转换，缺失时自动降级到 `fastq-dump`（单线程）
- 支持单个 `.sra` 文件或整个目录批量转换
- 双端数据自动拆分（`_1` / `_2`），输出可选 gzip 压缩
- 断点续传：已转换的 SRA 自动跳过（见 FAQ）
- 可选最小读长过滤、剪切 adapter
- 输出转换摘要报告与日志

## 快速开始 | Quick Start

```bash
biopytools sra2fastq -i sra_dir/ -o fastq_output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| SRA | NCBI 的测序数据归档格式，把原始数据打包压缩成 `.sra` 文件，不能直接用于比对 |
| FASTQ | 测序读段的标准文本格式，含序列和每个碱基的质量，是下游分析的标准输入 |
| 单端/双端 | 只测一端叫单端；从片段两端各测一次叫双端，会得到 `_1` 和 `_2` 两个文件 |
| gzip 压缩 | 把大文件压成 `.gz`，省空间，下游工具都能直接读 |

## 输入 | Input

- 单个 `.sra` 文件，或含 `.sra` 文件的目录。
- 目录下若无 `.sra`，会退而取该目录下所有非隐藏文件。
- 单文件可不带 `.sra` 后缀。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 告诉工具「读哪些 SRA、写到哪」。

- `-i, --input`：输入 SRA 文件或文件夹路径。

### 转换与性能 | Conversion & performance

**通俗理解|In plain words:** 控制转换用多少资源、临时文件放哪。线程越多越快；`--tmpdir` 指定高速存储可显著加速（SRA 转换会产生大量临时数据）。

- `-o, --output`：输出目录（默认 `./fastq_output`）。
- `-t, --threads`：线程数（默认 12）。
- `--tmpdir`：临时目录（可设为高速存储以加速）。

### 输出格式 | Output format

**通俗理解|In plain words:** 控制出来的 FASTQ 长什么样。压缩省空间（下游都支持 `.gz`）；拆分决定双端数据是拆成 `_1`/`_2` 两个文件还是混在一起。**一般保持默认（都开启）即可。**

- `--compress / --no-compress`：是否压缩输出为 `.gz`（默认压缩）。
- `--split / --no-split`：是否拆分双端读段（默认拆分）。

### 过滤 | Filtering

**通俗理解|In plain words:** 在转换的同时做轻量清洗。`--min-len` 丢掉太短的读段；`--clip` 去掉 adapter 序列。**一般不用动。**

- `--min-len`：最小读长过滤（默认 0，不过滤）。
- `--clip`：剪切 adapters。

## 输出 | Output

```text
fastq_output/
├── SRR123_1.fq.gz           # 双端 Read1（压缩时）
├── SRR123_2.fq.gz           # 双端 Read2
├── SRR123.fq.gz             # 单端数据（无 _1/_2）
├── conversion_summary.txt   # 转换摘要报告
└── sra_conversion.log       # 运行日志
```

- 文件命名以 SRA 登录号（如 `SRR123`）为基础，双端拆分为 `_1` / `_2`，统一用 `.fq.gz` 后缀（转换出的 `.fastq.gz` 会被自动重命名为 `.fq.gz`）。
- `conversion_summary.txt` 里汇总了每个文件成功/跳过/失败情况及耗时。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 主要看摘要报告里「成功 / 跳过 / 失败」三行数字，失败的文件单独列出。

- 成功后输出目录出现对应的 `.fq.gz` 文件，可直接喂给 fastp、比对等下游工具。
- 跳过（already done）说明之前已经转过，本次未重复转换。
- 若某个 SRA 失败，报告会列出文件名，日志里有具体报错（常见是网络/磁盘/工具缺失）。

## 参数选择建议 | Parameter Guidance

- **`-t`**：机器核数允许时可调大（如 16、24）加速；超算上别超过分配核数。
- **`--tmpdir`**：磁盘写瓶颈明显时，把临时目录指到高速盘（如本地 NVMe）能明显提速。
- **`--no-compress`**：下游工具若明确只吃未压缩 FASTQ、又不想二次解压时再用，一般保持压缩。
- **`--min-len`**：想顺手丢掉短读段时设个阈值（如 30），否则保持 0。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入SRA文件或文件夹路径｜Input SRA file or folder path |
| `--output, -o` | `./fastq_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--tmpdir` | — | Path | 临时目录(用于加速)｜Temporary directory (for acceleration) |
| `--compress/--no-compress` | `True` |  | 压缩输出为.gz格式｜Compress output to .gz format |
| `--split/--no-split` | `True` |  | 拆分双端测序文件｜Split paired-end reads |
| `--min-len` | `0` | int | 最小读长过滤｜Minimum read length filter |
| `--clip` | — |  | 剪切adapters｜Clip adapters |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入SRA文件或文件夹路径｜Input SRA file or folder path |
| `-o, --output` | `./fastq_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-c, --compress` | `True` | store_true | 压缩输出为.gz格式｜Compress output to .gz format |
| `--no-compress` | — | store_false | 不压缩输出｜Do not compress output |
| `--split` | `True` | store_true | 拆分双端测序文件｜Split paired-end reads |
| `--no-split` | — | store_false | 不拆分文件｜Do not split files |
| `--tmpdir` | — | str | 临时目录 (可以设置为高速存储以加速)｜Temporary directory (can be set to fast storage for acceleration) |
| `--min-len` | `0` | int | 最小读长过滤｜Minimum read length filter |
| `--clip` | — | store_true | 剪切adapters｜Clip adapters |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `parallel-fastq-dump`（推荐，多线程）或 `fastq-dump`（SRA Toolkit，备选），两者至少其一
- `fastq-dump` 从 **misc 域环境**解析（sra-tools 3.4.1）；`parallel-fastq-dump` 目前仅安装于旧环境 `sratoolkit_v.2.5.7`，默认路径锁定到该环境（避免跨环境扫描漂移），将来装入 misc 后自动切换域环境
- 可用环境变量 `PARALLEL_FASTQ_DUMP_PATH` / `FASTQ_DUMP_PATH` 或 `~/.config/biopytools/config.yml` 覆盖
- 程序启动时自动检测：有 `parallel-fastq-dump` 优先用它，否则回退 `fastq-dump`（misc 域环境版本）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。转换前检查输出目录是否已有该 SRA 对应的 FASTQ（`{base}.fq.gz` / `{base}_1.fq.gz` / `{base}_2.fq.gz` 及 `_pass` 变体；未压缩的 `.fastq`/`.fq` 同样计入），已存在则跳过。按**精确文件名**匹配，前缀相近的样本（如 `SRR12` 与 `SRR123`）不会互相误判。

**Q2：为什么输出是 `.fq.gz` 而不是 `.fastq.gz`？**
工具在转换完成后会把 `.fastq.gz` 统一重命名为 `.fq.gz`，保持输出命名一致；两者内容完全一样。

**Q3：没装 parallel-fastq-dump 能用吗？**
能，会自动降级到 `fastq-dump`（单线程，速度慢），并在日志里给出提示。两个都没有则报错退出。

**Q4：转换时总是跳过技术读段吗？**
是的，`--skip-technical` 在内部固定开启（未暴露到 CLI），会跳过技术序列（如接头引物）。

**Q5：目录里同时有 `.sra` 和其他文件会怎样？**
优先只处理 `.sra` 文件；目录里没有 `.sra` 时才会把其他非隐藏文件当作输入。

**Q6：部分文件转换失败时程序的退出码是？**
非零（1）。只要有文件失败，即使其余全部成功，程序也以退出码 1 结束——方便 HPC 上游脚本判断；全部成功（含跳过）时退出码 0。
