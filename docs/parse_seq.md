# parse_seq - 按区间提取序列 | Sequence Extraction by Region

一句话理解：**给一个 FASTA（基因组或蛋白库）+ 一份「哪条序列的哪一段」的区间清单，把每一段序列抠出来存成新的 FASTA 文件，还能顺手做反向互补、翻译成蛋白。**

## 功能概述 | Overview { #overview }

- 按 BED 风格区间文件从 FASTA 中批量提取子序列
- DNA 用 samtools faidx 提取（自动建 .fai 索引），蛋白用 Biopython/基础方法提取
- 支持负链自动反向互补、DNA 翻译为蛋白、合并或分开输出
- 1-based 闭区间坐标，第四列可选 +/- 指定链方向
- 分批提取（每批 100 个区间），避免命令行过长

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools parse-seq -s genome.fasta -r regions.bed -o extracted.fasta
```

最小输入：一个序列 FASTA + 一个区间文件（至少 3 列：染色体 起始 终止）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| FASTA | 纯文本存序列的格式，> 开头是名字，下面是序列字母 |
| 区间(region) | 「哪条染色体、从第几个碱基到第几个碱基」的一段坐标，像门牌号范围 |
| 1-based 闭区间 | 坐标从 1 开始数，起点终点都包含在内（与常见 BED 的 0-based 半开不同） |
| 链(strand) | 序列读的方向，+ 正链 / - 负链；负链要取反向互补 |
| 反向互补(reverse complement) | 把 DNA 倒过来并逐字母配对（A-T、G-C），得到负链方向的实际序列 |
| 翻译(translate) | 按遗传密码把 DNA 三联密码子换成氨基酸，得到蛋白序列 |
| .fai 索引 | samtools 为 FASTA 建的「目录」，让按区间随机取序列更快 |

## 输入 | Input { #input }

- 序列文件(-s)：FASTA 格式，DNA 或蛋白（由 --type 指定）。
- 区间文件(-r)：每行至少 3 列「染色体 起始 终止」，制表符或空格分隔，可用 # 开头注释；第 4 列若为 + 或 - 表示链方向（默认 +）。坐标 1-based、起始不大于终止、均为正整数。

示例：

```text
chr1	100	250
chr1	500	900	-
chr2	10	60
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** -s 序列文件、-r 区间文件、-o 输出文件，三者必填。-o 是输出文件名（合并模式），或作为分开输出时的文件前缀。

### 序列类型与处理 | Type & processing

**通俗理解|In plain words:** --type 选 dna（默认）还是 protein。--reverse-complement 把提取出的 DNA 全部反向互补（覆盖第 4 列的链信息）；--translate 把 DNA 翻译成蛋白（终止密码子记为 *）。蛋白序列不支持反向互补和翻译（会报错）。一般只需按数据实际类型给 --type，其余默认。

### 输出控制 | Output control

**通俗理解|In plain words:** 默认合并输出到一个文件；加 --separate-output 则每个区间一个文件（文件名带编号和坐标）。--no-headers 让序列名不包含区间信息（用 region_1、region_2 之类）；--line-width 控制 FASTA 每行字符数（默认 80）。一般不用动。

### 运行参数 | Runtime

**通俗理解|In plain words:** -t 线程数；--samtools-path 指定 samtools 路径（DNA 提取用，默认自动查找 samtools）。DNA 提取会自动按需建立 .fai 索引。

## 分析流程 | Pipeline { #pipeline }

```text
解析区间文件(校验坐标、识别 +/- 链)
    |
    v
DNA：若无 .fai 先 samtools faidx 建索引
    |
    v
按区间提取：
  DNA  -> samtools faidx（每 100 个区间一批）
  蛋白 -> Biopython / 基础方法切片
    |
    v
后处理：负链反向互补、--reverse-complement、--translate
    |
    v
合并或分开写入输出 FASTA（按 --line-width 换行）
```

## 输出 | Output { #output }

合并模式（默认）：

```text
extracted.fasta          # 提取出的序列
seq_extraction.log       # 日志（与输出文件同目录）
```

分开模式（--separate-output）：

```text
extracted_0001_chr1_100_250.fasta
extracted_0002_chr1_500_900_rc.fasta   # _rc 表示负链反向互补
...
```

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 打开输出 FASTA，每条序列名带「染色体:起始-终止」坐标，方便追溯来源；序列本体就是你要抠的那一段。

- 序列名格式（含区间信息时）：原名字_染色体:起始-终止；若原名字就是染色体名则直接写 染色体:起始-终止；负链额外加 rc_ 前缀
- 不存在的染色体或越界区间：会在日志里警告并跳过，不影响其余区间
- 翻译结果中的 * 是终止密码子，属正常

好坏判据：输出序列条数应约等于区间文件有效行数（除非有跳过）；每条序列长度 = 终止-起始+1（翻译/负链处理前）。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 基因组 DNA 区间提取：默认 --type dna，其余全默认
- 需要负链正确序列：第 4 列写 +/-，或加 --reverse-complement 全部反向互补
- 要蛋白序列：--type protein（注意区间坐标是蛋白残基位置）
- 要把 DNA 直接翻译成蛋白：--type dna 加 --translate
- 下游希望每段单独一个文件：加 --separate-output

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--sequence-file, -s` | 必填 |  | 输入序列文件(FASTA格式)｜Input sequence file (FASTA format) |
| `--regions-file, -r` | 必填 |  | 区域文件(BED格式:染色体 起始 终止)｜Regions file (BED format: chromosome start end) |
| `--output-file, -o` | 必填 | Path | 输出序列文件｜Output sequence file |
| `--type, --sequence-type` | `dna` | dna/protein | 序列类型｜Sequence type (default: dna) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--merge-output` | `True` |  | 合并输出到一个文件｜Merge output to one file |
| `--separate-output` | — |  | 分别输出到多个文件｜Output to separate files |
| `--no-headers` | — |  | 不包含区域信息在序列名中｜Do not include region info in sequence names |
| `--reverse-complement` | — |  | 反向互补DNA序列｜Reverse complement DNA sequences |
| `--translate` | — |  | 将DNA翻译为蛋白质｜Translate DNA to protein |
| `--line-width` | `80` | int | FASTA序列每行字符数｜Characters per line in FASTA sequence (default: 80) |
| `--verbose, -v` | `True` |  | 详细输出｜Verbose output |
| `--quiet, -q` | — |  | 安静模式｜Quiet mode |
| `--samtools-path` | `samtools` |  | samtools程序路径｜samtools program path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --sequence-file` | 必填 |  | 输入序列文件(FASTA格式)｜Input sequence file (FASTA format) |
| `-r, --regions-file` | 必填 |  | 区域文件(BED格式:染色体 起始 终止 [链])｜Regions file (BED format: chromosome start end [strand]) |
| `-o, --output-file` | 必填 |  | 输出序列文件｜Output sequence file |
| `--type, --sequence-type` | `dna` | dna/protein | 序列类型｜Sequence type |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--merge-output` | `True` | store_true | 合并输出到一个文件｜Merge output to one file |
| `--separate-output` | — | store_true | 分别输出到多个文件｜Output to separate files |
| `--no-headers` | — | store_true | 不包含区域信息在序列名中｜Do not include region info in sequence names |
| `--reverse-complement` | — | store_true | 反向互补DNA序列(会覆盖第4列的链信息)｜Reverse complement DNA sequences (overrides 4th column strand info) |
| `--translate` | — | store_true | 将DNA翻译为蛋白质｜Translate DNA to protein |
| `--line-width` | `80` | int | FASTA序列每行字符数｜Characters per line in FASTA sequence |
| `-v, --verbose` | `True` | store_true | 详细输出｜Verbose output |
| `-q, --quiet` | — | store_true | 安静模式｜Quiet mode |
| `--samtools-path` | `samtools` |  | samtools程序路径｜samtools program path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- samtools（DNA 提取与建索引，自动检测 conda 环境并用 conda run 调用）
- Biopython（蛋白提取/翻译优先用它，缺失时回退到基础方法）
- Python 3

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持整体续传；但 .fai 索引若已存在会跳过重建（复用）。每次运行会重新提取并覆盖输出文件。

Q2：坐标是 0-based 还是 1-based？
本工具按 1-based 闭区间处理（起始、终止都包含），与 samtools faidx 一致，但与常见 BED 的 0-based 半开不同，写区间时注意不要照搬 0-based BED 坐标。

Q3：线程数默认是多少？
通过 biopytools parse-seq 进入时默认 12；直接 python -m biopytools.parse_seq.main 进入时默认 88，两者不一致，注意按实际入口确认。

Q4：负链怎么处理？
区间第 4 列为 - 时自动反向互补，序列名前加 rc_；也可用 --reverse-complement 强制全部反向互补（此时覆盖第 4 列的链信息）。

Q5：蛋白序列能反向互补/翻译吗？
不能。--type protein 同时给 --reverse-complement 或 --translate 会直接报错。

Q6：区间越界会怎样？
对应区间在日志中警告并跳过，不影响其它区间；最终输出条数可能少于区间行数。