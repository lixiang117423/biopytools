# seq2genome - 序列到基因组比对 | Sequence to Genome Alignment

一句话理解：**把一批查询序列（DNA 或蛋白）定位到基因组上，告诉你每条序列比对到哪条染色体的哪个位置、有多像，并顺手导出 PAF/GFF3/BED/统计报告和比对区间的基因组序列。**

## 功能概述 | Overview { #overview }

- 自动检测查询序列类型：DNA 用 minimap2、蛋白用 miniprot（也可手动指定）
- 输出 PAF 原始比对、统计报告、GFF3、BED、提取的比对区间 FASTA
- 统计比对数、一致性、覆盖度、Mapping Quality 等
- 各输出可用 --no-gff3 / --no-bed / --no-statistics 单独关闭

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools seq2genome --genome genome.fa --query sequences.fa -o results
```

最小输入：一个基因组 FASTA + 一个查询序列 FASTA（DNA 或蛋白，默认自动识别）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 参考基因组(genome) | 拼接好的物种全基因组，比对的目标「地图」 |
| 查询序列(query) | 要定位到基因组上的序列，像要找位置的「人名」 |
| 比对(alignment) | 把查询序列放到基因组上找最匹配的位置 |
| DNA 序列 | 只含 ATCGN 的核酸序列 |
| 蛋白序列 | 20 种氨基酸组成的序列，比对时更看重「相似」而非「相同」 |
| minimap2 | DNA 比对到 DNA 的常用快速工具 |
| miniprot | 蛋白比对到基因组的专用工具 |
| PAF | 比对的通用文本格式，一行一条比对记录 |
| GFF3 / BED | 两种基因组注释格式，方便在基因组浏览器里看 |
| identity / 覆盖度 | 一致性（有多像）/ 覆盖度（query 有多少比例被比对到） |

## 输入 | Input { #input }

- 基因组 FASTA(--genome / -g)：参考基因组。
- 查询 FASTA(--query / -q)：DNA 或蛋白序列；默认自动检测类型（DNA 判据 ATCGN 占比大于 90%，蛋白判据氨基酸字符占比大于 30%）。

示例：

```text
>gene1
ATGGCTAGCTAGCTAGCTA
>gene2
ATGGCTAGCGAGCTAGCTA
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** -g 基因组、-q 查询序列、-o 输出目录（默认 ./seq2genome_output），三者必填。

### 序列类型 | Sequence type

**通俗理解|In plain words:** --query-type 选 dna / protein / auto（默认 auto 自动检测）。自动检测偶尔会判断错，判断不出来会报错让你手动指定。一般保持 auto 即可；确定类型时手动指定更稳。

### 运行参数 | Runtime

**通俗理解|In plain words:** -t 线程数（默认 12）；--minimap2-path / --miniprot-path 指定工具路径（默认直接调用 minimap2 / miniprot）。只有工具不在 PATH 里时才需要改路径。

### 输出开关 | Output toggles

**通俗理解|In plain words:** 默认导出全部（PAF、统计、GFF3、BED、提取序列）。不需要哪个就加 --no-gff3 / --no-bed / --no-statistics 关掉。注意：提取比对区间序列（alignment.fa）依赖 BED 输出，关闭 BED 就不会提取。

## 分析流程 | Pipeline { #pipeline }

```text
自动检测查询序列类型（DNA / 蛋白）
    |
    v
比对：DNA -> minimap2 -x asm5；蛋白 -> miniprot
    |
    v
解析 PAF 比对结果
    |
    v
生成统计报告 alignment_statistics.txt
    |
    v
导出 GFF3、BED
    |
    v
用 seqkit 按 BED 提取基因组序列 alignment.fa
```

## 输出 | Output { #output }

```text
results/
|-- alignment.paf                    # 原始比对结果（PAF）
|-- alignment_statistics.txt         # 统计报告
|-- alignment.gff3                   # GFF3 注释
|-- alignment.bed                    # BED 注释（BED6）
|-- alignment.fa                     # 比对区间的基因组序列
-- pep2genome_processing_<时间戳>.log  # 日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

**通俗理解|In plain words:** 先看统计报告把握全局（比上多少条、平均一致性多高），再看 PAF/GFF3/BED 定位到具体位置，最后用 alignment.fa 拿序列。

- alignment.paf：每条查询序列的比对记录，PAF 前 12 列含 query/target 名称、起止、链方向、匹配数、比对长度、Mapping Quality 等
- alignment_statistics.txt：分 5 节——总体统计、查询序列统计、目标序列统计、比对质量统计、覆盖度统计；核心看平均 identity、Mapping Quality、覆盖度
- alignment.gff3：把比对画成基因组上的区间（feature 为 protein_match），属性里带 Identity、覆盖度
- alignment.bed：BED6 格式（染色体、起止、名字、得分、链），可喂给浏览器或下游工具
- alignment.fa：从基因组抠出的比对区间序列

好坏判据：比对成功则 alignment.paf 非空；identity 平均较高（DNA 通常 95% 以上，蛋白看物种差异）、Mapping Quality 较高说明定位可信。多条序列比对到同一位置提示可能存在多拷贝/重复。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规用法：全部默认（auto 检测类型、全部输出）
- 只要比对和统计、不要注释：加 --no-gff3 --no-bed
- 确定是蛋白/ DNA 时：手动 --query-type protein / dna 更稳
- 工具不在 PATH：用 --minimap2-path / --miniprot-path 指定完整路径

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--query, -q` | 必填 |  | 查询序列FASTA文件（DNA或蛋白质）｜Query sequence FASTA file (DNA or protein) |
| `--query-type` | `auto` | dna/protein/auto | 序列类型（默认：auto自动检测）｜Sequence type (default: auto for auto-detection) |
| `--output-dir, -o` | `./seq2genome_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--minimap2-path` | `minimap2` |  | Minimap2工具路径｜Minimap2 tool path |
| `--miniprot-path` | `miniprot` |  | Miniprot工具路径｜Miniprot tool path |
| `--no-gff3` | — |  | 不导出GFF3格式｜Do not export GFF3 format |
| `--no-bed` | — |  | 不导出BED格式｜Do not export BED format |
| `--no-statistics` | — |  | 不生成统计报告｜Do not generate statistics report |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--query` | 必填 |  | 查询序列FASTA文件（DNA或蛋白质）｜Query sequence FASTA file (DNA or protein) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--query-type` | `auto` | dna/protein/auto | 查询序列类型（默认：auto自动检测）｜Query sequence type (default: auto for auto-detection) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minimap2-path` | `minimap2` |  | Minimap2工具路径｜Minimap2 tool path (default: minimap2) |
| `--miniprot-path` | `miniprot` |  | Miniprot工具路径｜Miniprot tool path (default: miniprot) |
| `--no-gff3` | — | store_false | 不导出GFF3格式｜Do not export GFF3 format |
| `--no-bed` | — | store_false | 不导出BED格式｜Do not export BED format |
| `--no-statistics` | — | store_false | 不生成统计报告｜Do not generate statistics report |
| `--no-extract` | — | store_false | 不提取基因组序列｜Do not extract genome sequences |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- minimap2（DNA 查询，默认命令 minimap2）
- miniprot（蛋白查询，默认命令 miniprot）
- seqkit（提取比对区间序列，命令 seqkit subseq）
- Biopython（序列类型自动检测）、numpy（统计报告）
- Python 3

## 常见问题 | FAQ { #faq }

Q1：支持断点续传吗？
不支持。每次运行都全量重新比对并覆盖输出。

Q2：DNA 和蛋白怎么区分？
默认 auto 自动检测：ATCGN 占比大于 90% 判 DNA，氨基酸字符占比大于 30% 判蛋白；都判不出会报错，此时用 --query-type 手动指定。

Q3：为什么没有 alignment.fa？
提取该文件用 seqkit subseq 完成，若没装 seqkit，这一步会失败但不影响其它输出（PAF/GFF3/BED 照常生成）。装上 seqkit 即可。

Q4：GFF3 里 source 和 type 为什么都写 miniprot / protein_match？
这是当前实现的已知小瑕疵：GFF3 导出器把 source 固定写成 miniprot、type 固定写成 protein_match，即使 DNA 用 minimap2 比对也是如此。不影响区间坐标的准确性。

Q5：比对结果为空怎么办？
若 PAF 解析出 0 条记录，程序会警告并退出。检查查询序列与基因组是否同物种、类型是否选对（DNA/蛋白弄反最常见）。

Q6：想关掉序列提取可以吗？
biopytools seq2genome 包装器未暴露 --no-extract（该参数只在直接 python -m biopytools.seq2genome.main 时可用），因此经 biopytools 入口时提取步骤默认总是执行，需装 seqkit。