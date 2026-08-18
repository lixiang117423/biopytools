# 染色体重命名 | Chromosome Rename (minimap2)

一句话理解：**把新组装基因组里「group1、group2…」这种看不懂的序列名，通过和参考基因组做全基因组比对，自动换成标准的 Chr1、Chr2… 命名，并按参考染色体的顺序排好输出**。

## 功能概述 | Overview

- 基于 minimap2 全基因组比对（支持 asm5/asm10/asm20 三种预设）
- 把待重命名基因组的每条序列（通常是 group/contig 命名）对应到参考基因组的标准染色体名
- 长序列优先处理、每条参考染色体只被映射一次，避免重复占用
- 输出映射表、汇总报告、以及重命名并排好序的 FASTA
- 未映射上的序列自动命名为 Contig_01、Contig_02…（按长度降序）

## 快速开始 | Quick Start

```bash
biopytools chr-rename -r reference.fa -q genome.fa -o output
```

最小输入：一个参考基因组 FASTA（标准染色体命名）+ 一个待重命名的基因组 FASTA（group 命名）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 参考基因组(reference) | 已经命名规范、当作「标准答案」的基因组，工具拿它当参照物 |
| 待重命名基因组(query) | 你新组装出来、序列名还是 group1/contig5 这种「临时编号」的基因组 |
| group 命名 | 组装软件默认起的编号名（如 group1、unitig_2），人看不懂哪条对应哪条染色体 |
| Chr 命名 | 约定俗成的标准名（如 Chr1、Chr2、ChrX），下游分析工具普遍认识 |
| 比对(alignment) | 把两条序列「对在一起」，看它们有多少地方长得一样 |
| PAF 文件 | minimap2 输出的比对结果表，一行记录一段比对，本工具据此判断对应关系 |
| 覆盖度(coverage) | 一条 query 序列「被参考序列覆盖的比例」，覆盖得越多说明越像同一条 |
| 一致性(identity) | 比对上的部分里「碱基完全一致」的比例，越高越可靠 |
| 嵌合/错误组装 | 一条序列其实是两条染色体「拼错了」接在一起，会表现为同时比对上多个参考 |

## 输入 | Input

两个 FASTA 文件（.fa / .fasta）：

- **参考基因组**：序列名应是标准染色体命名（Chr1、Chr2…），作为命名的参照；
- **待重命名基因组**：序列名是 group/contig 等临时编号，将被换成参考染色体的名字。

示例（参考基因组）：

```text
>Chr1
ACGTACGT...
>Chr2
TTGCAACG...
```

示例（待重命名基因组）：

```text
>group1
ACGTACGT...
>group2
TTGCAACG...
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这是「参照物」和「待改名对象」两个输入，缺一不可。工具会把 -q 里的每条序列拿去和 -r 比对，再照着 -r 的名字给 -q 改名。

### 输出与软件路径 | Output & software path

**通俗理解|In plain words:** -o 决定结果放哪，一般不用动；-a 是 minimap2 的路径，默认自动解析 align 功能域环境，缺失时回退 PATH 直接调用，装了但没进 PATH 才需要指定完整路径。

### 比对参数 | Alignment parameters

**通俗理解|In plain words:** -x（preset）决定比对时按「两条序列差异多大」来设参数：asm5 适合差异小于 5% 的近缘基因组，asm10/asm20 适合差异更大的。 -t 是线程数，越大越快、越吃 CPU，一般不用动。

### 过滤参数 | Filtering parameters

**通俗理解|In plain words:** -l（最小比对长度）决定「多长的比对才算数」——比对段太短（如几 kb）多半是局部相似或重复序列，不可靠，会被丢掉。-i（最小一致性）这个参数在配置里存在、也会写进汇总报告，但**当前映射判定实际不依赖它**（映射用的是固定的 20% 覆盖度阈值，见 FAQ），所以基本不用管它。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先全基因组比一遍，再按「谁长谁先挑、一个参考只配一条」的规则建立对应关系，最后照着对应关系改名并排序输出。

```text
参考 FASTA + 待重命名 FASTA
    │
    ▼
步骤1: minimap2 全基因组比对 → alignment.paf
    │
    ▼
步骤2: 解析 PAF，计算每条 query 对每个参考序列的合并覆盖度
    │
    ▼
步骤3: 建立映射（长序列优先）
    ├─ 按 query 长度降序处理
    ├─ 合并覆盖度 >= 20% 且总比对长度 >= 最小比对长度
    ├─ 每个参考序列只被映射一次
    └─ 选覆盖度最高的参考作为目标
    │
    ▼
步骤4: 检测一对多映射（多个 group 映射到同一参考 → 疑似碎片化组装）
    │
    ▼
步骤5: 写映射表 + 汇总报告
    │
    ▼
步骤6: 重命名输出（按参考染色体顺序，未映射的命名为 Contig_01…）
```

## 输出 | Output

```text
output/
├── alignment.paf            # minimap2 原始比对结果
├── chromosome_mapping.tsv   # 映射表：原序列名 → 新染色体名 + 长度/覆盖度/一致性
├── renamed_genome.fa        # 重命名并排好序的基因组（最终产物）
├── rename_summary.txt       # 汇总报告（映射统计、一对多映射清单）
└── chr_rename.log           # 运行日志
```

## 结果解读 | Interpreting Results

### 1. 映射表（chromosome_mapping.tsv）

**通俗理解|In plain words:** 这是一张「旧名 → 新名」对照表，每行告诉你原来的 group 现在叫什么、对应关系有多可靠。

- Original_Chromosome：原来的序列名（group 命名）
- Renamed_Chromosome：新名字（参考染色体的标准名）
- Coverage：覆盖度，越接近 100% 说明越确定是同一条
- Identity：一致性，越高越可靠

### 2. 重命名后的 FASTA（renamed_genome.fa）

**通俗理解|In plain words:** 这是最终要用的文件——序列已按参考染色体顺序排好，未映射的序列排在最后并命名成 Contig_01、Contig_02…。

- 已映射序列：按参考基因组染色体顺序输出，名字换成参考染色体名；
- 未映射序列：按长度降序命名 Contig_01、Contig_02…。

### 3. 汇总报告（rename_summary.txt）

**通俗理解|In plain words:** 一份「体检小结」，重点看「成功映射了多少条」和「一对多映射」清单。

- 成功映射数：占 query 总序列的比例越高越好；
- **一对多映射**：多个 group 映射到同一条参考染色体，通常意味着**基因组组装碎片化**（一条真染色体被拼成了好几段），值得检查组装质量。

## 参数选择建议 | Parameter Guidance

- **-x, --preset**：同物种/近缘物种用默认 asm5；跨物种或差异较大时用 asm10/asm20
- **-l, --min-alignment-length**：默认 100000（100kb）通常合适；重复序列多的基因组可适当调大，碎片化严重的小 contig 多时可调小
- **-t, --threads**：机器核数多就调大，默认 12 通常够用

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref, -r` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--query, -q` | 必填 |  | 待重命名的基因组FASTA文件｜Query genome FASTA file to rename |
| `--output-dir, -o` | `./chr_rename_output` | Path | 输出目录｜Output directory |
| `--minimap2-path, -a` | `minimap2` |  | minimap2软件路径｜minimap2 software path |
| `--preset, -x` | `asm5` | asm5/asm10/asm20 | minimap2预设模式｜minimap2 preset mode |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-identity, -i` | `0.9` | float | 最小序列一致性阈值(0-1)｜Minimum identity threshold (0-1) |
| `--min-alignment-length, -l` | `100000` | int | 最小比对长度(bp)｜Minimum alignment length (bp) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --ref` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-q, --query` | 必填 |  | 待重命名的基因组FASTA文件｜Query genome FASTA file to rename |
| `-o, --output-dir` | `./chr_rename_output` |  | 输出目录｜Output directory |
| `-a, --minimap2-path` | — |  | minimap2软件路径(默认域环境自动解析)｜minimap2 software path (default: auto domain env) |
| `-x, --preset` | `asm5` | asm5/asm10/asm20 | minimap2预设模式｜minimap2 preset mode (asm5/asm10/asm20) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-i, --min-identity` | `0.9` | float | 最小序列一致性阈值(0-1)｜Minimum identity threshold (0-1) |
| `-l, --min-alignment-length` | `100000` | int | 最小比对长度(bp)｜Minimum alignment length (bp) |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- minimap2（自动解析 align 域环境并经 conda run 调用；可用 -a/--minimap2-path 或环境变量 MINIMAP2_PATH 覆盖；域环境缺失时回退 PATH 直接调用）

## 常见问题 | FAQ

**Q1：-i, --min-identity 设置了为什么好像没起作用？**
映射判定当前使用的是**固定 20% 覆盖度阈值**（合并覆盖度 >= 20% 且总比对长度 >= -l），并不依赖 -i 的一致性阈值。-i 目前只写进汇总报告，实际不参与筛选。

**Q2：为什么一条参考染色体被多个 group 对应上了？**
这是「一对多映射」，通常说明你的组装是**碎片化的**——一条真染色体被拼成了多个 group。工具会保留最长的那个映射，并在汇总报告里列出所有相关 group 供你检查。

**Q3：没映射上的序列去哪了？**
不会丢。它们会按长度降序命名成 Contig_01、Contig_02…，排在重命名 FASTA 的最后。

**Q4：支持断点续传吗？**
不支持。每次运行都会重新跑 minimap2 比对和全部后续步骤；想换参数重跑直接重新执行即可，旧输出会被覆盖。

**Q5：改完名之后序列内容变了吗？**
没变。工具只改序列名和排列顺序，序列本身原样保留。

