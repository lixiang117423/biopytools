# RagTag 参考基因组挂载 | RagTag Reference-based Scaffolding

一句话理解：**拿一条高质量的参考基因组当「模板」，把你的散乱 contig 按参考的顺序、方向和间距排好并连成染色体级序列**。
输入一条参考 FASTA、一条待挂载的 contig FASTA 和一个样品名，输出挂上参考的序列、没挂上的序列、以及两者的合并版。

## 功能概述 | Overview

- 参考基因组引导的 scaffolding：用参考基因组的染色体结构给 contig「排座位」
- 三种比对器可选：`minimap2`（默认）/ `unimap` / `nucmer`
- 自动重命名序列 ID 并分类：挂上参考的（chr/chrom 命名）与没挂上的（scaffold_N 命名）分开输出
- 支持把未定位 contig 合并成 `chr0`（`-C`）和推断 gap 大小（`-R`）
- 可选给所有输出序列 ID 加统一前缀（`-p/--prefix`）

## 快速开始 | Quick Start

```bash
biopytools ragtag -r ref.fa -q query.fa -s Sample1 -o output_dir/
```

最小输入：一条参考基因组 + 一条待挂载的 contig 组装 + 一个样品名（用于输出文件命名）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| contig | 一段连起来、没有缺口的序列，是最小的组装单元 |
| scaffold | 多段 contig 用「缺口(gap)」连起来的更大单元 |
| scaffolding | 把 contig 按正确顺序和方向连成 scaffold 的过程 |
| 参考基因组 | 一条已知、质量高的同类基因组，当「标准答案」用 |
| 未定位 contig | 没能挂到任何参考染色体上的 contig |
| gap | 序列里一段未知的空白（用 N 表示），长度可能由工具推断 |
| AGP | 描述「哪段序列由哪些片段、按什么顺序拼起来」的文本文件 |
| MAPQ | 比对质量分，越高说明这段比对越可信 |

## 输入 | Input

### 参考基因组 FASTA

用于「排座位」的高质量参考基因组。最好和待挂载物种亲缘关系近，同物种最佳。

### 查询(待挂载) FASTA

你的 contig 级组装，序列 ID 可以是任意形式，工具会自动识别并重命名。

### 样品名

必填，用于给输出文件命名（如 `Sample1_RagTag_scaffolded.fa`）。

## 参数说明 | Parameters

### 输入与命名 | Input & naming

**通俗理解|In plain words:** `-p/--prefix` 会在所有输出序列 ID 前面加个统一前缀（如 `SpeciesA_`），用来区分不同来源的组装，不需要区分就留空。

### 比对器 | Aligner

**通俗理解|In plain words:** `--aligner` 决定「contig 和参考怎么对上」。默认 `minimap2` 又快又好，**一般不用动**；`unimap` 更注重长读长的唯一比对；`nucmer` 是经典老牌比对器，对高度相似的近缘参考更稳，但慢一些。

### 挂载行为 | Scaffolding behavior

**通俗理解|In plain words:** `-C/--concatenate-unplaced` 把没挂上的 contig 都串成一条假染色体 `chr0`，方便有些下游工具「要求每条序列都算染色体」；`-R/--infer-gaps` 让工具按比对证据推测 gap 长度（否则 gap 用固定长度），想得到更真实的总长就开它。两者默认都关，**一般不用动**。

## 分析流程 | Pipeline

```text
参考 FASTA + 查询 contig FASTA + 样品名
    │
    ▼
步骤1: 运行 ragtag.py scaffold(minimap2 比对 → 排序 → 定向 → 连 gap)
        输出 ragtag.scaffold.fasta / .agp / .confidence.txt 等
    │
    ▼
步骤2: 解析输出,按原ID判断是否挂上参考,重命名并分类
        chr/chrom 开头 → scaffolded;ptg/tig 开头 → scaffold_N
    │
    ▼
步骤3: 写三个 FASTA(scaffolded / unscaffolded / combined)并重命名原始输出
```

## 输出 | Output

```text
output_dir/
├── Sample1_RagTag_scaffolded.fa      # 挂上参考的序列(chr/chrom 命名)
├── Sample1_RagTag_unscaffolded.fa    # 没挂上的序列(scaffold_1, scaffold_2...)
├── Sample1_RagTag_combined.fa        # 上面两者合并(完整组装)
├── Sample1_ragtag.scaffold.fasta     # RagTag 原始 scaffold 输出
├── Sample1_ragtag.scaffold.agp       # 挂载关系(AGP)
├── Sample1_ragtag.scaffold.confidence.txt  # 各 contig 挂载置信度
├── Sample1_ragtag.scaffold.stats     # 挂载统计
├── Sample1_ragtag.scaffold.asm.paf   # 比对记录
└── Sample1_ragtag.log                # 运行日志
```

- `Sample1_RagTag_combined.fa`：最终完整组装，下游分析一般用这个
- `Sample1_RagTag_scaffolded.fa`：只含挂上参考的序列，染色体级
- `Sample1_RagTag_unscaffolded.fa`：没挂上的 contig，别丢，它们往往含新序列
- `Sample1_ragtag.scaffold.confidence.txt`：每条 contig 挂到哪、置信度多高，复核挂载质量的关键文件

## 结果解读 | Interpreting Results

- **scaffolded 越多越好**：挂上参考的序列占比高，说明组装和参考结构一致
- **confidence.txt 里置信度低的位置要留意**：低置信度可能意味着 contig 放错位置或方向，建议人工复核或后续用 Hi-C 再校正
- **unscaffolded 里出现大段序列**：可能是参考里没有的插入/新片段，未必是坏事
- **ID 命名规则**：原 ID 以 `chr`/`chrom` 开头或带 `_RagTag` 后缀 → 视为已挂载；以 `ptg`/`tig` 开头 → 视为未挂载并改名为 `scaffold_N`

## 参数选择建议 | Parameter Guidance

- **默认直接用**：`minimap2` + 不额外开任何开关，最省心
- **参考很近缘、要最高精度**：可试 `--aligner nucmer`
- **下游工具要求每条都是染色体**：加 `-C`
- **想要更真实的总长和 gap**：加 `-R`
- **要区分多个来源**：用 `-p SpeciesA` 加前缀

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --reference` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-q, --query` | 必填 |  | 查询基因组FASTA文件｜Query genome FASTA file |
| `-s, --sample-name` | 必填 |  | 样品名称，用于输出文件命名｜Sample name for output file naming |
| `-p, --prefix` | — |  | 序列ID前缀，会添加到所有输出序列的ID前面｜Sequence ID prefix to add to all output sequence IDs |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-o, --output-dir` | `./ragtag_output` | Path | 输出目录｜Output directory |
| `--aligner` | `minimap2` | minimap2/unimap/nucmer | 比对器｜Aligner to use |
| `-C, --concatenate-unplaced` | — |  | 将未定位的contigs合并为chr0｜Concatenate unplaced contigs into chr0 |
| `-R, --infer-gaps` | — |  | 推断gap大小｜Infer gap sizes |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --reference` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-q, --query` | 必填 |  | 查询基因组FASTA文件｜Query genome FASTA file |
| `-s, --sample-name` | 必填 |  | 样品名称，用于输出文件命名｜Sample name for output file naming |
| `-p, --prefix` | — |  | 序列ID前缀，会添加到所有输出序列的ID前面｜Sequence ID prefix to add to all output sequence IDs |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-o, --output-dir` | `./ragtag_output` |  | 输出目录｜Output directory |
| `--aligner` | `minimap2` | minimap2/unimap/nucmer | 比对器｜Aligner to use |
| `-C, --concatenate-unplaced` | — | store_true | 将未定位的contigs合并为chr0｜Concatenate unplaced contigs into chr0 |
| `-R, --infer-gaps` | — | store_true | 推断gap大小｜Infer gap sizes |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- RagTag（通过命令行 `ragtag.py` 调用，需在 PATH 中，未绑定 conda 环境）
- 比对器：`minimap2`（默认）/ `unimap` / `nucmer`，需在 PATH 中

## 常见问题 | FAQ

**Q1：重跑会不会跳过已完成的步骤？**
不会。本模块**没有断点续传**，每次运行都会重新执行 RagTag scaffold 并覆盖旧产物。换参数请换输出目录。

**Q2：为什么输出里序列名变了？**
工具会按原 ID 判断是否挂载并重命名：`chr/chrom` 开头或带 `_RagTag` 后缀视为已挂载（去掉 `_RagTag` 后缀），`ptg/tig` 开头视为未挂载并改成 `scaffold_N`。若希望保持某类命名，先调整输入 FASTA 的序列 ID。

**Q3：unscaffolded 里的序列能要吗？**
能。没挂上参考不代表没用，可能是参考里没有的新片段。下游分析建议用合并后的 `combined.fa`，而不是只留 scaffolded。

**Q4：RagTag 报错说找不到 ragtag.py 怎么办？**
本模块直接调用 PATH 里的 `ragtag.py`。先确认 `which ragtag.py` 能返回路径；没有就先安装 RagTag 或把它的脚本目录加进 PATH。

**Q5：`-C` 生成的 chr0 是什么？**
是所有未定位 contig 被串成的一条「假染色体」，只为满足某些下游工具「每条序列都叫染色体」的要求，不代表真实的染色体。
