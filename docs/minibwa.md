# minibwa - 短读长比对分析（标准/Hi-C/BS-seq/长读） | Minibwa Short-read Alignment

一句话理解：**把一批 FASTQ 测序读段「对」到参考基因组上，得到每个样本的 BAM 比对文件，并自动统计比对率、计算覆盖度**。支持普通短读、Hi-C、甲基化（BS-seq）、长读四种模式。

## 功能概述 | Overview

- 用 minibwa 做比对（自动构建索引、批量比对），samtools 做排序/索引/统计/覆盖度
- 四种模式：standard（普通双端/单端短读）、hic（Hi-C）、meth（方向性 BS-seq）、long（准确长读）
- 自动构建参考基因组索引（BS-seq 模式额外建 .meth.mbw），多样本共享一次索引
- 每样本产出排序 BAM、比对统计（flagstat/stats）、覆盖度（逐位深度 + 滑窗）
- 可选标记/移除重复，支持断点续传（--resume 跳过已完成样本）

## 快速开始 | Quick Start

```bash
biopytools minibwa -g genome.fa -i fastq_dir -o results -t 16
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 参考基因组 | 该物种的「标准答案」序列，比对就是把读段贴到它上面找位置 |
| 比对（alignment） | 给每条读段在基因组上找到最像的落点，像查字典定位单词 |
| BAM | 比对结果的标准二进制文件，记录每条读段落在基因组哪里、怎么落的 |
| 索引 | 给基因组预先建的「目录」，让比对更快，相当于书的索引页 |
| Hi-C | 一种捕获染色体三维空间相邻关系的测序，读段间有特殊配对模式 |
| BS-seq | 甲基化测序，C→T 的转化使比对要特殊处理（--meth 模式） |
| 覆盖度 | 基因组每个位置被测了几遍，越高越可靠 |
| 重复（duplicate） | PCR 复制出来的相同读段，不算新信息，可标记或移除 |

## 输入 | Input

- 参考基因组 FASTA（-g）：比对的目标序列。
- FASTQ 输入目录（-i）：存放测序读段，自动按 R1 模式配对。

默认 R1 匹配模式为 _1.fq.gz（-p 可改），程序会把 _1 换成 _2 找配对文件；找不到 R2 的样本按单端处理。

## 参数说明 | Parameters

### 必需与模式参数 | Required and mode

**通俗理解|In plain words:** -g 基因组、-i FASTQ 目录、-p R1 命名模式是三个基本输入；--mode 选比对模式（standard/hic/meth/long），按你的实验类型选，默认 standard。普通转录组/WGS 用 standard，Hi-C 用 hic，甲基化用 meth，长读用 long。

### 性能参数 | Performance

**通俗理解|In plain words:** -t 线程数默认 12，机器核多可调大加速。其余 minibwa 的预设/种子/带宽/打分等参数（--preset、-k、-c、-w、-A、-B 等）都是比对的精细旋钮，默认值经实践验证，一般不用动；除非你明确知道要调灵敏度/速度。

### 打分与 IO 参数 | Scoring and I/O

**通俗理解|In plain words:** -A/-B/-O/-E 是匹配得分、错配罚分、gap 罚分，决定比对「严不严」，默认值通用。--read-group 给 BAM 加样本组标签（合规要求时用）；-u 不输出未比对读段（省空间）；--outn 控制输出次要比对数量。这些都用默认即可。

### 后处理参数 | Post-processing

**通俗理解|In plain words:** --markdup 标记重复读段（不改数据、只加标签），--remove-dup 直接删掉重复（需同时 --markdup）。做变异检测建议开 --markdup；做表达量/覆盖度通常不删。

### 覆盖度参数 | Coverage

**通俗理解|In plain words:** 默认会算每个样本的逐位深度和滑窗覆盖度（窗口 1M、步长 100k）。--skip-coverage 跳过（不需要覆盖度时省时间）；--min-base-quality/--min-mapping-quality 只统计质量达标的碱基/比对；--max-depth 限制统计的最大深度（0=不限）。一般不用动。

### 工具路径与运行控制 | Tool paths and run control

**通俗理解|In plain words:** --minibwa-path/--samtools-path 是软件位置，默认路径一般已配好。--resume 开启断点续传——重跑时跳过已有 BAM+索引的样本，只补失败的。

## 分析流程 | Pipeline

```text
参考基因组 FASTA + FASTQ 目录
    │
    ▼
步骤1：检查依赖（minibwa、samtools）
    │
    ▼
步骤2：构建/复用基因组索引（01_index/，多样本共享）
    │
    ▼
步骤3：查找 FASTQ 配对样本
    │
    ▼
步骤4：逐样本 minibwa map | samtools sort → BAM → （可选 markdup）→ 索引
        → flagstat/stats → 覆盖度（depth + 滑窗）
    │
    ▼
步骤5：汇总报告 + 软件版本/参数元数据
```

## 输出 | Output

```text
results/
├── 00_pipeline_info/
│   ├── software_versions.yml        # 软件版本
│   ├── pipeline_params.yaml         # 运行参数
│   └── alignment_summary.txt        # 多样本汇总报告
├── 01_index/                        # 参考基因组索引（ref.l2b/.mbw，meth 加 .meth.mbw）
├── <样本名>/
│   ├── <样本>.minibwa.bam           # 排序比对结果（+ .bam.bai 索引）
│   ├── <样本>.minibwa.markdup.bam   # 标记重复后（仅 --markdup）
│   ├── <样本>.minibwa.flagstat.txt  # 比对统计（flagstat）
│   ├── <样本>.minibwa.stats.txt     # 详细统计（samtools stats）
│   ├── <样本>.minibwa.depth.txt     # 逐位覆盖深度
│   └── <样本>.minibwa.window.bed    # 滑窗平均覆盖度
└── 99_logs/minibwa_pipeline.log     # 运行日志
```

- <样本>.minibwa.bam：核心比对结果，后续变异检测/定量都从它出发。
- <样本>.minibwa.flagstat.txt：比对率、配对率等关键指标（多少读段比对上了）。
- <样本>.minibwa.window.bed：四列（染色体、窗口起、窗口止、平均深度），便于画覆盖度曲线。
- alignment_summary.txt：所有成功样本的清单与关键参数汇总。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看 flagstat 里的比对率（越高说明数据越干净、和参考基因组越匹配），再看覆盖度是否均匀。

- 比对率：flagstat 中 mapped 占比通常应 >90%；偏低提示数据污染、参考基因组不对、或接头未去除。
- 配对率：properly paired 占比，双端数据应较高，异常低提示插入片段或文库问题。
- 覆盖度：depth.txt 逐位深度与 window.bed 平均深度，用来判断测序是否均匀；窗口覆盖度波动大可能提示重复或 GC 偏差。
- 重复率：若开了 --markdup，markdup_stats 会给出重复比例，过高提示建库 PCR 过度。
- 失败样本：日志里标记比对/统计/覆盖度失败的样本，汇总报告只列成功样本。

## 参数选择建议 | Parameter Guidance

- 常规 WGS/RNA-seq：默认 --mode standard 即可。
- Hi-C 数据：必须 --mode hic；甲基化数据：--mode meth；准确长读：--mode long。
- 做变异检测：加 --markdup（必要时 --remove-dup）；只做覆盖度/定量：默认即可。
- 机器核多：调大 -t；内存/磁盘紧张：可 --skip-coverage 减少输出。
- 大批样本中断后重跑：加 --resume，只补未完成样本。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `-i, --input` | 必填 |  | FASTQ输入目录｜FASTQ input directory |
| `-p, --pattern` | `_1.fq.gz` |  | R1匹配模式｜R1 pattern |
| `-o, --output-dir` | `./minibwa_output` |  | 输出目录｜Output directory |
| `--mode` | `standard` | standard/hic/meth/long | 比对模式｜Alignment mode |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--preset` | `adap` |  | -x 预设｜preset (sr/lr/adap) |
| `-k, --min-seed` | `19` | int | -k 最小种子长度｜min seed length |
| `-c, --max-occ` | `250` | int | -c 最大种子出现次数｜max seed occurrences |
| `--max-gap` | `100` | int | minibwa -g 最大gap｜max gap size |
| `-w, --bandwidth` | `100` | int | -w 带宽｜bandwidth |
| `-W, --bandwidth-long` | `20000` | int | -W 长读带宽｜long bandwidth |
| `-m, --min-chain-score` | `25` | int | minibwa -m 最小链接分数｜min chaining score |
| `--sec-ratio` | `0.5` | float | -p 次要/主要得分比｜secondary-to-primary ratio |
| `-N, --max-sec` | `50` | int | -N 保留次要比对数｜retain N secondary alignments |
| `-s, --min-dp-score` | `30` | int | -s 最小DP得分｜min DP score |
| `-A, --match-score` | `2` | int | -A 匹配得分｜matching score |
| `-B, --mismatch-penalty` | `8` | int | -B 错配罚分｜mismatch penalty |
| `-O, --gap-open` | `12,23` |  | -O gap开放罚分｜gap open penalty |
| `-E, --gap-ext` | `2,1` |  | -E gap延伸罚分｜gap extension penalty |
| `-R, --read-group` | — |  | -R SAM read group｜read group line |
| `-u, --no-unmap` | — |  | -u 不输出未比对read｜do not output unmapped |
| `--outn` | `0` | int | --outn 输出次要比对上限｜max secondary to output |
| `-y, --copy-comment` | — |  | -y 复制FASTA/Q注释｜copy comments |
| `-Y, --soft-clip-supp` | — |  | -Y 软剪切补充比对｜soft clip supplementary |
| `-K, --batch-size` | `100m,1g` |  | -K 批处理大小｜batch size |
| `--markdup` | — |  | 标记重复｜Mark duplicates |
| `--remove-dup` | — |  | 移除重复（需--markdup）｜Remove duplicates (requires --markdup) |
| `--skip-coverage` | — |  | 跳过覆盖度分析｜Skip coverage analysis |
| `--min-base-quality` | `0` | int | 最小碱基质量｜Min base quality |
| `--min-mapping-quality` | `0` | int | 最小比对质量｜Min mapping quality |
| `--max-depth` | `0` | int | 最大深度限制(0=无限)｜Max depth (0=unlimited) |
| `--window-size` | `1000000` | int | 窗口大小｜Window size bp |
| `--step-size` | `100000` | int | 步长｜Step size bp |
| `--minibwa-path` | `~/software/minibwa/minibwa` |  | minibwa二进制路径｜minibwa binary path |
| `--samtools-path` | `~/.local/bin/samtools` |  | samtools二进制路径｜samtools binary path |
| `--resume` | — |  | 断点续传｜Resume (skip completed samples) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `-i, --input` | 必填 |  | FASTQ输入目录｜FASTQ input directory |
| `-p, --pattern` | `_1.fq.gz` |  | R1匹配模式｜R1 pattern |
| `-o, --output-dir` | `./minibwa_output` |  | 输出目录｜Output directory |
| `--mode` | — |  | 比对模式｜Alignment mode: standard=短读PE/SE, hic=Hi-C, meth=BS-seq, long=长读 |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--preset` | `adap` |  | -x 预设: sr/lr/adap｜preset |
| `-k, --min-seed` | `19` | int | -k 最小种子长度｜min seed length |
| `-c, --max-occ` | `250` | int | -c 最大种子出现次数｜max seed occurrences |
| `--max-gap` | `100` | int | minibwa -g 最大gap｜max gap size |
| `-w, --bandwidth` | `100` | int | -w 带宽｜bandwidth |
| `-W, --bandwidth-long` | `20000` | int | -W 长读带宽｜long bandwidth |
| `-m, --min-chain-score` | `25` | int | minibwa -m 最小链接分数｜min chaining score |
| `--sec-ratio` | `0.5` | float | minibwa -p 次要/主要得分比｜secondary-to-primary ratio |
| `-N, --max-sec` | `50` | int | -N 保留次要比对数｜retain N secondary alignments |
| `-s, --min-dp-score` | `30` | int | -s 最小DP得分*{-A}｜min DP score |
| `-A, --match-score` | `2` | int | -A 匹配得分｜matching score |
| `-B, --mismatch-penalty` | `8` | int | -B 错配罚分｜mismatch penalty |
| `-O, --gap-open` | `12,23` |  | -O gap开放罚分｜gap open penalty |
| `-E, --gap-ext` | `2,1` |  | -E gap延伸罚分｜gap extension penalty |
| `-R, --read-group` | — |  | -R SAM read group line｜read group line |
| `-u, --no-unmap` | — | store_true | -u 不输出未比对read｜do not output unmapped |
| `--outn` | `0` | int | --outn 输出次要比对上限｜max secondary to output |
| `-y, --copy-comment` | — | store_true | -y 复制FASTA/Q注释｜copy comments |
| `-Y, --soft-clip-supp` | — | store_true | -Y 软剪切补充比对｜soft clip supplementary |
| `-K, --batch-size` | `100m,1g` |  | -K 批处理大小｜batch size |
| `--markdup` | — | store_true | 标记重复｜Mark duplicates |
| `--remove-dup` | — | store_true | 移除重复（需--markdup）｜Remove duplicates (requires --markdup) |
| `--skip-coverage` | — | store_true | 跳过覆盖度分析｜Skip coverage analysis |
| `--min-base-quality` | `0` | int | 最小碱基质量｜Min base quality |
| `--min-mapping-quality` | `0` | int | 最小比对质量｜Min mapping quality |
| `--max-depth` | `0` | int | 最大深度限制(0=无限)｜Max depth (0=unlimited) |
| `--window-size` | `1000000` | int | 窗口大小｜Window size bp |
| `--step-size` | `100000` | int | 步长｜Step size bp |
| `--minibwa-path` | `~/software/minibwa/minibwa` |  | minibwa二进制路径｜minibwa binary path |
| `--samtools-path` | `~/.local/bin/samtools` |  | samtools二进制路径｜samtools binary path |
| `--resume` | — | store_true | 断点续传（跳过已完成样品）｜Resume (skip completed samples) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- minibwa（默认 ~/software/minibwa/minibwa，可用 MINIBWA_PATH 环境变量覆盖）
- samtools（默认 ~/.local/bin/samtools，可用 SAMTOOLS_PATH 覆盖）
- 两者路径支持 get_tool_path 三级解析（环境变量 > 配置文件 > 默认）

## 常见问题 | FAQ

**Q1：断点续传怎么用？**
加 --resume 后，重跑会检查每个样本是否已有 <样本>.minibwa.bam 和 .bam.bai，都有则跳过该样本。索引也一样：已存在 .l2b/.mbw 就不再重建。

**Q2：BS-seq（meth）模式为什么多建一次索引？**
BS-seq 比对需要标准索引再加一个 .meth.mbw 索引，程序会自动检测并按需补建，无需手动处理。

**Q3：--remove-dup 为什么要求 --markdup？**
移除重复必须先标记重复（markdup 计算哪些是重复），所以 --remove-dup 需与 --markdup 同时指定，否则报错。

**Q4：提示找不到 R2，按单端处理正常吗？**
若某样本只有 R1、没有配对 R2，程序会警告并把它当单端处理，这是容错行为；请确认是不是 R2 文件名与 -p 模式不匹配。

**Q5：换比对参数重跑为什么结果没变？**
索引和已完成的样本会被复用。换 --preset/-k 等比对参数后，需删除旧索引和旧 BAM（或换输出目录）重跑，否则会沿用旧结果。
