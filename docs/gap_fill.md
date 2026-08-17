# Gap 填充 | Gap Filling (TGS-GapCloser)

一句话理解：**用三代长读长测序数据（ONT/PacBio/HiFi）把基因组组装里那些「NNNN…」的空洞（gap）填上，让基因组更连续、更完整**。

## 功能概述 | Overview

- 封装 TGS-GapCloser2，用长读段跨 gap 比对来填充 N 区段
- 支持三种 TGS 类型：ont / pb / hifi（自动匹配对应参数与同一性阈值）
- 可选纠错模式：none（不纠错，默认）/ racon（长读纠错）/ pilon（NGS 短读精修）
- 内置第 2 轮 quarTeT gapfiller：用 hifiasm unitig/contig 填第 1 轮没填上的残余 gap
- 支持断点续传（-f 强制重跑）
- 输出填充前后对比报告 + 逐个 gap 明细表

## 快速开始 | Quick Start

```bash
biopytools gap-fill -s scaffolds.fa -t ont -ir ont_reads.fq.gz -o gapclosed
```

最小输入：一个 scaffold 级 FASTA + 一个三代长读段文件（FASTA/FASTQ，可 .gz）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| gap / N 区段 | 组装时没测通、用一串 N 占位的「空洞」，基因组里像打了补丁的坑 |
| scaffold | 已经排好相对位置、但中间还有空洞的序列（比 contig 更完整，比染色体还差一步） |
| 三代测序(TGS) | 一次能读几 kb 到上百 kb 的长读段，能横跨 gap 两端，看清空洞里是什么 |
| unitig | hifiasm 组装中间产物，一段高可信的连续序列，用来给 gap「补料」 |
| 纠错 | 把读段里的错误碱基改对；racon 用长读纠错，pilon 用 NGS 短读精修 |
| flanking 序列 | gap 两侧的「扶手」，用来锚定填充序列该接在哪 |

## 输入 | Input

- **scaffold 文件**：FASTA，含 N 区段（gap）的组装序列；
- **TGS reads 文件**：FASTA 或 FASTQ，支持 .gz 压缩；
- **unitig 文件（可选）**：hifiasm 的 unitig/contig（FASTA 或 .gfa），提供第 2 轮填充；
- **NGS reads（pilon 模式）**：-m pilon 时需另给 NGS 短读。

示例（scaffold，gap 用 N 表示）：

```text
>scaffold1
ACGTACGT...NNNNNNNNNN...TTGCAACG
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 四个缺一不可。-s 是「要填的基因组」，-ir 是「用来填的长读段」，-t 告诉程序读段是哪种三代技术（决定比对参数），-o 是输出文件名的前缀（最终产物是 前缀.gapcloser.fa）。

### 纠错模式 | Error correction mode

**通俗理解|In plain words:** -m 决定「填完要不要顺手纠错」。none 最省事、最快，一般先用它；racon 用长读把填充区域的错误改掉，质量更高；pilon 需要额外的 NGS 短读和 Java/Pilon 环境，最慢。**追求速度用 none，追求精度再上 racon/pilon。**

### 软件路径 | Software paths

**通俗理解|In plain words:** 这些路径程序都会自动探测（环境变量 > 配置文件 > 默认路径），一般不用手动指定。只有软件装在不常见位置时才需要显式传。

### 过滤参数 | Filtering parameters

**通俗理解|In plain words:** -idy（最小同一性）和 -l（最小匹配长度）控制「多像、多长的读段才用来填」。**不传时会按 tgstype 自动设**：ont 是 0.3 / 300，pb 和 hifi 是 0.2 / 200。一般不用动。

### 性能参数 | Performance

**通俗理解|In plain words:** -threads 线程数越大越快；-chunk 是把输入切成几块并行处理，块越多越并行但也越吃内存，默认 3 通常够用。

### 高级筛选 | Advanced filtering

**通俗理解|In plain words:** 这些是给「填充太激进或太保守」时微调用的：-min-nread/-max-nread 限制每个 gap 用多少条读段；-max-candidate 限制候选数；-g-check 开启 gap 大小差异检查；-minmap-arg 传自定义 minimap2 参数。**绝大多数情况不需要动。**

### 第 2 轮填充（quarTeT）| Round-2 filling

**通俗理解|In plain words:** 第 1 轮填完还有残余 gap 时，用 -ug 传入 hifiasm unitig/contig 做第 2 轮补充填充。-fl 是 gap 两侧取多长的「扶手」，-al/-ai 是第 2 轮比对的最短长度/最低同一性，-mfl 限制单次最多填多长，-mgl 是第 2 轮只处理多长的 gap。**不传 -ug 就只做第 1 轮。**

### 流程控制 | Pipeline control

**通俗理解|In plain words:** -f 忽略断点续传、从头强制重跑；--dry-run 只打印会执行的命令、不真跑（用来检查命令对不对）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先跑 TGS-GapCloser2 填一遍，再看还剩多少空洞；如果给了 unitig 且还有残余空洞，就用 quarTeT 再补一轮。

```text
scaffold FASTA + TGS reads
    │
    ▼
第1轮: TGS-GapCloser2 填充 → {prefix}.gapcloser.fa
    │
    ▼
检查剩余 gap（是否还有 >= min_gap_length 的 N 区段）
    │
    ├─ 无 gap 或未给 -ug → 结束
    │
    └─ 有 gap 且给了 -ug → 第2轮: quarTeT gapfiller
          ├─ 提取 gap 两侧 flanking
          ├─ minimap2 把 unitig 比对到 flanking
          ├─ 找同时跨两侧的 unitig 填充
          └─ 生成最终 {prefix}.gapcloser.fa
    │
    ▼
生成 gap 前后对比报告 + 明细表
```

## 输出 | Output

```text
（输出到 -o 前缀所在目录）
├── gapclosed.gapcloser.fa      # 填充后的最终基因组（核心产物）
├── gapclosed.gap_report.txt    # 填充前后 gap 对比报告
├── gapclosed.gap_table.tsv     # 逐个 gap 明细（填充状态/残留长度）
├── gapclosed.gapcloser.fa.round1  # 第1轮结果备份（仅跑第2轮时生成）
└── tgsgapcloser.log            # 运行日志
```

注：文件名前缀随 -o/--output-prefix 改变，上面以 -o gapclosed 为例。

## 结果解读 | Interpreting Results

### 1. 填充后的基因组（{prefix}.gapcloser.fa）

**通俗理解|In plain words:** 这是最终要用的文件——原来的 N 空洞被真实序列替换了（或缩短了）。可以重新跑一遍 gap-stat 或 seqkit stats 对比 N 比例。

### 2. gap 对比报告（{prefix}.gap_report.txt）

**通俗理解|In plain words:** 一份「填洞前后对比成绩单」。

- 处理前/后：总 gap 数、总 gap 长度、平均/最大 gap 长度、长度分布；
- 填充效果：填充了多少个 gap、多少碱基、按数量的填充率（%）；
- **填充率越高越好**；剩余 gap 越少、越短越好。

### 3. gap 明细表（{prefix}.gap_table.tsv）

**通俗理解|In plain words:** 逐个 gap 的「命运清单」，每个 gap 一行。

- status 列：Filled（已填）/ Remaining（残留）/ Unknown（无法定位）；
- 可快速定位哪些 gap 没填上、残留在哪条序列哪段坐标。

## 参数选择建议 | Parameter Guidance

- **-t, --tgstype**：按你的读段选——ONT 用 ont，PacBio CLR 用 pb，PacBio HiFi 用 hifi
- **-m, --mode**：追求速度用默认 none；填充质量要求高且是长读用 racon；有 NGS 短读可再上 pilon
- **-ug, --unitig-file**：第 1 轮填完还有不少残余 gap 时，提供 hifiasm 的 unitig 文件做第 2 轮，通常能再降一批 gap
- **-threads**：核多就调大，加速明显

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --scaff-file` | 必填 |  | 输入scaffold文件｜Input scaffold file |
| `-t, --tgstype` | 必填 | ont/pb/hifi | TGS类型｜TGS type (ont/pb/hifi) |
| `-ir, --reads-file` | 必填 |  | 输入TGS reads文件｜Input TGS reads file |
| `-o, --output-prefix` | 必填 |  | 输出前缀｜Output prefix |
| `-m, --mode` | `none` | none/racon/pilon | 纠错模式｜Error correction mode (default: none) |
| `--tgsgapcloser-path` | — |  | TGS-GapCloser路径｜TGS-GapCloser path (default: auto-detect) |
| `-idy, --min-idy` | — | float | 最小同一性｜Min identity (auto-set) |
| `-l, --min-match` | — | int | 最小匹配长度｜Min match length (auto-set) |
| `-threads, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-chunk` | `3` | int | 分块数量｜Chunk count (default: 3) |
| `-g-check` | — |  | 启用Gap大小差异检查｜Enable gap size difference check |
| `-min-nread` | `1` | int | 最小reads数量｜Min read count (default: 1) |
| `-max-nread` | `-1` | int | 最大reads数量｜Max read count (default: -1) |
| `-max-candidate` | `200` | int | 最大候选数｜Max candidates (default: 200) |
| `-racon, --racon-path` | — |  | Racon路径｜Racon path |
| `-racon-round` | `3` | int | Racon轮数｜Racon rounds (default: 3) |
| `-pilon, --pilon-path` | — |  | Pilon路径｜Pilon path |
| `-ngs, --ngs-file` | — |  | NGS reads文件｜NGS reads file |
| `-java, --java-path` | — |  | Java路径｜Java path |
| `-samtools, --samtools-path` | — |  | Samtools路径｜Samtools path |
| `-pilon-mem` | `300G` |  | Pilon内存｜Pilon memory (default: 300G) |
| `-pilon-round` | `3` | int | Pilon轮数｜Pilon rounds (default: 3) |
| `-minmap-arg` | — |  | 自定义minimap2参数｜Custom minimap2 arguments |
| `-ug, --unitig-file` | — |  | hifiasm unitig/contig文件（第2轮填充）｜hifiasm unitig/contig file (round 2) |
| `-fl, --flanking-len` | `5000` | int | Flanking序列长度（bp）｜Flanking sequence length (bp) (default: 5000) |
| `-al, --min-align-len` | `1000` | int | 最小比对长度（bp）｜Min alignment length (bp) (default: 1000) |
| `-ai, --min-identity` | `40` | int | 最小比对同一性（%）｜Min alignment identity (%) (default: 40) |
| `-mfl, --max-filling-len` | `1000000` | int | 最大填充长度（bp）｜Max filling length (bp) (default: 1000000) |
| `-mgl, --min-gap-length` | `100` | int | 第2轮识别/填充的最小gap长度(bp)｜Min gap length (bp) for round-2 (default: 100) |
| `-f, --force` | — |  | 忽略断点续传强制重跑｜Force rerun, ignore checkpoint |
| `--dry-run` | — |  | 只打印命令不执行｜Dry run, print commands only |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- TGS-GapCloser2（tgsgapcloser2，默认自动探测 ~/software/TGS-GapCloser2/.../tgsgapcloser2，可用环境变量 TGSGAPCLOSER_PATH 覆盖）
- minimap2（仅第 2 轮 quarTeT 填充需要）
- racon（仅 -m racon）
- pilon + Java + samtools + NGS reads（仅 -m pilon）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。第 1 轮输出 {prefix}.gapcloser.fa 已存在时会跳过第 1 轮；第 1 轮 + 第 2 轮都完成（存在 round1 备份）时会整体跳过。想强制重跑用 -f。

**Q2：-idy/-l 不传会怎样？**
会按 tgstype 自动设置：ont → 同一性 0.3、匹配长度 300；pb/hifi → 同一性 0.2、匹配长度 200。通常无需手动指定。

**Q3：第 2 轮填充一定要吗？**
不是。只有同时满足「第 1 轮后仍有剩余 gap」和「提供了 -ug unitig 文件」才会跑第 2 轮。没有 unitig 就只做第 1 轮。

**Q4：--dry-run 有什么用？**
只打印即将执行的命令、不真跑，适合第一次用或换参数前先核对命令是否正确。

**Q5：为什么运行前会删 done_step* 文件？**
TGS-GapCloser 会用 done_step* 文件做内部断点，程序会在运行前清理它们、避免并发或上次残留影响本次填充，属正常行为。

