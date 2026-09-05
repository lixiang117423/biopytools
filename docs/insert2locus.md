# insert2locus — 转基因插入位点提取 | Transgenic Insertion Locus Extraction

一句话理解：**把「外源片段插到基因组哪个位置」这件事自动找出来**——从转基因材料的二代测序数据出发，围绕插入序列（载体 + 插入片段）反复「钓鱼」向外扩展，最终重构出完整的插入位点序列（左边基因组侧翼 + 插入片段 + 右边侧翼），并用测序覆盖验证结果可信度，输出一张带分级的 HTML 报告。

## 功能概述 | Overview

- 用 soft-clip 与 mate-unmapped 两种信号「钓鱼」，抓取跨越插入边界的 reads
- 迭代步移（walking），从已知序列一路向外招募新 reads，像「摸着石头过河」般扩展未知侧翼
- 用 SPAdes 重组装，配合 T-DNA 作骨架拼出完整 locus（含内部串联重复）
- 把 WGS 重新比回重构出的 locus，做覆盖验证并给出 PASS / WARN / FAIL 三级分级
- 单一整合 HTML 报告，多样本自动带顶部导航
- 断点续传（按步骤跳过已完成部分），单样本失败不中断其余样本

## 快速开始 | Quick Start

```bash
biopytools insert2locus -i fq_dir/ -f insert.fasta -o output/
```

最小输入：一个 fastq 目录（或单个 R1 文件）+ 一个插入序列 fasta（载体 + 片段的完整构建）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表，后面的参数说明都会用到：

| 术语 | 通俗理解 |
|------|----------|
| 转基因 / 插入片段 | 人为转入生物体的外源 DNA，如「载体 + T-DNA」 |
| 构建 (construct) | 载体骨架 + 插入片段连在一起的完整人工序列，像「拼装好的一整段说明书」 |
| soft-clip | 一条 read 比对到参考时「一头对得上、另一头对不上」，那半截悬空部分；常是跨边界的信号 |
| 步移 (walking) | 从已知序列出发，反复招募新 reads 向外延伸，像「踩着已知的石头找下一块」 |
| LB / RB 边界 | 插入片段与基因组左右相连的两个「接口」 |
| locus | 完整插入位点 = 左侧基因组侧翼 + 插入片段 + 右侧侧翼 |
| MAPQ | 比对质量分，越高越可信；0 表示多位置都能比上 |
| 覆盖深度 | 某个位置被多少条 reads 盖住，越高越有把握 |

## 输入 | Input

### 测序数据（-i, --input）

fastq 目录或单个 R1 文件。目录模式下自动识别配对样本，支持的配对后缀（按优先级）：
`_1_clean.fq.gz` / `_2_clean.fq.gz`，其次 `_1.fq.gz` / `_2.fq.gz`。

### 插入序列（-f, --insert-fasta）

完整构建（载体 + 插入片段）的 fasta，作为比对「伪基因组」和步移的起点。单条或多条 record 均可。

```text
>construct
ACGT...（载体 + T-DNA 的完整序列）
```

### T-DNA 序列（--tdna-fasta，可选）

单独的插入片段（T-DNA）fasta。给了它，程序能区分「插入区」与「载体骨架」，骨架区覆盖度会被单独统计，拼装时也会把 T-DNA 当作可信骨架穿过内部串联重复。不给则把 `-f` 整体当作插入片段。

## 参数说明 | Parameters

### 必需输入 | Required

**通俗理解|In plain words:** 三个必填项：测序数据、插入序列、输出目录。插入序列是「起点」，测序数据是「证据」，输出目录是「归宿」。没有可省的。

### 样本与资源 | Sample & resources

**通俗理解|In plain words:** 这几个管「跑多快、吃多少内存」。线程数越大越快但越占 CPU；`--sort-mem` 是排序时给每个线程的内存，数据大时可适当调大；`--read1-suffix` 是当你的文件命名不是标准后缀时手动指定 R1 后缀（一般自动检测就够，不用动）。

### 步移参数 | Walking parameters

**通俗理解|In plain words:** 步移是「从已知序列往外一层层扩展」的过程。`--max-rounds` 是扩展最多走几轮；`--min-softclip` 和 `--min-unmapped` 是「什么样的悬空信号才配当诱饵」的最短长度门槛（太短的多是噪声）；`--min-growth` 是「这轮新长了多少才算有进展」的收敛判据（太小会让它过早停，太大又会空转）；`--mapq-min` 是招募 reads 的最低质量分；`--repeat-cap` 是单轮新招 reads 的上限（撞上基因组重复区时防止爆炸）。**绝大多数情况用默认值即可，只在步移过早收敛或撞重复区时才需要微调。**

### 边界与报告 | Junction & report

**通俗理解|In plain words:** `--junction-flank` 是边界报告要求的最短侧翼长度，侧翼太短就降级为 WARN；`--target-flank` 是 LB/RB 想走多远的目标（默认不限，能走多远走多远，靠自然收敛刹住）。

### 运行控制 | Execution control

**通俗理解|In plain words:** `--force` 忽略断点全部重跑（换参数重跑时用）；`--log-level` 是日志详细程度，一般不用动。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 六步流水线：先把测序数据比到插入序列上找「半截比对」的信号，再据此钓出跨边界的 reads，然后一步步向外扩展，拼出完整插入位点，最后用全基因组数据回验一遍并打分。

```text
输入 fastq + 插入序列
    │
    ▼
01 比对：WGS 比对到 insert 伪基因组（bwa + samtools）
    │
    ▼
02 junction 钓取：soft-clip + mate-unmapped 抓跨边界 reads
    │
    ▼
03 迭代步移：逐轮招募新 reads 向外扩展（SPAdes 组装 + 再比对）
    │
    ▼
04 完整 locus：招募全部覆盖 reads 重组装，锚定 LB/insert/RB
    │
    ▼
05 覆盖验证：WGS 比回 locus，分段覆盖 + 跨界 reads 分级
    │
    ▼
06 单一 HTML 报告（多样本带导航）
```

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml            # 工具版本与运行参数
├── 01_mapping/
│   ├── {sample}_vs_insert_sorted.bam    # 比对到插入序列的 BAM
│   ├── {sample}_vs_insert_flagstat.txt  # 比对统计
│   └── {sample}_insert_coverage.tsv     # 逐 record / 分区覆盖统计
├── 02_junction_reads/
│   ├── {sample}_softclip.fastq          # 悬空半截 reads
│   ├── {sample}_mate_unmapped.bam       # 另一端未比对的 reads
│   └── {sample}_flank_candidates_R{1,2}.fastq  # 候选侧翼 reads
├── 03_walking/
│   └── rounds/{sample}/                 # 步移中间文件（recruited/bait 等）
├── 04_locus/
│   ├── {sample}_contigs.fasta           # SPAdes 组装出的 contigs
│   ├── {sample}_junction_contigs.fasta  # 跨界 contig 序列
│   ├── {sample}_junction_report.tsv     # 边界报告（核心）
│   └── {sample}_complete_locus.fasta    # 完整 locus 序列（核心）
├── 05_verify/
│   ├── {sample}_verification_summary.tsv  # 分级结果（核心）
│   ├── {sample}_coverage.tsv            # LB/insert/RB 分段覆盖
│   └── {sample}_vs_locus_sorted.bam     # 比回 locus 的 BAM
├── 99_logs/
│   └── insert2locus_YYYYMMDD_HHMMSS.log
└── insert2locus_report.html             # 整合 HTML 报告（核心）
```

### 关键文件说明 | Key files

- `{sample}_complete_locus.fasta`：重构出的完整插入位点，FASTA 头里带 `LB=xxxbp insert=xxxbp RB=xxxbp` 三段长度，是最核心的结果
- `{sample}_junction_report.tsv`：边界报告，列出每个跨界 contig 的 border（L/R/LR/unmapped）、在插入序列上的位置、锚定长度与侧翼长度
- `{sample}_verification_summary.tsv`：分级（grade）、左右跨界 reads 数、LB/insert/RB 三段平均深度、左右侧翼是否为植物来源
- `insert2locus_report.html`：所有样本汇总的可视化报告，含分级徽章、locus 结构示意图、序列卡片

## 结果解读 | Interpreting Results

### 分级（grade）| Grade

**通俗理解|In plain words:** 相当于一份「可信度体检报告」。PASS 是「证据充分」，WARN 是「有结果但某个证据略弱」，FAIL 是「没拼出完整位点」。

| 级别<br>Grade | 含义<br>Meaning |
|------|------|
| PASS | 三段全覆盖、左右侧翼足够长、左右跨界 reads 都有、侧翼不是载体来源 |
| WARN | 任一段有零覆盖窗口 / 侧翼不足 50bp / 某侧跨界 reads 为 0 / 侧翼匹配到构建（载体骨架来源） |
| FAIL | 未重构出完整 locus |

- 落到 WARN 的原因会体现在 `verification_summary.tsv` 各列里：`lb_depth / insert_depth / rb_depth` 看哪段覆盖不足，`lb_junction_reads / rb_junction_reads` 看哪侧跨界证据缺
- `lb_plant / rb_plant` 为 False 表示该侧翼匹配到了构建序列（当年 LB 误判的坑），提示侧翼其实是载体骨架而非植物基因组

### 覆盖统计（coverage.tsv）| Coverage

每段（LB / insert / RB）给出 start、end、length、mean_depth、min_depth、zero_windows。`zero_windows > 0` 表示该段存在完全没被测到的窗口，会触发 WARN。

## 参数选择建议 | Parameter Guidance

- **默认参数直接跑**：大多数样本用默认值即可拿到可靠结果
- **步移过早收敛**（拼出的侧翼很短却标成完成）：适当提高 `--max-rounds`，或降低 `--min-growth`（比如 50 降到 20）
- **撞上重复区**（新招 reads 暴涨却无进展）：调低 `--repeat-cap` 或调高 `--min-softclip` / `--min-unmapped` 收紧诱饵质量
- **有单独 T-DNA 序列**：务必给 `--tdna-fasta`，能区分插入区与骨架、拼通内部串联重复，显著提升成功率
- **换参数重跑**：加 `--force`，否则断点续传会复用旧结果

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq目录或R1文件(自动识别样本)｜fastq dir or R1 file (auto-detect samples) |
| `-f, --insert-fasta` | 必填 |  | 插入序列fasta(载体+片段)｜Insert fasta (vector+fragment) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--sort-mem` | `2G` |  | samtools sort单线程内存｜samtools sort per-thread memory |
| `--read1-suffix` | — |  | R1后缀(默认自动检测)｜R1 suffix (auto-detect by default) |
| `--max-rounds` | `30` | int | 步移最大轮数｜Max walking rounds |
| `--min-softclip` | `25` | int | 诱饵softclip最短长度｜Min softclip for bait |
| `--min-unmapped` | `400` | int | 未比对contig作诱饵最短长度｜Min unmapped length for bait |
| `--min-growth` | `50` | int | 收敛判定增量阈值｜Growth threshold for convergence |
| `--mapq-min` | `1` | int | 招募最低MAPQ｜Min MAPQ for recruitment |
| `--repeat-cap` | `10000` | int | 单轮新招reads上限(撞重复区)｜Per-round new-read cap |
| `--junction-flank` | `50` | int | 边界报告最短侧翼｜Min flank for border report |
| `--tdna-fasta` | — |  | 单独插入序列fasta(区分insert与载体骨架)｜Standalone T-DNA fasta |
| `--target-flank` | — | int | LB/RB目标侧翼长度(默认不限,尽可能走远)｜Target flank length (default: walk as far as possible) |
| `--force` | `False` |  | 忽略断点全部重跑｜Ignore checkpoints and rerun |
| `--log-level` | `INFO` |  | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq目录或R1文件｜fastq dir or R1 file |
| `-f, --insert-fasta` | 必填 |  | 插入序列fasta(载体+片段)｜Insert fasta (vector+fragment) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `-t, --threads` | `12` | int | 线程数｜Threads (default 12) |
| `--sort-mem` | `2G` |  | samtools sort单线程内存｜per-thread sort memory |
| `--read1-suffix` | — |  | R1后缀(默认自动检测)｜R1 suffix (auto-detect by default) |
| `--max-rounds` | `30` | int | 步移最大轮数｜Max walking rounds |
| `--min-softclip` | `25` | int | 诱饵softclip最短长度｜Min softclip for bait |
| `--min-unmapped` | `400` | int | 未比对contig作诱饵最短长度｜Min unmapped length for bait |
| `--min-growth` | `50` | int | 收敛判定增量阈值｜Growth threshold for convergence |
| `--mapq-min` | `1` | int | 招募最低MAPQ｜Min MAPQ for recruitment |
| `--repeat-cap` | `10000` | int | 单轮新招reads上限(撞重复区)｜Per-round new-read cap |
| `--junction-flank` | `50` | int | 边界报告最短侧翼｜Min flank for border report |
| `--tdna-fasta` | — |  | 单独插入序列fasta(区分insert与载体骨架;不给则-f整体当insert)｜Standalone T-DNA fasta (separates insert from backbone) |
| `--target-flank` | — | int | LB/RB目标侧翼长度(默认不限,尽可能走远到自然收敛;到达后小增量即收敛)｜Target flank length (default: walk as far as possible) |
| `--force` | — | store_true | 忽略断点全部重跑｜Ignore checkpoints and rerun |
| `--log-level` | `INFO` |  | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- bwa（`~/miniforge3/envs/align/bin/bwa`，可由 `BWA_PATH` 覆盖）
- samtools（`~/miniforge3/envs/align/bin/samtools`）
- seqkit（`~/miniforge3/envs/misc/bin/seqkit`）
- SPAdes（`~/miniforge3/envs/asm/bin/spades.py`）

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
断点续传按输出文件存在性判断。换参数后要加 `--force`，否则会跳过已完成的步骤复用旧结果。

**Q2：为什么不能用 MAPQ 过滤掉 0 分 reads？**
拼装时内部串联重复区必须靠 MAPQ0 的多比对 reads 才能贯穿。代码实测：用全部 3911 条 mapped reads 才拼出贯穿 contig，只保留 MAPQ≥1 的 2025 条拼不出。所以招募阶段故意不加 `-q` 过滤。

**Q3：WARN 里说「侧翼匹配构建」是什么意思？**
程序把 locus 两端侧翼重新比回构建序列，若匹配 ≥50bp 说明该侧翼其实来自载体骨架而不是植物基因组（当年的 LB 误判坑），此时会降级为 WARN 提醒人工复核。

**Q4：拼出的 contig 从中间断成两半怎么办？**
通常是插入片段内部的串联重复把组装劈开了。给 `--tdna-fasta` 让程序把 T-DNA 当作可信骨架（trusted-contigs），就能跨过重复区拼通。

**Q5：多样本中某个样本失败了怎么办？**
单样本失败不中断其余样本，失败的样本会在报告里标为 FAIL，其余样本正常出结果。
