# XP-CLR 跨群体选择信号扫描 | XP-CLR Cross-Population Selection Scan

一句话理解：**比较两个群体在每个基因组位置的基因频率差异，给全基因组逐段打分，找出其中一个群体可能正在被自然选择/人工选择"改造"的基因组区域**。

## 功能概述 | Overview

- 输入一张 VCF 变异表 + 两份样本名单（群体A=疑似受选择的群体，群体B=参照群体），输出全基因组每个窗口的**选择信号分数**和一张**Top 候选区域排行榜**
- 全流程：读 VCF 头部校验样本 → **逐染色体**调用 xpclr 计算（支持断点续传、单条染色体失败不拖垮整体）→ 合并成全基因组表 → 按分数排出 Top 候选窗口 → 记录软件版本
- XP-CLR 是选择信号分析的经典方法（Chen et al. 2010），特别适合**两个亲缘群体**（如驯化种 vs 野生种、抗病品系 vs 感病品系）之间的受选择区域定位
- 核心产出：`02_merged/*.xpclr.genome.tsv`（全基因组逐窗口分数表）和 `03_top/*.xpclr.top50.tsv`（最值得优先关注的候选窗口）

## 快速开始 | Quick Start

```bash
biopytools xpclr -i pop.vcf.gz -a popA.txt -b popB.txt -o out_dir/
```

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>In plain words |
|---|---|
| XP-CLR | 给"两地口音差异"打分的方法：如果某段基因组里群体A和群体B的基因频率差异特别大，就像两个村子的人说话口音突然大不一样，说明这段区域可能在其中一个群体里被选择压力"带着跑"了 |
| 群体A / 群体B | 群体A=你怀疑受选择的群体（如驯化品种、抗病品系）；群体B=参照群体（如野生近缘、感病品系）。分数越高=A 相对 B 在这段越"特殊" |
| 选择信号(selective sweep) | 自然/人工选择"清扫"某段基因组的痕迹：被青睐的基因版本在群体里迅速铺开，把旁边一整片基因都带成了同一个版本 |
| 连锁不平衡(LD) | "邻居消息重复"：基因组上挨得近的位点往往一起遗传，消息互相重复。XP-CLR 计算时会考虑这种重复，把冗余信息的权重降下来，避免"一段话数成十票" |
| 窗口(window) | 把染色体切成一段一段（默认每段 20kb），逐段打分。像把一条长路按 20 米一段画格子，逐格检查路况 |
| XP-CLR 分数(xpclr) | 每个窗口的原始得分，来自"假设受选择的模型"比"假设没受选择的模型"好多少，越大越像受选择 |
| 标准化分数(xpclr_norm) | 把原始分换算成"全班排名"式的标准分（z-score），消除不同染色体之间分数尺度的差异，跨染色体比较时看它 |
| bgzip / tabix | bgzip 是给 VCF 文件"打包压缩"的格式，tabix 是给压缩包建"目录索引"，有了索引程序才能只抽取指定染色体的数据而不必通读全文件 |

## 输入 | Input

- `-i` VCF 文件：必须是 **bgzip 压缩**且带 **tabix 索引**（`.tbi` 或 `.csi`）。若是普通文本 VCF，先转换：

```bash
bgzip -c in.vcf > in.vcf.gz && tabix -p vcf in.vcf.gz
```

- `-a` / `-b` 样本名单：纯文本，**每行一个样本 ID**（须与 VCF 最后一行表头的样本名完全一致，含大小写）；空行自动忽略。两份名单**不能有同一样本**。

```text
# popA.txt 示例|example
SAMPLE_01
SAMPLE_02
SAMPLE_03
```

- 可选 `--chroms`：逗号分隔，只跑指定染色体（如 `--chroms 3L,3R`）；不给则跑 VCF 里全部染色体。
- 样本列表**路径里不能含逗号**（xpclr 工具用逗号区分"文件"和"内联名单"，会误解析）。

## 参数说明 | Parameters

**通俗理解|In plain words:** 绝大多数参数**一般不用动**——窗口、LD、重组率的默认值就是 xpclr 工具的官方默认，适合大多数重测序数据。真正值得看的是 `--size`/`--step`（数据密度极端时才调）和 `--top-n`（想要更长/更短的候选清单时才调）。

**输入与输出组|Input & output：** 管"喂什么、产到哪"。`-i/-a/-b` 三个必填项决定分析对象；`-o` 决定结果目录（默认在当前目录建 `xpclr_out/`）；`--label` 只是结果文件名前缀，改它纯粹为了好认，不影响任何计算。`--chroms` 调大调小只影响跑多少条染色体、跑多久，不影响每条染色体内部的结果。

**窗口组|Window：** 管"格子画多大、每隔几米画一个"。`--size` 是格子宽度，`--step` 是每次往前挪的距离——`step < size` 时格子会重叠（看得细、算得慢），默认两者相等即不重叠。`--maxsnps`/`--minsnps` 是每个格子里 SNP 数的上下限：太多会把远处无关位点混进来，太少则算不出可靠分数。一般都不用动。

**模型组|Model：** 管"打分时对重复信息的态度"。`--ld` 是 LD 截断，越接近 1 越"宽容"地保留重复信息；`--phased` 只在数据已做单倍型定相时打开（r2 估计更准）；`--rrate` 是每碱基重组率，动它的情况极少。默认值即文献常规用法。

**汇总与运行组|Summary & runtime：** `--top-n` 决定候选排行榜取前几名（默认 50）；`--xpclr-path`/`--log-level` 属于环境调试项，除非换环境或排查问题，一般不动。

完整参数速查表见下方自动生成区块|Full parameter reference: see the auto-generated block below.

## 分析流程 | Pipeline

```text
VCF(.gz+.tbi) + popA.txt + popB.txt
        |
        v
[1] 校验: 样本在VCF中存在 / A与B不重叠 / 染色体列表合法
        |
        v
[2] 逐染色体循环: xpclr --chr <每条一次>   (已完成染色体自动跳过=断点续传;
        |                                 单条失败记WARNING继续跑其余)
        v
[3] 合并: 所有成功染色体 -> 02_merged/<label>.xpclr.genome.tsv
        |
        v
[4] Top候选: 按xpclr_norm降序取前N -> 03_top/<label>.xpclr.top<N>.tsv
        |
        v
[5] 记录: 00_pipeline_info/software_versions.yml + 99_logs/xpclr.log
```

## 输出 | Output

```text
out_dir/
├── 00_pipeline_info/
│   └── software_versions.yml     # 软件版本+关键参数快照
├── 01_xpclr/
│   └── <chrom>.xpclr.tsv         # 每条染色体一个(断点续传的粒度)
├── 02_merged/
│   └── <label>.xpclr.genome.tsv  # 全基因组逐窗口合并表(主结果)
├── 03_top/
│   └── <label>.xpclr.top<N>.tsv  # Top N 候选窗口
└── 99_logs/
    └── xpclr.log                 # 运行日志(含每条命令完整记录)
```

`01_xpclr/` 与 `02_merged/` 表格列含义相同（来自 xpclr 工具）：

| 列<br>Column | 含义<br>Meaning |
|---|---|
| `id` | 窗口唯一标识 `<chrom>_<start>_<stop>` |
| `chrom` / `start` / `stop` | 染色体与名义窗口坐标（bp） |
| `pos_start` / `pos_stop` | 实际落到 SNP 上的窗口边界（真实数据范围） |
| `modelL` / `nullL` | 受选择模型 / 中性模型的对数似然 |
| `sel_coef` | 选择系数估计 |
| `nSNPs` | 参与计算的 SNP 数；**过小（接近 minsnps）的窗口结论要慎用** |
| `nSNPs_avail` | 窗口内可用 SNP 数 |
| `xpclr` | 原始分数 = `2*(modelL - nullL)`，越大越像受选择 |
| `xpclr_norm` | 标准化分数（z-score）；跨染色体比较、排 Top 都以它为准 |

## 结果解读 | Interpreting Results

- **看什么：** 主要看 `03_top/` 排行榜和 `02_merged/` 表里的 `xpclr_norm` 列。数值越高=该窗口两群体差异越偏离中性预期，受选择的可能性越大。
- **经验阈值：** `xpclr_norm` 没有绝对及格线，惯例做法是取全基因组前 1%（或 Top 50）作为候选，再结合生物学背景判断。
- **好坏判据|Good vs suspicious：**
  - 好信号：候选窗口的 `nSNPs` 充足（明显大于 `minsnps`），且**相邻连续几个窗口分数都高**——成片的高分比孤峰可信得多
  - 慎报：单窗口孤立高峰、或 `nSNPs` 贴着下限的"高分"，多半是噪声或标记密度伪象
  - 若日志末尾出现"染色体失败,结果为部分基因组"WARNING，说明有染色体没算成，Top 表只覆盖了部分基因组，横向比较可能有偏

## 参数选择建议 | Parameter Guidance

- **默认 20kb 窗口适合大多数数据**：重测序群体遗传的常见密度下，20kb 内通常落 20-200 个 SNP，正是 XP-CLR 设计的舒适区。第一次跑直接用默认。
- **SNP 密度特别高**（如 WGS 深度极高或大群体）：单窗口 SNP 容易爆到 `--maxsnps` 上限之外被截断，可加大 `--size`（如 50kb）或减小 `--maxsnps`，让窗口内信息量适中。
- **SNP 密度很低**（如低深度或分型芯片）：可减小 `--size` 并把 `--minsnps` 降到 5 左右，避免大量窗口因 SNP 不足而无效。
- **样本量：每个群体 <15 个个体时噪声明显偏大**，XP-CLR 分数波动剧烈，建议优先加样本；实在不行把 `--top-n` 看成"提示清单"而非结论。
- **想看得更细：** `--step` 设为 `--size` 的一半（如 20kb 窗 10kb 步），窗口重叠能更精确定位受选择区间边界，代价是计算量约翻倍。
- **先试跑：** 大基因组想快速验证流程，用 `--chroms` 先跑一条染色体（如 `--chroms 3L`）看输出是否符合预期。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | bgzip VCF(须带.tbi/.csi索引)｜bgzipped VCF with tabix index |
| `--samples-a, -a` | 必填 |  | 群体A样本列表(每行一个ID)｜Pop A sample list (one ID per line) |
| `--samples-b, -b` | 必填 |  | 群体B样本列表文件｜Pop B sample list |
| `--output-dir, -o` | `xpclr_out` |  | 输出目录｜Output directory |
| `--label` | `popA_vs_popB` |  | 结果文件前缀｜Result file prefix |
| `--chroms` | — |  | 逗号分隔染色体列表(默认VCF全部contig)｜Comma-separated chroms (default all) |
| `--size` | `20000` | int | 窗口大小bp｜Window size bp |
| `--step` | `20000` | int | 滑窗步长bp｜Step size bp |
| `--maxsnps` | `200` | int | 窗口最大SNP数｜Max SNPs per window |
| `--minsnps` | `10` | int | 窗口最小SNP数｜Min SNPs per window |
| `--ld` | `0.95` | float | LD加权截断｜LD cutoff |
| `--phased` | `False` |  | 数据已phased(更精确r2)｜Data phased |
| `--rrate` | `1e-08` | float | 每碱基重组率｜Recombination rate per base |
| `--top-n` | `50` | int | Top候选窗口数｜Top candidate windows |
| `--xpclr-path` | `~/miniforge3/envs/selective_sweep/bin/xpclr` |  | xpclr可执行路径｜xpclr executable path |
| `--log-level` | `INFO` |  | 日志级别(DEBUG/INFO/WARNING/ERROR)｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | bgzip VCF(须带 .tbi/.csi 索引)｜bgzipped VCF with tabix index |
| `-a, --samples-a` | 必填 |  | 群体A样本列表文件(每行一个ID)｜Pop A sample list (one ID per line) |
| `-b, --samples-b` | 必填 |  | 群体B样本列表文件｜Pop B sample list |
| `-o, --output-dir` | `xpclr_out` |  | 输出目录｜Output directory (default: xpclr_out) |
| `--label` | `popA_vs_popB` |  | 结果文件前缀｜Result file prefix (default: popA_vs_popB) |
| `--chroms` | — |  | 逗号分隔染色体列表(默认VCF全部contig)｜Comma-separated chromosomes (default: all contigs in VCF) |
| `--size` | `20000` | int | 窗口大小bp｜Window size bp (default: 20000) |
| `--step` | `20000` | int | 滑窗步长bp｜Step size bp (default: 20000) |
| `--maxsnps` | `200` | int | 窗口最大SNP数｜Max SNPs per window (default: 200) |
| `--minsnps` | `10` | int | 窗口最小SNP数｜Min SNPs per window (default: 10) |
| `--ld` | `0.95` | float | LD加权截断｜LD cutoff (default: 0.95) |
| `--phased` | `False` | store_true | 数据已phased(更精确r2)｜Data phased (more precise r2) |
| `--rrate` | `1e-08` | float | 每碱基重组率｜Recombination rate per base (default: 1e-8) |
| `--top-n` | `50` | int | Top候选窗口数｜Top candidate windows (default: 50) |
| `--xpclr-path` | `~/miniforge3/envs/selective_sweep/bin/xpclr` |  | xpclr可执行路径｜xpclr executable path |
| `--log-level` | `INFO` |  | 日志级别｜Log level (default: INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **xpclr 1.1.2**（python-xpclr 发行包）：安装在 conda 环境 `selective_sweep`；因上游包与 scipy 1.18/numpy 2.x 不兼容，环境内已打兼容补丁（romberg backport、quad 结果 `.item()` 包装、size-1 数组强转），**补丁源码持久化在 `~/software/xpclr`**，重建环境时须重做补丁
- conda（工具通过 `conda run` 包装调用，自动检测环境）
- pandas / PyYAML（biopytools 侧表格合并与版本记录）

## 常见问题 | FAQ

**Q1: 报错"VCF缺少tabix索引"怎么办？**
输入必须是 bgzip 压缩 + tabix 索引的 VCF（`.vcf.gz` + `.vcf.gz.tbi`）。普通 gzip 压缩的不行，转换命令：`bgzip -c in.vcf > in.vcf.gz && tabix -p vcf in.vcf.gz`。

**Q2: 某条染色体失败了，整个流程会挂吗？**
不会。单条染色体失败记 WARNING 并跳过，其余染色体照常完成；全部染色体都失败才中止退出。结果只覆盖成功的染色体，日志末尾会有"部分基因组"WARNING 提醒。

**Q3: 群体A和B的样本能重叠吗？**
不能。同一个样本出现在两份名单里会让"差异比较"失去意义，程序会在启动时直接报错拦下。同理，名单里的样本ID必须能在 VCF 表头找到，否则也会启动报错（一次性列出所有问题）。

**Q4: 中断后重跑，之前算完的会重算吗？**
不会。断点续传按染色体粒度生效：`01_xpclr/` 里已存在且非空的染色体输出会被跳过（日志显示"跳过已完成染色体"），只补算缺失的部分。

**Q5: 结果只有部分基因组时，Top 表还能用吗？**
能用但要留心：有染色体失败时，Top 50 只是"已算完部分"的前 50，跨全基因组比较的排名可能偏高或偏低，建议先解决失败染色体（看日志 stderr 摘要）再定论。

**Q6: 样本列表路径里带逗号为什么报错？**
xpclr 工具约定：`--samplesA/B` 的值含逗号=内联样本名单，不含逗号=文件路径。路径带逗号会被误当成内联名单解析，故本模块直接在参数校验时拦下。
