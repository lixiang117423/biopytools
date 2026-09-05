# 组装间结构变异检测 | Assembly-to-Assembly SV Calling

一句话理解：**把一个参考基因组和其它每个基因组两两比对，自动找出它们之间的大片段结构差异**（缺失、插入、倒位、重复、易位），并合并成一张群体层面的结构变异清单，解决「不同个体/品系的基因组之间到底差在哪」的问题。

## 功能概述 | Overview { #overview }

- 输入一个样本清单(fof)，指定一个参考样本，其余样本作为查询逐一分析
- 每个查询样本：minimap2 比对到参考(asm5/asm10/asm20 预设) → svim-asm 检测结构变异(SV) → SURVIVOR 把各样本的 SV 合并
- 合并后的群体 SV 清单输出为 `pan_sv_survivor.vcf`，并生成每样本 SV 类型统计表 `sv_summary.tsv`
- 同时输出每条 SV 的代表序列 `pan_sv_sequences.fa` 和样本×SV 的 PAV 矩阵(`pav_matrix.tsv` / `pav_binary.tsv`)
- 输出每条 SV 上下游侧翼序列(`06_sv_flanks/`，默认 ±300bp 可调) + 每样本基因型表(TSV/FASTA/XLSX)，可直接用于引物设计
- 全程断点续传：比对、SV 调用、SURVIVOR 合并已完成的步骤重跑时自动跳过
- 六个依赖工具(minimap2 / samtools / svim-asm / bcftools / bedtools / SURVIVOR)统一在 `align` 域环境中调用(比对管道要求单环境)

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools genome2sv -i samples.fof -r ref -o results/
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| 结构变异(SV) | 基因组上「一大段」序列的变化，像书里整页被撕掉(缺失)、多塞进一页(插入)、一页颠倒(倒位)、一页被复印成两份(重复)、两页互换位置(易位) |
| fof 清单 | 一个「点名表」文件，每行是 `样本名 TAB 基因组文件路径`，告诉程序去分析谁 |
| 参考/查询样本 | 把谁当「标准答案」(参考)，谁拿去和它比(查询) |
| asm5/asm10/asm20 | minimap2 比对时对「两个基因组有多像」的预期档位：asm5=很接近(约≤5%差异)，asm10/asm20=差异更大时用 |
| svim-asm | 专门从「两个组装好的基因组比对」里找 SV 的工具 |
| SURVIVOR 合并 | 多个样本各自报出的 SV，把「位置接近、类型一致」的合并成同一个，避免重复计数 |
| 断点(breakpoint) | 一段 SV 的两端坐标，像「被撕掉那页」的开头和结尾页码 |
| 侧翼序列(flank) | 以 SV 断点为中心向两边各多取一段(默认 300bp)的参考序列，像「被撕那页前后各多复印 300 个字」，让 PCR 引物能落在变异之外的正常区域上 |
| GT/LN/QV | 每个样本对这条 SV 的三个值：GT=基因型(如 1/1 表示两套拷贝都有此变异，./.=没测出)，LN=变异长度，QV=质量值(越高越可信，`.`=未提供) |

## 输入 | Input { #input }

### fof 样本清单

每行 `样本名<TAB>基因组FASTA路径`；`#` 开头的行为注释，空行忽略；非注释行必须含 TAB 分隔符。

```text
# 样本名<TAB>基因组路径
ref     /data/genomes/ref.fa
sampleA /data/genomes/sampleA.fa
sampleB /data/genomes/sampleB.fa
```

### 基因组 FASTA

每个样本(含参考)一个 FASTA 文件；`-r` 指定的参考样本名必须与 fof 第一列中的某个名字完全一致，且 fof 中除参考外至少还要有一个查询样本。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 三个必填项分别是「点名表在哪、谁当参考、结果放哪」。`-r` 指定的参考样本名必须能在 fof 里找到，否则直接报错。

### 比对参数 | Alignment

**通俗理解|In plain words:** `--preset` 决定比对时按「两个基因组有多像」选档，值越小越严格。默认 `asm5` 适合同物种内差异很小的组装；跨亚种/近缘种可试 `asm10`，差异更大用 `asm20`。`-t` 是每个样本比对时用的线程数，调大跑得快但更吃 CPU，一般用默认 12 即可。

### SV 调用参数 | SV calling

**通俗理解|In plain words:** `--svim-mode` 决定按「单套基因组」(haploid) 还是「双套」(diploid) 来报 SV。普通参考基因组比对用默认 `haploid` 就行；只有当你明确在比对双倍体/两套单倍型时才改 `diploid`。一般不用动。

### 合并参数 | SURVIVOR merge

**通俗理解|In plain words:** 这一组决定「什么样的 SV 才算同一个、哪些值得保留」。`--max-dist`(默认 1000)是允许两个断点最多差多少 bp 仍视为同一个 SV，调大会更宽松地合并；`--min-sv-length`(默认 50)只保留足够长的 SV；`--min-support`(默认 1)是「至少几个样本都报出这个 SV 才保留」，调大能滤掉只在单个样本里出现的噪音，但也可能丢掉真实稀有变异；`--survivor-type` 和 `--survivor-strand`(默认都为 1)控制合并时是否要求类型/链方向一致。**绝大多数项目用默认值即可，一般不用动。**

### 侧翼参数 | Flank

**通俗理解|In plain words:** `--flank`(默认 300)决定步骤 6 给每条 SV 的断点向两边各多取多少 bp 参考序列。管什么：引物设计时让引物落在 SV 之外的正常区域，需要「变异位点 + 两侧足够的落脚点」。调大：序列更长，离变异更远，适合设计更长片段(如长距离 PCR、探针)；调小：序列更紧凑，适合小片段扩增。300bp 是常规 PCR 引物设计的常用跨度，**一般不用动**；要 SV 本体序列(不含侧翼)可设为 0。

### 运行参数 | Runtime

**通俗理解|In plain words:** `--log-level` 控制日志详细程度，默认 INFO，排查问题时才调到 DEBUG，一般不用动。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先把参考准备好(建索引)，再逐个样本「比对→找SV」，最后把所有样本的 SV 汇总合并、统计。

```text
输入 fof 清单 + 参考样本名
    |
    ▼
步骤0: 参考预处理(软链到 reference/ + samtools faidx 建索引)
    |
    ▼
步骤1: 读取查询样本(fof 中除参考外)
    |
    ▼
步骤2: 逐个样本比对(minimap2 -ax asm5 | samtools sort → _sorted.bam + 索引)
    |
    ▼
步骤3: 逐个样本 SV 调用(svim-asm → 02_svim/<样本>/)
    |
    ▼
步骤4: SURVIVOR 合并 + bcftools stats + SV 类型统计表
    |
    ▼
步骤5: SV 序列提取 + PAV 矩阵(从合并 VCF 生成)
    |
    ▼
步骤6: SV ±flank 侧翼序列 + 每样本 GT/LN/QV 基因型表(TSV/FASTA/XLSX)
    |
    ▼
输出 pan_sv_survivor.vcf + pan_sv_sequences.fa + pav_matrix.tsv + sv_flank300bp.* + sv_summary.tsv
```

## 输出 | Output { #output }

```text
results/
├── reference/                  # 参考基因组软链 + .fai 索引
├── 01_alignment/               # 比对结果
│   ├── sampleA_sorted.bam
│   └── sampleA_sorted.bam.bai
├── 02_svim/                    # 每个查询样本的 SV 调用输出
│   └── sampleA/                # svim-asm 产出的 VCF
├── 03_merged/
│   ├── survivor_input.txt      # SURVIVOR 输入清单(各样本 VCF 绝对路径)
│   └── pan_sv_survivor.vcf     # 合并后的群体 SV 清单(核心结果)
├── 04_stats/
│   ├── sv_summary.tsv          # 每样本(含 merged)的 SV 类型统计表
│   ├── pan_sv_survivor.stats   # bcftools stats 输出
│   ├── pav_matrix.tsv          # PAV 矩阵(SV 元数据 + 每样本 0/1)
│   └── pav_binary.tsv          # 纯 0/1 PAV 矩阵(R 可直接 as.matrix)
├── 05_sv_sequences/
│   └── pan_sv_sequences.fa     # 每条 SV 一条代表序列(FASTA)
├── 06_sv_flanks/
│   ├── sv_flank300bp.fa        # 每条 SV ±300bp 侧翼序列(FASTA)
│   ├── sv_flank300bp.tsv       # 侧翼坐标 + 每样本 GT/LN/QV + 序列(程序友好)
│   └── sv_flank300bp.xlsx      # 同 TSV,Excel 查看(超长序列只入 fasta)
├── 00_pipeline_info/
│   └── software_versions.yml   # 软件版本与运行参数
└── 99_logs/
    └── genome2sv.log           # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. pan_sv_survivor.vcf(核心结果)

**通俗理解|In plain words:** 这是「群体里都有哪些结构变异」的最终清单，一行一个 SV，可拿去做下游分析或画图。

- 每行是一个合并后的 SV；`INFO` 字段里的 `SUPP` 表示有多少个样本支持这个 SV(受 `--min-support` 过滤)
- `SVTYPE` 标注类型：`DEL`(缺失)、`INS`(插入)、`INV`(倒位)、`DUP`(重复)、`BND`(易位)

### 2. sv_summary.tsv(统计表)

**通俗理解|In plain words:** 一张「谁有多少种 SV」的汇总表，一行一个样本，最后一行是合并结果。

```text
sample    total    INS    DEL    INV    DUP    BND    OTHER
sampleA   120      30     60     10     8      10     2
merged    150      35     80     15     10     8      2
```

- `total` 是该样本(或合并后)检测到的 SV 总数，后面几列是各类型计数，`OTHER` 是其它类型
- 合并后的 `merged` 行 SV 数一般不少于单个样本，因为汇总了所有样本的 SV

### 3. pan_sv_sequences.fa(SV 序列)

**通俗理解|In plain words:** 把每条 SV 对应的那段 DNA 序列抠出来存成 FASTA，可直接拿去 blast、设计引物或做注释。

- 序列来源按类型：INS 是插入序列本身(来自 ALT 字段)、DEL 是被删除的参考序列(含断点首碱基)、INV 是倒位区序列(反向互补)、DUP 是重复单元参考序列
- `BND`(易位)没有明确的单一区间，不输出序列(日志会计数)
- 序列名形如 `pan_sv.INS.00001`,header 里带坐标/长度/来源/支持样本;`sv_id` 与 PAV 矩阵共用，可交叉对照
- 若与 PAV 矩阵行数对不上，差额就是 BND 或无法提取的记录(日志有"跳过|skipped"计数)

### 4. pav_matrix.tsv / pav_binary.tsv(PAV 矩阵)

**通俗理解|In plain words:** 一张「哪个样本有哪条 SV」的 0/1 表——1 表示有(present)、0 表示没有(absent)，是单倍型/品种间 SV 差异分析的核心输入。

- `pav_matrix.tsv` 带元数据列(sv_id/chrom/pos/end/svtype/svlen),便于按染色体/类型筛选
- `pav_binary.tsv` 只有 0/1,方便 R 里 `read.table` 后直接 `as.matrix` 做聚类/PCA
- 判定规则：样本基因型含 allele 1(如 `1/1`、`0/1`)记 1,缺失(`./.`)记 0
- 行数与合并 VCF 的 SV 总数一致(含 BND);`sv_id` 与序列 FASTA 对应

### 5. 06_sv_flanks/(SV 侧翼序列与基因型表)

**通俗理解|In plain words:** 给每条 SV 拍一张「连周边环境一起」的照片——以断点为中心，两边各多取 300bp 参考序列，并附上每个样本在这条 SV 上的基因型(_GT)、长度(_LN)和质量(_QV)。拿到手就能直接挑「亲本间有差异」的 SV 去设计 PCR 引物：引物落在两侧翼上，扩增产物长度随 SV 有无而不同，从而把不同基因型区分开。

- 序列区间 = 参考 `[min(pos,end)-flank, max(pos,end)+flank]`(1-based 闭区间)；INS 的 pos=end，插入位点正好在序列中部；靠近染色体首尾时区间会截断(日志给"边界截断"计数)
- `TRA/BND`(易位)没有单一参考区间，不输出(日志计数)；不做额外长度过滤，保留哪些 SV 由 `--min-sv-length` 在合并时决定
- `sv_id`(如 `pan_sv.DEL.00001`)与 `pav_matrix.tsv`、`pan_sv_sequences.fa` 完全一致，可用 vlookup 互相对照；FASTA header 还带 `flank=起点-终点`
- 三个文件内容一致、各有所长：`.fa` 给引物设计软件，`.tsv` 给脚本/R 批量处理，`.xlsx` 给人筛选浏览
- Excel 单元格最多存 32767 字符，超长的序列在 `.xlsx` 里只写占位提示，全长以 `.fa` 为准

### 6. 好坏判据

- 合并后 SV 数量合理、各样本统计表都能正常生成，说明流程正常
- 若某样本 SV 数异常为 0 或极少，检查该样本比对是否失败(看日志)，或它与参考差异是否过大(考虑换 `--preset`)

## 参数选择建议 | Parameter Guidance { #guidance }

- **同物种内差异小**：默认 `--preset asm5` 即可，最常用
- **近缘种/差异较大**：`--preset asm10` 或 `asm20`，比对更宽松
- **只要高可信的群体 SV**：调大 `--min-support`(如 2)，滤掉单样本噪音
- **想保留稀有变异**：`--min-support 1`(默认)保留所有样本各自报出的 SV
- **双倍体组装**：`--svim-mode diploid`
- **要设计 PCR/INDEL 标记引物**：默认 `--flank 300` 即可；要更长同源臂(如长片段 PCR)调大到 500–1000；只要 SV 本体序列设 `--flank 0`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 样本清单 fof(name<TAB>path)｜Sample manifest fof |
| `-r, --ref` | 必填 |  | 参考样本名(fof 第一列)｜Reference sample name |
| `-o, --output-dir` | `./genome2sv_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数(每样本,默认12)｜Threads per sample (default 12) |
| `--preset` | `asm5` | asm5/asm10/asm20 | minimap2 预设(默认asm5)｜minimap2 preset (default asm5) |
| `--svim-mode` | `haploid` | haploid/diploid | svim-asm 模式(默认haploid)｜svim-asm mode (default haploid) |
| `--max-dist` | `1000` | int | SURVIVOR 断点最大距离bp(默认1000)｜SURVIVOR max breakpoint dist (default 1000) |
| `--min-sv-length` | `50` | int | SURVIVOR 最小SV长度bp(默认50)｜SURVIVOR min SV length (default 50) |
| `--survivor-type` | `1` | 0/1 | SV类型一致1/任意0(默认1)｜Require same type (default 1) |
| `--survivor-strand` | `1` | 0/1 | 链方向一致1/任意0(默认1)｜Require same strand (default 1) |
| `--min-support` | `1` | int | SURVIVOR 最小支持调用数(默认1)｜SURVIVOR min supporting callers (default 1) |
| `--flank` | `300` | int | SV 上下游侧翼长度bp(默认300)｜SV flank length bp (default 300) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 样本清单 fof(name<TAB>path)｜Sample manifest fof |
| `-r, --ref` | 必填 |  | 参考样本名(fof 第一列)｜Reference sample name |
| `-o, --output-dir` | `./genome2sv_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数(每样本)｜Threads per sample (default 12) |
| `--preset` | `asm5` | asm5/asm10/asm20 | minimap2 预设｜minimap2 preset (default asm5) |
| `--svim-mode` | `haploid` | haploid/diploid | svim-asm 模式｜svim-asm mode (default haploid) |
| `--max-dist` | `1000` | int | SURVIVOR 断点最大距离(bp)｜SURVIVOR max breakpoint dist (default 1000) |
| `--min-sv-length` | `50` | int | SURVIVOR 最小 SV 长度(bp)｜SURVIVOR min SV length (default 50) |
| `--survivor-type` | `1` | 0/1 | SV 类型一致(1)/任意(0)｜Require same type (default 1) |
| `--survivor-strand` | `1` | 0/1 | 链方向一致(1)/任意(0)｜Require same strand (default 1) |
| `--min-support` | `1` | int | SURVIVOR 最小支持调用数｜SURVIVOR min supporting callers (default 1) |
| `--flank` | `300` | int | SV 上下游侧翼长度 bp(默认300)｜SV flank length bp (default 300) |
| `--log-level` | `INFO` |  | 日志级别｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- minimap2(比对)
- samtools(排序/索引/faidx)
- svim-asm(SV 检测)
- bcftools(合并后统计)
- bedtools(版本探测)
- SURVIVOR(SV 合并)

以上工具统一在 `align` 域环境中调用(无需手动指定，也可用 `MINIMAP2_PATH` 等环境变量覆盖)。

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持。比对步骤按 `_sorted.bam` 和 `.bam.bai` 是否都存在判断跳过；SV 调用按 `02_svim/<样本>/` 下是否已有 VCF 判断跳过。换参数重跑前需先删除对应的旧产物。

**Q2：SURVIVOR 显示「缺失工具」但明明装了？**
SURVIVOR 是子命令式命令行(没有 `--version`)，依赖检查按「能否执行」而非「退出码为 0」判定，属正常。确认 SURVIVOR 在 `align` 环境里即可。

**Q3：fof 报「格式错误」？**
fof 每行必须是 `样本名<TAB>路径`(TAB 分隔)，非注释行缺 TAB 会报错。路径里带 `~` 或环境变量会自动展开。

**Q4：minimap2 和 samtools 不在同一个 conda 环境怎么办？**
默认路径已把两者固定在 `align` 环境，不会出现不一致。若用 `MINIMAP2_PATH`/`SAMTOOLS_PATH` 环境变量覆盖后两者分属不同环境，比对管道会以 samtools 所在环境为准并在日志给出 WARNING——此时请确保该环境里同时装有 minimap2，否则会找不到命令。

**Q5：换 `--flank` 重跑为什么没有重新生成侧翼文件？**
侧翼输出按「三个文件是否都已存在」判断断点续传。改了 `--flank` 后文件名(如 `sv_flank300bp.fa` → `sv_flank500bp.fa`)会变，新生成的文件名不同，不影响；但重跑同 flank 时旧文件不会被覆盖，需要换 flank 长度或先删除 `06_sv_flanks/` 旧产物。

**Q6：为什么部分 SV 没有 xlsx 序列/没有进 FASTA？**
易位(`TRA/BND`)没有单一的参考坐标区间，侧翼序列无从取起，日志有跳过计数；这是 SV 类型本身的特点，不是错误。
