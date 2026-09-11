# anno_curate - GSA提交前注释校正 | Pre-GSA Annotation Curation

一句话理解：**把基因注释文件里的结构错误自动修好，再给你一份「提交 GSA 前该人工核对哪些转录本」的清单**。
输入一个 GFF3/GTF 注释文件（可选带上基因组 FASTA 和转录组 BAM 作为证据），输出修复后的规范 GFF3 + 一张逐转录本的问题诊断表；跑一遍、照着清单人工修订、再跑一遍，直到清单干净再提交。

## 功能概述 | Overview

- **单命令一键完成「修复 → 诊断」**：先跑 11 步结构修复，再对修复结果做质量诊断，一步到位
- 修复项：坐标颠倒交换、环状 Parent 解开、gene/mRNA/exon/UTR 补建、重复子特征 ID 去重、悬空 mRNA 隔离、CDS 相位校正、按层级规范排序
- 诊断项：相位缺失/错乱、CDS 过短/过长、UTR 占比过高、CPC2 编码潜能（判 noncoding）、外显子 RNA 零覆盖
- GFF3/GTF 输入自动识别，**输出统一为规范 GFF3**（GTF 输入也会被转换）
- 与 GSAman（TBtools GUI 工具）规则对齐的独立 Python 重实现（行为差异见 FAQ），超算命令行可直接调用
- 两级断点续传：换诊断参数重跑只补诊断、不重跑修复（`--force` 全部重跑）

## 快速开始 | Quick Start

```bash
biopytools anno_curate -i anno.gff3 -o out_dir/
```

最小输入就是一个注释文件。想加上两个可选证据源（编码潜能 + RNA 覆盖）时再加 `-g` 和 `--rna-bam`（见参数说明）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表，后面的参数说明都会用到：

| 术语<br>Term | 通俗理解<br>In plain words |
|------|----------|
| GFF3/GTF | 基因注释文件，记录「每条基因在哪、由哪些零件组成」的清单；GFF3 和 GTF 是同一件事的两种记账格式 |
| GFF3 层级 | 省-市-区三级：gene（省）→ mRNA（市）→ exon/CDS/UTR（区）。哪一级缺了、接错了，文件就是坏的 |
| 转录本<br>transcript | 一条基因的「成品说明书版本」；一条基因可以有多个版本（可变剪接） |
| CDS | 说明书里真正「翻译成蛋白」的正文段落 |
| CDS 相位<br>phase | 「从第几个字开始读句子」：密码子 3 个字一组，段落开头错位 1-2 个字，后面整段全读歪。GFF3 里每个 CDS 段都要标 0/1/2 |
| UTR | 正文前后的「开场白（5'UTR）和收尾（3'UTR）」：不翻译，但开场白/收尾占整条转录本比例反常地大，多半是注释边界画错了 |
| 悬空 mRNA<br>dangling mRNA | 「有户口没房子」：登记了 mRNA 行，名下却没有任何 exon/CDS/UTR 零件，多半是残缺记录 |
| 零覆盖外显子<br>zero-coverage exon | 「没人路过的房间」：注释预测出了这个外显子，但所有转录组数据里一条 reads 都没落到这里——房间可能根本不存在（预测错误） |
| CPC2 | 一个判断「一段序列是不是真的编码蛋白」的工具；把 CDS 序列喂给它，它回答 coding / noncoding |
| GSA | 国家基因组科学数据中心，组装+注释完成后要提交入库的官方仓库；注释格式有硬性要求 |

## 输入 | Input

### 注释文件（必需）

GFF3 或 GTF，按属性列写法自动识别，支持 `.gff` / `.gff3` / `.gtf` 后缀：

- 必须每行 9 列（tab 分隔）；列数不对或坐标不是数字会**收集全部问题行一次性报错退出**
- 起点大于终点的颠倒坐标会在解析时自动交换
- GTF 输入会自动转换（子特征行 `transcript_id` → `Parent`；transcript 行自身 `transcript_id`/`gene_id` → `ID`/`Parent`），输出统一为 GFF3
- 文件末尾的 `##FASTA` 区块会原样保留到输出尾部

注意两类**硬失败**（程序不会替你猜）：

- **重复的 transcript ID**：直接报错退出，列出全部重复 ID 与次数，须先人工改掉再跑
- 列数/坐标非法：报错退出并列出全部问题行号

### 基因组 FASTA（可选，`-g`）

提供后触发 CPC2 编码潜能诊断。要求 FASTA 中的序列名与注释文件的 seqid 一致（不一致的转录本会逐条 WARNING，提取不到任何 CDS 时该项降级跳过）；已有 `.fai` 索引会直接复用，没有则纯 Python 扫描自建，不改动你的 FASTA。

### 转录组 BAM（可选，`--rna-bam`）

提供后触发外显子零覆盖诊断。二代 RNA-seq 与三代 Iso-Seq/ONT 转录组 BAM 可混合传入（可一次多个）；建议**二三代都给**，不同技术互相补盲区。BAM 缺索引会自动在同目录建（`samtools index`），建不了（如无写权限）则该 BAM 跳过并 WARNING。

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** 告诉程序读哪个注释文件、结果写到哪。`-o` 的目录不存在会自动创建；输出文件名里的「样品名」取自输入文件名（剥掉全部后缀层，如 `anno.gff3` → `anno`）。**一般只需要给这两个。**

| 参数<br>Option | 默认值<br>Default | 类型<br>Type | 说明<br>Description |
|------|--------|------|------|
| `-i, --input` | 必填<br>required | str | 输入 GFF3/GTF 注释文件 |
| `-o, --output-dir` | `./anno_curate_output` | str | 输出目录 |

### 诊断阈值 | Diagnosis thresholds

**通俗理解|In plain words:** 这组参数决定「什么样的转录本要被点名」。调严（阈值变小/变短）点名的人变多、误报也变多；调松则漏报增多。**默认值全部来自 GSAman 的经验值，不是 GSA 官方标准，一般不用动**——被点名的条目也未必真错，报警只是「请人工看一眼」。

| 参数<br>Option | 默认值<br>Default | 类型<br>Type | 说明<br>Description |
|------|--------|------|------|
| `--min-cds-nt` | `90` | int | CDS 总长低于此值（nt）标记 `CDS too short` |
| `--max-cds-nt` | `10500` | int | CDS 总长高于此值（nt）标记 `CDS too long` |
| `--utr5-frac-max` | `0.15` | float | 5'UTR 占转录本跨度比例基准，超过有效阈值标记 |
| `--utr3-frac-max` | `0.30` | float | 3'UTR 占转录本跨度比例基准，超过有效阈值标记 |
| `--relax` | `0.5` | float | 松弛系数：UTR 有效阈值 = 基准 × (1 + relax)，即默认放到 0.225 / 0.45；负值自动截为 0 |
| `--no-check-utr` | 关<br>off（默认开检查） | flag | 关闭 UTR 比例检查 |

UTR 检查默认开启（GSAman GUI 默认关）：UTR 很多是修复步骤自己补建的，比例检查正是验证补建是否合理；不想看这类报警就加 `--no-check-utr`。

### RNA 证据 | RNA evidence

**通俗理解|In plain words:** 给程序「目击证人」——转录组 reads。注释说这里有个外显子，如果所有 BAM 里一条 reads 都没有，这个外显子就很可疑。判据只有「零覆盖 / 有覆盖」两档、不做深度阈值（三代 reads 长而稀疏，按深度卡必然误报）。**有转录组数据就建议给，二三代一起给最好；没有就不给，跳过这项。**

| 参数<br>Option | 默认值<br>Default | 类型<br>Type | 说明<br>Description |
|------|--------|------|------|
| `--rna-bam` | 无<br>none | str (可多值) | 转录组 BAM，一次传多个；触发外显子零覆盖诊断 |
| `--samtools-path` | 自动检测 | str | samtools 路径覆盖（默认查 align 域环境，仅 RNA 校验使用） |

### CPC2 编码潜能 | CPC2 coding potential

**通俗理解|In plain words:** 让 CPC2 当「测谎仪」，判断每条 CDS 是不是真像一段会编码蛋白的序列。给了 `-g` 才启用；CPC2 缺失或运行失败只 WARNING 降级跳过，不阻塞其余诊断。**一般不用动路径参数**（默认从 annot 域环境自动定位）。

| 参数<br>Option | 默认值<br>Default | 类型<br>Type | 说明<br>Description |
|------|--------|------|------|
| `-g, --genome` | 无<br>none | str | 基因组 FASTA，触发 CPC2 编码潜能诊断 |
| `--cpc2-path` | 自动检测 | str | CPC2 入口路径覆盖（默认查 annot 域环境，支持 `~`/环境变量/config.yml） |

### 运行控制 | Run control

**通俗理解|In plain words:** 管「跑多少、跑多吵」。`--skip-diagnosis` 用于只想要修好的 GFF3、不要诊断清单的快跑场景；`--force` 用于强制全部重算。**平时都不用给。**

| 参数<br>Option | 默认值<br>Default | 类型<br>Type | 说明<br>Description |
|------|--------|------|------|
| `--skip-diagnosis` | `False` | flag | 只修不诊（与 `-g`/`--rna-bam` 互斥，给了会直接报错） |
| `--force` | `False` | flag | 无视两级断点，修复+诊断全部重跑；上一轮的条件附加文件（悬空/相位 GFF3）与旧诊断目录一并清掉，不残留陈旧副本 |
| `--log-level` | `INFO` | str | 日志级别：DEBUG / INFO / WARNING / ERROR |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先把文件读成内存里的「家谱树」，修好每一层缺胳膊少腿的记录，再用同一棵树做体检，最后一次性写出规范 GFF3 和诊断清单。

```text
输入 GFF3/GTF
    │
    ▼
解析(单遍读入;GTF转GFF3属性;坐标交换;列数/坐标错误收集后报错)
    │
    ▼
11步结构修复(fixer)
    ├─ ② ID/Parent 冲突:环状Parent解开;重复transcript ID→硬失败(人工先解决)
    ├─ ③ gene 行清理:全文件无 mRNA 时删全部 gene
    ├─ ④ mRNA 补建:有零件没妈妈的,按 Parent 聚合补建(source=Recall)
    ├─ ⑤ UTR 补建:exon−CDS 拆分 / mRNA−CDS 差集(按链方向定 5'/3')
    ├─ ⑥ exon 补建:仅无 exon 行的转录本,CDS+UTR 排序合并
    ├─ ⑦ gene 补建:无 gene 的 mRNA 按 Parent 回退链新建
    ├─ ⑧ 子特征 ID 去重:重复 ID 加 .u<N> 后缀
    ├─ ⑨ 悬空 mRNA 隔离 → .fixed.dangling.gff3
    └─ ⑩ CDS 相位校正(总长非3倍数不修只报)
    │
    ▼
排序输出(⑪):染色体首现序 + 起始位置 + 层级(gene→mRNA→exon/UTR/CDS)
    │
    ▼
质量诊断(diagnosis,--skip-diagnosis 可跳过)
    ├─ 结构检查:相位(NO/INVALID)、CDS长度、UTR比例
    ├─ CPC2 编码潜能(给了 -g):提取CDS → CPC2 → noncoding 标记
    └─ RNA 零覆盖(给了 --rna-bam):外显子BED → samtools bedcov → 两级判定
    │
    ▼
输出:01_fixed/ + 02_diagnosis/ + 00_pipeline_info/ + 99_logs/
```

相位算法只信 CDS 长度、不信原 phase：CDS 总长能被 3 整除时一律按公式重算覆盖（`phase[0]=0`，其后 `phase[i] = cumLen%3==0 ? 0 : 3−cumLen%3`）；`+` 链按起点升序、`-` 链降序累加。总长非 3 倍数属于结构本身错了，**不修只报**——自动「改正」反而会制造错误注释。

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml           # 模块版本;启用时附 CPC2 / samtools 版本
├── 01_fixed/
│   ├── {sample}.fixed.gff3             # 主输出:修复后的规范 GFF3
│   ├── {sample}.fixed.dangling.gff3    # 悬空 mRNA 及其父 gene(有内容才生成)
│   ├── {sample}_phase_corrected.gff3   # 被相位校正的转录本全记录+父gene(有内容才生成)
│   └── {sample}_phase_problematic.gff3 # 总长非3倍数的转录本(有内容才生成)
├── 02_diagnosis/                       # --skip-diagnosis 时整目录不生成
│   ├── {sample}.diagnosis.tsv          # 逐转录本问题清单(核心结果)
│   ├── {sample}.diagnosis.summary.txt  # 总数/各类计数/证据源状态
│   └── .diagnosis_meta.json            # 诊断参数指纹(隐藏文件,断点续传用,勿手删)
└── 99_logs/
    ├── anno_curate.log                 # 全量日志
    ├── anno_curate.out.log             # INFO 及以下(stdout)
    └── anno_curate.err.log             # WARNING 及以上(stderr)
```

- `{sample}` = 输入文件名剥掉全部后缀层（`anno.gff3` → `anno`；`anno.gff3.gz` 同样剥干净）
- 附加文件「有内容才生成」：没有悬空 mRNA / 没有被校正 / 没有 problematic 转录本时，对应文件不出现，**空文件不落盘**
- 中间文件（CDS fasta、外显子 BED、bedcov 输出）写在 `output_dir/tmp/`，诊断完成自动清理

**断点续传两级粒度**：

1. **修复级**：`{sample}.fixed.gff3` 已存在 → 跳过整个修复，直接复用
2. **诊断级**：诊断输出存在**且** `.diagnosis_meta.json` 记录的参数指纹与本次一致 → 跳过诊断；指纹不一致（换了 `--rna-bam`/阈值/`-g` 等）→ **只重跑诊断**，修复输出复用

所以「换个 --rna-bam 重新诊断」不用连修复一起重跑；`--force` 无视两级全部重算。

## 结果解读 | Interpreting Results {#interpreting-results}

### 1. 诊断清单（`{sample}.diagnosis.tsv`）

**通俗理解|In plain words:** 这就是「人工核对工作清单」——每个被点名的转录本一行，tags 列写明所有问题。**提交 GSA 前把这张表过一遍：NO RNA COVERAGE 的优先看，其余按量力抽查。**

```text
#Found totally 3 transcripts with potential issues.
transcript_id    chr    strand    cds_len    zero_cov_exons    tags
evm.TU1    chr2    -    84    0/3    CDS too short
evm.TU2    chr3    +    1200    5/5    NO RNA COVERAGE
evm.TU3    chr4    +    990    1/6    ZERO-COV EXON 1/6
```

- 只有**有问题的转录本**才进这张表；`zero_cov_exons` 列在未提供 `--rna-bam` 时留空
- 各标签的读法与优先级：

| 标签<br>Tag | 含义<br>Meaning | 怎么处理<br>Action |
|------|------|------|
| `NO RNA COVERAGE` | 全部外显子在所有 BAM 里零 reads，**最可疑**（很可能是假基因模型） | **最优先人工核对**：IGV 里看该位点有无表达；确认是假的就删 |
| `ZERO-COV EXON n/m` | m 个外显子里有 n 个零覆盖 | 人工核对定位错误外显子/过度预测（内含子保留、嵌合预测常见） |
| `NO PHASE` | 有 CDS 段缺相位（`.`） | 修复步骤本应已补好；修复后仍出现说明结构问题更深（见 FAQ Q2） |
| `INVALID PHASE` | 有 CDS 段相位与推算不符 | 同上；修复后仍出现提示该转录本结构异常，须人工看 |
| `CDS 总长非3倍数|total CDS length not multiple of 3` | CDS 拼起来不是 3 的倍数，密码子读不成句 | **只报不修**（防自动改正制造错误注释）：人工核对 CDS 边界 |
| `CDS too short` / `CDS too long` | CDS 总长越过长度阈值（经验值） | 抽查：小基因家族/长基因可能天然越线 |
| `5'UTR fraction x exceeds y` / `3'UTR fraction x exceeds y` | UTR 占比超过有效阈值（基准×(1+relax)） | 抽查边界画错/漏剪接位点；误报多可调 `--relax` |
| `noncoding` | CPC2 判该 CDS 不像编码序列 | 与 RNA 证据交叉看；ncRNA 注释、短 CDS 易误判 |

> **注意|Warning:** 阈值均为经验值、非 GSA 官方标准，**报警不等于必错**——诊断清单是「请人工看一眼」的排序线索，不是删改依据。修复是自动的，但「要不要信这个基因模型」永远由人决定。

### 2. 摘要（`{sample}.diagnosis.summary.txt`）

**通俗理解|In plain words:** 清单的「抬头」——总共多少转录本、多少有问题、每类问题几条、两个证据源跑没跑成。先看这张再决定花多少力气逐条核对。

```text
total transcripts: 45307
transcripts with issues: 812 (1.8%)
5'UTR fraction 0.26 exceeds 0.22: 18
CDS too short: 156
NO RNA COVERAGE: 12
ZERO-COV: 97
...
CPC2: 完成|done (noncoding 23/45307)
RNA coverage: 完成|done (2 BAMs, 97 transcripts with zero-cov exons)
```

- 相位/长度/NO RNA COVERAGE 等标签按标签名计数；UTR 与「CDS 总长非 3 倍数」类标签带具体数值，**每档数值各计一行**（如 `5'UTR fraction 0.26 exceeds 0.22: 18`），看总量时按前缀归并相加即可
- `CPC2` / `RNA coverage` 行显示启用状态：`完成|done(...)`、`跳过|skipped(...)`（未给输入）或 `失败降级|degraded(...)`（工具失败但不阻塞其余诊断）

### 3. 修复输出（`01_fixed/`）

- `{sample}.fixed.gff3`：主产物，可直接作为 GSA 提交的注释文件；被补建的特征 source 标为 `anno_curate`（补建 mRNA 标 `Recall`），便于 grep 审计
- `{sample}_phase_corrected.gff3`：本次被相位校正的转录本**完整记录**（含全部子特征与父 gene），供人工抽查校正是否合理
- `{sample}_phase_problematic.gff3`：总长非 3 倍数、相位救不了的转录本，**必须人工处理**
- `{sample}.fixed.dangling.gff3`：被隔离的悬空 mRNA（连同独占的父 gene），不进主输出也不进诊断，决定保留就人工补全零件后重跑

## 参数选择建议 | Parameter Guidance

- **阈值组（`--min-cds-nt` 等）：默认值即 GSAman 经验值，一般不动**。真要调的场景：注释对象是极小基因组/大基因家族（如 NLR 抗病基因 CDS 普遍偏长）时，可按物种实际情况放宽 `--max-cds-nt`
- **`--relax`**：UTR 报警刷屏、抽查后多数是真实长 UTR 时调大（如 1.0 → 阈值放到基准×2）；确认阈值没吃掉真问题前不要盲目调
- **`--rna-bam`**：有转录组数据就给，**建议二代+三代一起给**——聚合口径是「任一 BAM 有 reads 即算覆盖」，多给一个 BAM 只会减少误报（前提是样本/组织匹配）
- **`--skip-diagnosis`**：只需要修好的 GFF3（比如上轮已核对过清单）时用；注意它与 `-g`/`--rna-bam` 互斥
- **`--force`**：换输入文件内容但沿用了旧输出目录、或怀疑断点状态不对时用；日常重跑不需要（换参数会自动按指纹判断）
- 无 `--threads` 参数：纯 Python 单遍文本变换 + CPC2 单线程，没有可并行处

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入GFF3/GTF｜Input GFF3/GTF |
| `--output-dir, -o` | `./anno_curate_output` | Path | 输出目录｜Output directory |
| `--genome, -g` | — | Path | 基因组FASTA(触发CPC2)｜Genome FASTA (enables CPC2) |
| `--rna-bam` | — |  | 转录组BAM,可多次传(触发零覆盖校验)｜RNA BAM, repeatable (enables zero-coverage check) |
| `--skip-diagnosis` | `False` |  | 只修不诊｜Skip diagnosis |
| `--check-utr/--no-check-utr` | `True` |  | UTR比例检查｜UTR fraction check |
| `--min-cds-nt` | `90` | int |  |
| `--max-cds-nt` | `10500` | int |  |
| `--utr5-frac-max` | `0.15` | float |  |
| `--utr3-frac-max` | `0.3` | float |  |
| `--relax` | `0.5` | float | 松弛系数,有效阈值=基准×(1+relax)｜Relax factor |
| `--cpc2-path` | — |  | CPC2路径｜CPC2 path |
| `--samtools-path` | — |  | samtools路径｜samtools |
| `--force` | `False` |  | 忽略断点重跑｜Rerun all |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR |  |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GFF3/GTF｜Input GFF3/GTF |
| `-o, --output-dir` | `./anno_curate_output` |  | 输出目录｜Output directory |
| `-g, --genome` | — |  | 基因组FASTA(触发CPC2编码潜能校验)｜Genome FASTA (enables CPC2 check) |
| `--rna-bam` | — |  | 转录组BAM(可多个,触发外显子零覆盖校验)｜RNA BAMs (multiple allowed, enables zero-coverage check) |
| `--skip-diagnosis` | `False` | store_true | 只修不诊｜Skip diagnosis |
| `--no-check-utr` | `True` | store_false | 关闭UTR比例检查｜Disable UTR fraction check |
| `--min-cds-nt` | `90` | int | CDS过短阈值(nt)｜CDS too-short threshold |
| `--max-cds-nt` | `10500` | int | CDS过长阈值(nt)｜CDS too-long threshold |
| `--utr5-frac-max` | `0.15` | float | 5'UTR比例基准｜5'UTR fraction baseline |
| `--utr3-frac-max` | `0.3` | float | 3'UTR比例基准｜3'UTR fraction baseline |
| `--relax` | `0.5` | float | 松弛系数(有效阈值=基准×(1+relax))｜Relax factor |
| `--cpc2-path` | `` |  | CPC2路径(默认查annot域)｜CPC2 path (defaults to annot domain env) |
| `--samtools-path` | `` |  | samtools路径(默认查align域)｜samtools path |
| `--force` | `False` | store_true | 忽略断点全部重跑｜Rerun everything |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3 标准库（**零第三方 Python 依赖**，解析/修复/诊断全内置）
- samtools：`align` 域环境（默认 `~/miniforge3/envs/align/bin/samtools`），仅 `--rna-bam` 时使用（bedcov / index）
- CPC2：`annot` 域环境（默认 `~/miniforge3/envs/annot/bin/CPC2.py`），仅 `-g` 时使用；可用 `--cpc2-path` 或环境变量 `CPC2_PATH` / `~/.config/biopytools/config.yml` 覆盖
- 无其他外部工具（排序/解析全部内存完成）

## 常见问题 | FAQ

**Q1：报「重复转录本ID|duplicate transcript ID」直接退出了，怎么办？**
这是故意的硬失败：同一 ID 对应多条不同记录时，程序无法替你决定哪条是对的。错误信息会列出**全部**重复 ID 和出现次数（GSAman 只报第一个），按清单人工把 ID 改唯一后重跑即可。

**Q2：修复后诊断里仍有 INVALID PHASE / CDS 总长非3倍数？**
设计如此：CDS 总长非 3 倍数说明结构本身错了（缺/多一段 CDS），自动「改正」只会制造错误注释，所以**只报不修**（进 `_phase_problematic.gff3` 与诊断清单）。到 IGV 里人工核对 CDS 边界，改完再跑一遍本模块验证。INVALID PHASE 在正常流程中已被相位校正覆盖（进 `_phase_corrected.gff3`），修复后仍出现即说明该转录本结构异常，同样须人工处理。

**Q3：输入是 GTF，输出是什么格式？**
输出统一为规范 GFF3（`##gff-version 3` 头 + 排序后的记录）；GTF 的 `gene_id`/`transcript_id` 在解析时即转为 GFF3 的 `Parent`/`ID` 语义。

**Q4：和 GSAman（TBtools）是什么关系？**
独立的 Python 重实现：算法逻辑/阈值/触发条件逐条对齐 GSAman 的注释校正功能（GXF Fix / Quick Diagnosis / CDS Phase Corrector 三件套），**未复制任何其源码**（TBtools 许可禁止逆向再分发）。为超算命令行场景重写，并做了少量行为改进：

| # | 差异<br>Difference | 理由<br>Reason |
|---|------|------|
| 1 | 染色体顺序 = 输入首次出现顺序（GSAman 按外部排序） | 不打乱用户 scaffold 原序 |
| 2 | 补建特征 source=`anno_curate`（mRNA 补建用 `Recall`） | 不冒名 TBtools，便于审计 |
| 3 | UTR 比例检查默认开（GSAman GUI 默认关） | CLI 无勾选摩擦；UTR 是补建产物需验证 |
| 4 | 重复 transcript ID 报全部明细（GSAman 只抛第一个） | 一次改完少来回 |
| 5 | CDS 匹配统一大小写敏感 | 消除原版内部不一致 |
| 6 | 列数/坐标非法解析即硬失败 | 报错更早更明确 |
| 7 | 诊断 TSV 增加 chr/strand/cds_len/zero_cov_exons 列 | 方便筛选排序 |
| 8 | **新增 RNA 零覆盖校验**（`--rna-bam`） | GSAman 没有；人工核对最有力的证据源 |
| 9 | tRNA/rRNA 等非 mRNA 转录本显式按层级排序 | 修复原版排序隐患 |
| 10 | feature 类型大小写变体（如小写 `cds`）WARNING | 防静默失效（见 Q5） |

**Q5：日志里 WARNING「feature类型大小写变体|case-variant feature type」是什么？**
输入里存在小写 `cds`、`mrna` 这类非规范大小写的行。主流程的类型匹配是大小写敏感的，这些行**不会参与相位/UTR/exon 逻辑**；为防你毫无感知，程序会按类型汇总出现次数与行号样本打 WARNING，但**不改写输入**（不越修复授权）。请人工把类型名统一为规范大小写后重跑。

**Q6：换了个 `--rna-bam` 重跑，为什么修复没有重跑？**
两级断点：修复输出存在即跳过修复；诊断按 `.diagnosis_meta.json` 参数指纹判断，参数变了就**只重跑诊断**（复用修复结果）。想全部重来加 `--force`。

**Q7：`--skip-diagnosis` 和 `-g`/`--rna-bam` 一起给了怎么报错？**
两者互斥：`--skip-diagnosis` 意为「只要修好的 GFF3」，而 `-g`/`--rna-bam` 只影响诊断，给了没有意义。去掉其一即可。

**Q8：悬空 mRNA 去哪了？**
连同其独占的父 gene 一起隔离进 `01_fixed/{sample}.fixed.dangling.gff3`，不进主输出也不进诊断清单。想保留这些模型就人工补全 exon/CDS 后重新提交；不想要就无视。

---

> **参数表说明|Parameter tables note:** 本文参数说明各表由 `biopytools/anno_curate/main.py` 的 argparse 定义自动提取生成（生成器一次性产物不入仓库），**禁止手写、参数有变须重新提取**，保证文档与代码不漂移。
