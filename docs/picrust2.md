# PICRUSt2 - 微生物群落功能丰度预测 | PICRUSt2 Functional Abundance Prediction

一句话理解：**只知道「样本里有哪些细菌、各有多少」还不够，PICRUSt2 据此推断「这些细菌整体上具备哪些酶和代谢通路」**，把 16S 的物种丰度表升级成功能丰度表。

## 功能概述 | Overview { #overview }

- 输入代表序列 FASTA + 特征表（自动识别 BIOM/TSV/Excel/Mothur）
- 序列放置（epa-ng/sepp）挂到参考树，隐状态预测（HSP）推断基因家族丰度
- 默认预测 EC 与 KO 两类功能，可选 GO/PFAM/BIGG/CAZY
- 默认做 MetaCyc 代谢通路推断（MinPath + gap filling），可 --no-pathways 跳过
- 预处理贴心：Excel 自动转 TSV、BIOM 行列方向自动纠正、FASTA 头自动清理
- 断点续传：检测到核心产物即跳过整个流程

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools picrust2 -s study_seqs.fna -i seqabun.biom -o picrust2_out
```

最小输入：一个代表序列 FASTA（`-s`）+ 一个特征丰度表（`-i`）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 16S rRNA | 细菌的「身份证」标记基因，测它就知道样本里有哪些菌 |
| 功能预测 | 由「你是谁」推断「你会干什么」（有哪些酶/通路） |
| 参考树 | 已知物种的亲缘关系「家谱」，是预测的骨架 |
| 序列放置(placement) | 把你的新序列挂到参考树的正确位置 |
| NSTI | 放置的「离家距离」，越小说明在参考树里找得到近亲、预测越可信 |
| HSP(隐状态预测) | 借用树上近亲的已知功能，推断新序列的功能 |
| EC / KO | 酶(EC)与基因(KO)的功能分类编号 |
| MetaCyc 通路 | 代谢通路数据库，PICRUSt2 用它推「整条代谢路径」 |
| 分层表(stratified) | 细到「哪个物种贡献了多少功能」的表 |

## 输入 | Input { #input }

- **代表序列 FASTA**（`-s`）：每条 ASV/OTU 一条序列，header 的 ID 须与特征表的观测 ID 一致（不一致会自动清理）。
- **特征表**（`-i`）：行=序列(ASV/OTU)，列=样本，值为计数。支持 BIOM、TSV、Excel（.xlsx/.xls）、Mothur shared。

TSV 特征表示例：

```text
OTU_ID    S1    S2    S3
OTU1      120   45    0
OTU2      8     300   12
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** `-s` 代表序列、`-i` 特征表、`-o` 输出目录（默认 ./picrust2_output）、`-t` 线程数（默认 12）。

### 功能数据库 | Functional databases

**通俗理解|In plain words:** `--in-traits` 决定预测哪些功能库（默认 EC,KO；可加 GO/PFAM/BIGG/CAZY）。**默认 EC+KO 覆盖绝大多数分析需求**；选得越多耗时越长。

### 放置与预测 | Placement & prediction

**通俗理解|In plain words:** `--placement-tool` 默认 epa-ng（快）；`--hsp-method` 默认 mp（多数投票）；`--edge-exponent` 是 HSP 的边权重指数（默认 0.5）。`--max-nsti` 是可信度上限（默认 2.0），NSTI 超过它的序列在宏基因组预测时会被排除。**这些是算法细节，一般不用动。**

### 过滤 | Filtering

**通俗理解|In plain words:** `--min-align`（默认 0.8）、`--min-reads`（默认 1）、`--min-samples`（默认 1）决定哪些序列参与预测。**默认最宽松，一般不用动**；数据噪声大时可适当提高 min-reads/min-samples。

### 通路与流程控制 | Pathway & pipeline control

**通俗理解|In plain words:** 默认会做通路推断；`--no-pathways` 跳过（只要 EC/KO 基因丰度时更快）；`--coverage` 额外算通路覆盖度；`--skip-minpath` 与 `--no-gap-fill` 是通路推断的微调；`--skip-norm` 跳过标记基因拷贝数归一化。`--pipeline` 默认 auto 自动选 split（双域）或 single（单参考）。**日常用默认。**

## 分析流程 | Pipeline { #pipeline }

```text
输入 FASTA + 特征表
    │
    ▼
断点续传检查(path_abun/pred_metagenome 已存在则跳过)
    │
    ▼
预处理: Excel→TSV / BIOM方向自动纠正 / FASTA头清理
    │
    ▼
运行 PICRUSt2 pipeline(放置→HSP→宏基因组预测→通路推断)
    │
    ▼
重组输出到编号目录 + 给 EC/KO/通路表加描述列与均值列
    │
    ▼
生成 software_versions.yml + 清理 work/
```

## 输出 | Output { #output }

```text
picrust2_out/
├── 00_pipeline_info/
│   └── software_versions.yml
├── 01_placement/                     # 序列放置结果(.tre、place_seqs_*)
├── 02_hsp/                           # 标记基因预测与 NSTI(marker_predicted/nsti)
├── 03_metagenome/
│   ├── EC_pred_metagenome_unstrat.tsv    # EC 丰度(加描述列)
│   ├── KO_pred_metagenome_unstrat.tsv    # KO 丰度(加描述列)
│   ├── pred_metagenome_unstrat.tsv.gz    # 原始宏基因组预测表
│   ├── seqtab_norm.tsv.gz                # 归一化序列表
│   └── weighted_nsti.tsv.gz              # 每样本加权 NSTI
├── 04_pathway/
│   └── path_abun_unstrat.tsv             # MetaCyc 通路丰度(加描述列)
└── 99_logs/
    └── picrust2_pipeline.log
```

加 `--stratified` 会额外产出分层表（功能×物种×样本）；加 `--coverage` 产出通路覆盖度表。原始未注释表保留为 `*_raw`。

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. 03_metagenome/ 功能丰度表

**通俗理解|In plain words:** 行是 EC/KO 功能，列是样本，数值是预测丰度，额外带 `description`（功能描述）和 `Mean`（跨样本均值，已按均值降序）。直接拿去做差异分析、绘图。

### 2. 04_pathway/path_abun_unstrat.tsv

**通俗理解|In plain words:** 每个 MetaCyc 代谢通路在各样本的预测丰度。想看「整条代谢路径」而非单个酶时用它。

### 3. weighted_nsti.tsv.gz（可信度指标）

**通俗理解|In plain words:** 每个样本的加权 NSTI——预测可信度的「温度计」。**越低越好**；>2 的样本说明其序列在参考树里找不到近亲，功能预测要谨慎对待。

### 4. 好坏判据

- **weighted_nsti 普遍偏高（>2）**：样本含大量参考库未覆盖的新类群，预测可信度低。
- **分层表与不分层表对不上**：先核对 FASTA ID 与特征表观测 ID 是否一致（本工具已自动清理/转置，通常已解决）。
- **通路缺失**：低丰度或 NSTI 高的序列会被排除（max_nsti=2.0），属正常过滤。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **常规流程**：全部默认（EC,KO + 通路 + epa-ng）。
- **只要基因丰度、不要通路**：加 `--no-pathways`，显著提速。
- **要更细的功能归属**：加 `--stratified`。
- **数据量很大**：`--placement-tool epa-ng`（默认）+ 提高 `-t`。
- **样本 NSTI 普遍高**：可调高 `--max-nsti`（如 3.0）保留更多序列，但要意识到预测可信度下降。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --study-fasta` | 必填 |  | 代表序列FASTA文件｜FASTA of unaligned study sequences |
| `-i, --input` | 必填 |  | 特征表文件(自动识别BIOM/TSV/Excel/Mothur)｜Input table of sequence abundances |
| `-o, --output-dir` | `./picrust2_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--max-nsti` | `2.0` | float | 最大NSTI阈值｜Maximum NSTI value |
| `--stratified` | — |  | 生成分层输出表｜Generate stratified output tables |
| `--in-traits` | `EC,KO` |  | 功能数据库(EC,KO,GO,PFAM,BIGG,CAZY)｜Gene families to predict |
| `--placement-tool` | `epa-ng` | epa-ng/sepp | 序列放置工具｜Placement tool |
| `--hsp-method` | `mp` | mp/emp_prob/pic/scp/subtree_average | 隐状态预测方法｜HSP method |
| `--edge-exponent` | `0.5` | float | HSP edge exponent |
| `--pipeline` | `auto` | auto/split/single | 流程类型: auto/split/single｜Pipeline type |
| `--min-align` | `0.8` | float | 最小比对比例｜Minimum alignment proportion |
| `--min-reads` | `1` | int | 每ASV最小reads数｜Minimum reads per ASV |
| `--min-samples` | `1` | int | 每ASV最小样本数｜Minimum samples per ASV |
| `--no-pathways` | — |  | 跳过通路推断｜Skip pathway inference |
| `--coverage` | — |  | 计算通路覆盖度｜Calculate pathway coverages |
| `--skip-minpath` | — |  | 跳过MinPath｜Skip MinPath |
| `--no-gap-fill` | — |  | 跳过gap filling｜Skip gap filling |
| `--per-sequence-contrib` | — |  | 逐序列运行MinPath｜Run MinPath per sequence |
| `--skip-norm` | — |  | 跳过归一化｜Skip normalization |
| `--remove-intermediate` | — |  | 移除中间文件｜Remove intermediate files |
| `--verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --study-fasta` | 必填 |  | 代表序列FASTA文件｜FASTA of unaligned study sequences |
| `-i, --input` | 必填 |  | 特征表文件(自动识别BIOM/TSV/Excel/Mothur)｜Input table of sequence abundances (auto-detect BIOM/TSV/Excel/Mothur) |
| `-o, --output-dir` | `./picrust2_output` |  | 输出目录｜Output directory (default: ./picrust2_output) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--max-nsti` | `2.0` | float | 最大NSTI阈值｜Maximum NSTI value (default: 2.0) |
| `--stratified` | — | store_true | 生成分层输出表｜Generate stratified output tables |
| `--in-traits` | `EC,KO` |  | 功能数据库(EC,KO,GO,PFAM,BIGG,CAZY)｜Gene families to predict (default: EC,KO) |
| `--placement-tool` | `epa-ng` | epa-ng/sepp | 序列放置工具｜Placement tool (default: epa-ng) |
| `--hsp-method` | `mp` | mp/emp_prob/pic/scp/subtree_average | 隐状态预测方法｜HSP method (default: mp) |
| `--edge-exponent` | `0.5` | float | HSP edge exponent (default: 0.5) |
| `--pipeline` | `auto` | auto/split/single | 流程类型: auto自动检测, split双域, single单参考｜Pipeline type (default: auto) |
| `--min-align` | `0.8` | float | 最小比对比例｜Minimum alignment proportion (default: 0.8) |
| `--min-reads` | `1` | int | 每ASV最小reads数｜Minimum reads per ASV (default: 1) |
| `--min-samples` | `1` | int | 每ASV最小样本数｜Minimum samples per ASV (default: 1) |
| `--no-pathways` | — | store_true | 跳过通路推断｜Skip pathway inference |
| `--coverage` | — | store_true | 计算通路覆盖度｜Calculate pathway coverages |
| `--skip-minpath` | — | store_true | 跳过MinPath｜Skip MinPath |
| `--no-gap-fill` | — | store_true | 跳过gap filling｜Skip gap filling |
| `--per-sequence-contrib` | — | store_true | 逐序列运行MinPath｜Run MinPath per sequence |
| `--skip-norm` | — | store_true | 跳过marker gene拷贝数归一化｜Skip normalization by marker gene copies |
| `--remove-intermediate` | — | store_true | 移除中间文件｜Remove intermediate files |
| `--verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- PICRUSt2（conda 环境 `picrust_v.2.6.3`，脚本 `picrust2_pipeline.py` 与 `picrust2_pipeline_singleRef.py`，默认 `~/miniforge3/envs/picrust_v.2.6.3/bin/`）
- add_descriptions.py（给功能表加描述，走同一 conda 环境）
- Python 包：biom、pandas（输入预处理：BIOM 读取、Excel 转换）
- 参考数据库：PICRUSt2 自带的参考树与基因内容数据库（随环境安装）

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持。检测到 `04_pathway/path_abun_unstrat.tsv.gz`（或跳过通路时的 `03_metagenome/pred_metagenome_unstrat.tsv.gz`）即跳过整个流程。换参数重跑前先删旧产物。

**Q2：为什么报「BIOM 表与 FASTA ID 不匹配」？**
特征表行列方向放反了。本工具会自动尝试转置并比对 ID，转置后仍不匹配才报错——这时请检查 FASTA header 与特征表观测 ID 是否真的同一套命名。

**Q3：FASTA header 带描述信息行吗？**
可以。程序会自动清理，只保留第一个空格前的 ID，确保与特征表一致。

**Q4：Excel 特征表怎么传？**
直接 `-i 表.xlsx` 即可自动转 TSV；.xls 需额外装 xlrd 库（建议另存为 .xlsx）。

**Q5：max-nsti 是什么，要改吗？**
NSTI>2.0 的序列会被排除出宏基因组预测（它们没有可信的近亲）。样本整体 NSTI 偏高时可调大，但要接受预测精度下降。**默认 2.0 一般不动。**