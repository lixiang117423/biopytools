# Pathorepeat - 病原菌重复序列注释 | Pathogen Repeat Annotation

一句话理解:**给病原菌(卵菌、原生生物等非植物物种)的基因组自动注释转座子**——RepeatModeler 从头建重复库,RepeatMasker 软屏蔽(序列信息不丢),TEsorter 给每个重复家族做精确分类,一条命令拿到「库 + 屏蔽基因组 + 分类报告」。它是 EDTA(植物专用)的病原菌替代品。

## 功能概述 | Overview { #overview }

- RepeatModeler2 从头建库(默认 `-LTRStruct` 提高 LTR 分类准确性)
- RepeatMasker 默认 `-xsmall` 小写软屏蔽——下游 SignalP/RxLR motif 扫描、基因预测序列信息不丢
- TEsorter 蛋白域分类(家族→Order/Superfamily/Clade),补足 RepeatModeler 对 unknown 家族分类不足
- 可选 effector 候选区交叉检查:输出 repeat↔effector overlap 逐条清单,供人工核查 edge case
- 每类 repeat 的 bp/占比/GC% 统计——服务 two-speed genome 中 AT-rich 重复区识别
- 支持单基因组或文件夹批量;断点续传(样品×步骤粒度);TEsorter 失败自动降级为全 unknown 继续出报告

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools pathorepeat -i genome.fa -o out_dir/
```

最小输入:一个基因组 FASTA。批量:`biopytools pathorepeat -i genome_dir/ -o out_dir/`(文件夹下放多个 `.fa/.fna/.fasta`)。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 重复序列/转座子(TE) | 基因组里大量反复出现的"复制粘贴"片段,病原菌里常聚集成"重复丰富区" |
| 从头建库(de novo) | 不参考数据库,直接从这份基因组归纳出重复序列"代表模板库" |
| 软屏蔽(soft/xsmall) | 重复区写成小写字母,序列还在、信息不丢;硬屏蔽(N)则是把序列抹掉——病原菌必须软屏蔽,否则下游 SignalP/RxLR motif 扫描会丢序列 |
| -LTRStruct | RepeatModeler 的"LTR 结构精查"模式,像给 LTR 逆转座子做"指纹比对",分类更准,代价是更慢 |
| TEsorter | 给每个重复家族按"核心蛋白域"归类(像按发动机型号给汽车分类),比序列相似度分类更可靠 |
| two-speed genome | 病原菌常见现象:一部分染色体区域"跑得快"(重复多、AT 含量高、效应子基因聚集),另一部分保守;GC% 列就是帮你找这些"快区" |
| effector 交叉检查 | 把效应子候选区与重复区间求交集,列出落在重复区的效应子位点清单,供人工核查(不自动删改任何东西) |

## 输入 | Input { #input }

### 基因组 FASTA(单文件或文件夹)

```text
>chr1
ATGCGTACGTACGTAGCTA...
```

- 单文件:`.fa/.fna/.fasta`;文件夹批量:目录下放多个 FASTA(不递归子目录)
- 样品名 = 文件名去掉扩展名(`Psojae.fa` → `Psojae`);`.gz` 需先解压
- 批量模式下不同文件不能只有后缀不同(`a.fa` 与 `a.fasta` 视为重名报错)

### effector 候选区(可选,BED 或 GFF3)

BED(染色体、起点、终点,0-based)或 GFF3(取每个 feature 的区间,ID/Name 作为名称)。来自 `biopytools phyto_effector` 或人工整理均可。

## 参数说明 | Parameters { #parameters }

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 基因组FASTA或文件夹(批量)｜Genome FASTA or directory (batch) |
| `--output-dir, -o` | `./pathorepeat_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Thread count |
| `--masking-mode` | `xsmall` | xsmall/soft/hard/x | 屏蔽模式(xsmall=小写软屏蔽,病原菌默认)｜Masking mode (xsmall=lowercase soft mask, pathogen default) |
| `--ltr-struct/--no-ltr-struct` | `True` |  | RepeatModeler -LTRStruct(默认开)｜-LTRStruct (default on) |
| `--tesorter-db` | `rexdb` | gydb/rexdb/rexdb-plant/rexdb-metazoa/rexdb-v3/rexdb-plantv3/rexdb-metazoav3/rexdb-pnas/rexdb-line/sine | TEsorter数据库(REXdb植物/动物为主,卵菌可试gydb)｜TEsorter db |
| `--db-hmm` | — |  | 自定义TEsorter HMM文件(优先于--tesorter-db)｜Custom HMM file |
| `--famdb-dir` | — | Path | Dfam famdb数据目录(注入FAMDB_DIR启用RM2自带分类;不设则分类失败自动降级)｜Dfam famdb dir (injected as FAMDB_DIR; auto-degrades if unset) |
| `--effector-bed` | — |  | effector候选区BED(仅单文件模式)｜Effector BED (single-sample) |
| `--effector-gff` | — |  | effector候选区GFF3(仅单文件模式)｜Effector GFF3 (single-sample) |
| `--genome-name` | — |  | 输出前缀(仅单文件模式)｜Output prefix (single-sample only) |
| `--skip-completed/--no-skip-completed` | `True` |  | 断点续传(跳过已完成步骤)｜Resume (skip completed steps) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组FASTA或文件夹(批量)｜Genome FASTA or directory (batch) |
| `-o, --output-dir` | `./pathorepeat_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Thread count |
| `--masking-mode` | `xsmall` | xsmall/soft/hard/x | 屏蔽模式(xsmall=小写软屏蔽,保留序列信息,病原菌默认)｜Masking mode (xsmall=lowercase soft mask, default for pathogens) |
| `--ltr-struct` | `True` | store_true |  |
| `--no-ltr-struct` | — | store_false | 关闭RepeatModeler -LTRStruct(更快但LTR建库变差)｜Disable -LTRStruct (faster but worse LTR library) |
| `--tesorter-db` | `rexdb` |  | TEsorter数据库(REXdb植物/动物lineage为主,卵菌/原生生物可试gydb)｜TEsorter db (REXdb is plant/metazoa-heavy; try gydb for oomycetes/protists) |
| `--db-hmm` | — |  | 自定义TEsorter HMM文件(优先于--tesorter-db)｜Custom TEsorter HMM file (overrides --tesorter-db) |
| `--famdb-dir` | — |  | Dfam famdb数据目录(含famdb.py与*.h5;设置后注入FAMDB_DIR,启用RM2自带分类;不设则分类失败时自动降级)｜Dfam famdb dir (famdb.py + *.h5; injected as FAMDB_DIR to enable RM2 classification; auto-degrades if unset) |
| `--effector-bed` | — |  | effector候选区BED(仅单文件模式)｜Effector regions BED (single-sample mode only) |
| `--effector-gff` | — |  | effector候选区GFF3(仅单文件模式)｜Effector regions GFF3 (single-sample mode only) |
| `--genome-name` | — |  | 输出文件前缀(仅单文件模式,默认basename剥后缀)｜Output prefix (single-sample only) |
| `--skip-completed` | `True` | store_true |  |
| `--no-skip-completed` | — | store_false | 全部重跑(忽略已完成步骤)｜Rerun all (ignore done steps) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

**通俗理解|In plain words(分组导读):**

- 流程控制组:一般不用动。`--no-ltr-struct` 只在赶时间时用(牺牲 LTR 精度);`--no-skip-completed` 用于改参数后强制全部重跑。
- 屏蔽模式组:**病原菌保持默认 xsmall**。hard(N)只在你确定下游工具完全不用序列内容时才选。
- TEsorter 数据库组:`rexdb` 默认;**REXdb 以植物/动物 lineage 为主,卵菌/原生生物可能大量家族分不进去**——发现 unknown 比例高时先试 `gydb`,仍不行用 `--db-hmm` 提供自定义 HMM 库。
- effector 组:只接 BED 或 GFF3 其一;仅单基因组模式可用(批量时多基因组染色体名会撞车)。

## 输出 | Output { #output }

```text
out_dir/
├── 00_pipeline_info/software_versions.yml   # 软件版本
├── 01_modeler/                              # 重复库(核心)
│   └── {样品}_rm_run/{样品}_db-families.fa  # de novo 重复库
├── 02_masker/
│   ├── {样品}.masked.fa                     # 小写软屏蔽基因组(下游基因预测用这个)
│   ├── {样品}.out / {样品}.gff              # 逐条重复注释
│   └── {样品}.tbl                           # RepeatMasker 自带汇总
├── 03_tesorter/
│   └── {样品}_db-families.rexdb.cls.tsv     # 家族分类表(核心)
├── 04_summary/
│   ├── {样品}.repeat_summary.tsv            # 总体+逐类 bp/%/GC%/家族数
│   ├── {样品}.families_classified.tsv       # 每个家族的分类(未分类标 unknown)
│   ├── {样品}.effector_overlap.tsv          # effector 交叉检查(可选)
│   └── batch_summary.tsv                    # 批量模式每样品状态
├── 99_logs/                                 # 全量/stdout/stderr 三份日志
└── tmp/                                     # 运行结束清理
```

## 结果解读 | Interpreting Results { #results }

- **repeat_summary.tsv 的 overall 行**:masked_pct 是全基因组重复占比。病原菌参考:P. sojae ~10-20%,根肿菌等小型原生生物低复杂度区多、数值可能偏高;显著高于同属近缘种要警惕组装质量或污染。
- **per_class 行的 gc_pct 列**:某类 repeat GC 显著低于基因组平均 → 该类富集在 AT-rich 区,正是 two-speed genome 的"快区"信号;结合 effector_overlap 看效应子是否落在这些区。
- **families_classified.tsv 的 unknown 行比例**:unknown >50% 属正常(REXdb 覆盖所限),但此时分类结论要保守;先试 `--tesorter-db gydb` 重跑第 3 步。
- **effector_overlap.tsv**:每个落在重复区的 effector 位点都要人工过一遍——repeat 注释边界错了会连带效应子基因模型出错(你之前遇到的错分类就是这么来的);overlapped 比例高不一定是坏事(P. sojae 效应子本来就富集在重复丰富区);`distance_to_nearest_repeat` 列语义为**间隔碱数**(不重叠时两区间之间的碱基数,重叠时为 0)。
- **batch_summary.tsv 的 failed 行**:看 reason 列定位失败步骤,修复后原样重跑,已成功的样品自动跳过。

## 参数选择建议 | Parameter Guidance { #guidance }

| 场景 | 建议 |
|------|------|
| P. sojae(~70Mb,标准) | 全默认;关注 per_class GC% 与 effector_overlap |
| P. brassicae(~25Mb,低复杂度高) | 全默认;`--masking-mode xsmall` 下重复占比偏高属预期 |
| 只有 WGS 没有组装 | 本模块不适用(需要组装);reads 层面重复估计另行处理 |
| 效应子错分类复盘 | `--effector-bed` 喂入候选区,重点人工核查 overlap≥50% 的位点 |
| unknown 太多 | `--tesorter-db gydb` 重跑;或 `--db-hmm` 提供卵菌专用 HMM 库 |

## 依赖 | Dependencies

- `repeat` 域环境:RepeatModeler 2.0.9(含 BuildDatabase)、RepeatMasker 4.2.4、TEsorter 1.5.1
- 可选(启用 RM2 自带 Dfam 分类):Dfam famdb 数据两层配置,缺一不可——①数据目录:从 dfam.org `releases/current/families/FamDB/` 下载基础片+curated 四件(~2.2GB),`gunzip` 解压到同一目录(如 `~/database/dfam`),并软链环境内 `famdb.py` 进去;②站点配置:`famdb.conf`(位于 `~/miniforge3/envs/repeat/share/famdb-3.0.0/`)设 `FAMDB_DATA_DIR = <数据目录>`——`famdb.py` 靠它找 h5(RepeatClassifier 调用不带 `-i`,环境变量无效)。验证:`conda run -n repeat ~/database/dfam/famdb.py info` 应输出 `Version : 4.0`。配置好后默认即可用;`--famdb-dir` 仅用于显式指定非默认的 famdb.py 位置。未配置时分类失败会自动降级为未分类库继续流程
- 运行时自动探测版本写入 `software_versions.yml`;工具缺失时报错并列出可设置的 `*_PATH` 环境变量

## 常见问题 | FAQ { #faq }

- **为什么不用 EDTA?** EDTA 针对植物 TE 组成优化(数据库与分类规则偏植物);卵菌/原生生物 TE 谱系差异大,通用 RepeatModeler2+TEsorter 更稳。
- **为什么默认软屏蔽?** 硬屏蔽把序列换成 N,下游 SignalP/RxLR motif 扫描直接丢序列信息;软屏蔽只是变小写,基因预测器照样识别。
- **RepeatModeler 阶段很慢?** `-LTRStruct` 约占一半耗时;赶时间 `--no-ltr-struct`,但 LTR 家族分类会变粗。
- **重跑会不会推倒重来?** 不会。断点续传按"样品×步骤"跳过已完成步骤;换参数(如换 `--tesorter-db`)后想全部重跑需 `--no-skip-completed` 或删除对应步骤目录。
- **批量模式一个样品失败会影响其他吗?** 不会,失败样品记入 batch_summary.tsv 后继续;最终退出码非 0 提示有失败。
- **日志出现"RepeatModeler 分类失败…降级"或 "Could not determine FamDB version"?** RM2 分类步需要 Dfam famdb 数据,**两层配置缺一不可**:①数据目录(dfam.org `releases/current/families/FamDB/` 下载基础片+curated 四件 ~2.2GB,`gunzip` 到同一目录,软链环境内 `famdb.py` 进去);②`~/miniforge3/envs/repeat/share/famdb-3.0.0/famdb.conf` 里设 `FAMDB_DATA_DIR = <数据目录>`(famdb.py 靠它找 h5;RepeatClassifier 不带 `-i` 调用,此步必做,环境重建后需重设)。验证:`conda run -n repeat ~/database/dfam/famdb.py info` 出 `Version : 4.0` 即通。未配置时模块自动降级:用未分类 consensi.fa 继续 Masker/TEsorter 步骤,9 小时级建模成果不浪费;配置好后想补 RM2 自带分类,删除 `01_modeler/{样品}_rm_run/{样品}_db-families.fa` 重跑(会重跑建模)。
- **"No LTRs identified" 是真的没有 LTR 吗?** 未必。conda 环境下 LTR_retriever 对 RepeatMasker 的 RMblast 依赖检查可能误报("The RMblast engine is not installed"),导致 -LTRStruct 的 LTR 结构验证空手而归;主库与 TEsorter 分类不受影响。需要 LTR 精修可单独用 lai 模块(其 LTR_retriever 链路独立,不经该检查)。
