# PSVCP 线性泛基因组构建 | PSVCP Linear Pangenome Construction

一句话理解：把多个基因组的「独有插入片段」按顺序拼进一个参考基因组里，得到一条不断变长的「线性泛基因组」，并标出每个插入片段来自哪个基因组。

## 功能概述 | Overview { #overview }

- 输入一个参考基因组 + 若干 query 基因组（各配 FASTA + GFF 注释），逐轮并入
- 每轮用 nucmer（MUMmer）比对 + Assemblytics 检测结构变异，把 query 相对参考的 >50bp 插入片段拼进参考
- 输出一条线性泛基因组 `pan.fa` + 注释 `pan.gff` + 插入片段（PAV）清单
- 附带生成 PAV 信息表与「样本 × PAV」0/1 来源矩阵
- 支持整体与逐轮两级断点续传，`--force` 可强制重跑

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools psvcp -i genome_gff_dir/ -l genome_list.txt -o pangenome_out/
```

最小输入：一个含 `{name}.fa + {name}.gff/.gff3` 的目录，加一个 genome_list（第 1 行参考、其余 query）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 线性泛基因组 | 一条「代表全物种」的序列，把各基因组的独有片段都塞进同一条链里 |
| PAV（存在/缺失变异） | 某段序列「这个基因组有、那个基因组没有」的差异 |
| 参考基因组 / query | 第 1 个基因组当「底稿」（参考），其余基因组逐个「贡献」独有片段（query） |
| nucmer / MUMmer | 基因组两两比对的经典工具，找出两条序列哪里对应得上 |
| Assemblytics | 读 nucmer 结果，把结构变异（插入/缺失等）一条条列出来 |
| 插入（Insertion）/ 缺失（Deletion） | 相对参考，query 多出来的一段是「插入」，少掉的一段是「缺失」 |
| 断点续传 | 已算完的步骤下次重跑自动跳过，中断了不用从头再来 |

## 输入 | Input { #input }

### genome 目录（-i/--genome-dir）

一个目录，每个基因组放两个文件，且**文件名必须同名**（去 `.fa` 后缀）：

```text
genome_gff_dir/
├── ref.fa       # 参考基因组序列
├── ref.gff      # 参考注释（.gff 优先于 .gff3）
├── q1.fa
├── q1.gff3
├── q2.fa
└── q2.gff
```

### genome_list（-l/--genome-list）

纯文本，每行一个基因组文件名；**第 1 行是参考，其余按顺序作为 query 逐轮并入**，顺序即最终并入顺序：

```text
ref.fa
q1.fa
q2.fa
```

至少需要 2 个基因组（1 参考 + 至少 1 query），程序会校验每个名字都有对应的 `.fa` 与 `.gff/.gff3`。

## 参数说明 | Parameters { #parameters }

### 必需输入 | Required input

**通俗理解|In plain words:** `-i` 指「基因组文件放哪个目录」，`-l` 指「按什么顺序、谁当参考」。genome_list 的顺序很关键：第一行就是底稿，后面的会一个接一个往里面拼，换顺序结果会不同。

### 输出与运行 | Output & runtime

**通俗理解|In plain words:** `-o` 是结果放哪，默认 `~/psvcp_out`。`-t` 传给 nucmer 做并行，调大通常更快，默认 12 够用。`--force` 是「忽略断点续传、强制从头重跑」，只有换参数或想彻底重算时才用。`--log-file` 与 `-v` 管日志位置和详细程度，一般不用动。

### 内部固定阈值（不暴露）| Fixed internal thresholds

**通俗理解|In plain words:** 这些是沿用 PSVCP 原始脚本的固定值，不上命令行：只有 >50bp 的插入才会并入（`MIN_INSERTION_SIZE=50`）；Assemblytics 只报 50bp 到 10Mb 之间、且唯一长度 1000bp 以上的变异；nucmer 用 `--maxgap 500 --mincluster 1000 --diagdiff 20`。正常使用无需关心，记在 `software_versions.yml` 里备查。

## 分析流程 | Pipeline { #pipeline }

```text
参考 + query 列表
  -> 预检：assemblytics 依赖的 numpy 是否可用
  -> 初始化 ref0 = 参考（符号链接）
  -> 逐轮（对每个 query）：
       1. samtools faidx 建索引
       2. nucmer 比对 ref vs query
       3. Assemblytics 检测结构变异
       4. 提取插入片段坐标
       5. 用 query 序列拼接插入片段进 ref
       6. 更新 GFF 注释与 PAV 清单
       7. ref{N-1} + query -> ref{N}
  -> 终化：pan.fa / pan.gff / pan.pav.gff 链接 + 排序
  -> 生成 PAV 信息表与 0/1 矩阵
```

## 输出 | Output { #output }

```text
output_dir/
├── pan.fa                    # 最终线性泛基因组序列
├── pan.gff                   # 最终注释
├── pan.pav.gff               # 插入片段（PAV）清单，未排序
├── pan.pav.sorted.gff        # 按位置排序后的 PAV
├── pan.pav.info.tsv          # PAV 信息表：pan 区间/来源/原基因组位置/长度
├── pan.pav.matrix.tsv        # 样本 × PAV 的 0/1 来源矩阵
├── 00_pipeline_info/
│   └── software_versions.yml # 软件版本与运行参数记录
├── pan_dir_result/           # 中间工作目录
│   ├── ref0.fa / ref0.gff    # 参考基因组（符号链接）
│   ├── ref{N}.fa/.gff/.pav.gff  # 每轮并入后的累积参考
│   ├── {ref}{q}.fa/.gff/.pav.gff # 每轮链式产物
│   └── {ref}_{q}/             # 每轮工作子目录（nucmer/assemblytics 中间文件）
└── psvcp.log                 # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. pan.fa / pan.gff

**通俗理解|In plain words:** 这是最终成品——一条包含所有并入片段的线性泛基因组序列及其注释。下游分析（比对、注释、比较）都从这里出发。

### 2. pan.pav.gff（及其 sorted 版）

**通俗理解|In plain words:** 泛基因组里每个「独有插入片段」的清单，一行一个 PAV，记录它在 pan 上的区间和来自哪个样本。第 9 列 `ID=` 形如 `样本_染色体_起点`。

### 3. pan.pav.info.tsv

每行一个 PAV 的信息表，列：`pav_id`（pan 上的区间）、`pan_chr/pan_start/pan_end`、`length_bp`（长度）、`source`（来源样本）、`orig_chr/orig_start/orig_end`（在原基因组上的位置）。

### 4. pan.pav.matrix.tsv

样本 × PAV 的 0/1 矩阵。**注意语义**：1 = 该样本贡献了这个插入（来源），0 = 其余样本（含参考及未参与该轮比对的 query）。这是「来源标注」矩阵，不是样本间的真实共享性（PSVCP 原始流程只做 query vs 参考两两比对，未做样本两两比较）。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- genome_list 顺序即并入顺序，把最完整、质量最高的基因组放第 1 行当参考
- 换参数（如换了 query 顺序）重跑旧输出目录前，用 `--force`，否则会因 `pan.fa` 已存在而直接跳过
- `-t` 主要影响 nucmer 速度，核多可调大
- 内部阈值（插入 50bp、Assemblytics 50~10Mb 等）沿用原始设定，一般不改

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-dir` | 必填 | Path | genome 目录(含 {name}.fa + {name}.gff/.gff3)｜genome dir with {name}.fa + {name}.gff/.gff3 |
| `-l, --genome-list` | 必填 | Path | genome_list 文本(行1=ref,其余=query)｜genome_list (line1=ref, rest=queries) |
| `-o, --output-dir` | `~/psvcp_out` | Path | 输出目录｜output directory |
| `-t, --threads` | `12` | int | 线程数｜threads |
| `--force` | — |  | 忽略断点续传,强制重跑｜ignore checkpoint, force rerun |
| `--log-file` | — | Path | 日志文件路径(默认 output_dir/psvcp.log)｜log file path |
| `-v, --verbose` | — |  | 详细输出｜verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome-dir` | 必填 |  | genome 目录(含 {name}.fa + {name}.gff/.gff3)｜genome dir with {name}.fa + {name}.gff/.gff3 |
| `-l, --genome-list` | 必填 |  | genome_list 文本(行1=ref,其余=query)｜genome_list (line1=ref, rest=queries) |
| `-o, --output-dir` | `~/psvcp_out` |  | 输出目录｜output directory (default: ~/psvcp_out) |
| `-t, --threads` | `12` | int | 线程数｜threads (default: 12) |
| `--force` | — | store_true | 忽略断点续传,强制重跑｜ignore checkpoint, force rerun |
| `--log-file` | — |  | 日志文件路径(默认 output_dir/psvcp.log)｜log file path |
| `-v, --verbose` | — | store_true | 详细输出｜verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- nucmer（MUMmer，conda 环境 `pan`）
- assemblytics（conda 环境 `psvcp_v.1.0.1`，依赖 numpy）
- bedtools、samtools（conda 环境 `align`）
- Rscript、python3（conda 环境 `psvcp_v.1.0.1`，供 vendored helper 使用）

## 常见问题 | FAQ { #faq }

**Q1：报「env python cannot import numpy」？**
assemblytics 依赖 numpy，若 `psvcp_v.1.0.1` 环境的 python 被 GraalPy 顶掉会导致 numpy 导入失败。程序会明确报错并给修复命令：`conda remove -n psvcp_v.1.0.1 graalpy && conda install -n psvcp_v.1.0.1 python=3.10 numpy`。

**Q2：断点续传怎么判断？**
两级：整体上若 `pan.fa` 已存在且非空则跳过整个流程；逐轮上若该轮 `ref{N}.fa` 已存在则跳过该轮。重跑想全部重算用 `--force`。

**Q3：为什么报「genome_list 至少需要 2 个基因组」？**
必须 1 个参考 + 至少 1 个 query 才有东西可并入，缺 query 时流程无意义，程序在校验阶段报错。

**Q4：为什么报某个基因组「gff 不存在」？**
每个 genome_list 里的名字都要在 genome 目录下有同名 `.fa` 和 `.gff`（或 `.gff3`）。注意文件名要完全一致（仅后缀不同），且程序优先找 `.gff`。

**Q5：pan.pav.matrix.tsv 里的 0 表示「没有这个片段」吗？**
不完全是。矩阵是「来源标注」：只有贡献该插入的样本标 1，其余一律 0（含参考，以及没参与该轮比对的 query）。它不等同于样本间真实的 presence/absence。
