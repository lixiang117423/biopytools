# Swave 结构变异检测 | Swave Structural Variant Detection

一句话理解：**从泛基因组图（GFA）里找出样本之间的结构变异（大片段缺失/插入/倒位/重复等），并把图上用编号表示的路径还原成真实碱基序列**，供下游分析使用。

## 功能概述 | Overview { #overview }

- 主命令 `call`：从 minigraph / cactus / pggb 构建的泛基因组图检测结构变异（SV）
- 检测后自动把 VCF 中的图路径编号转成真实 REF/ALT 碱基序列（`*.converted.vcf`）
- 支持多种 GFA 来源，cactus/pggb 需额外提供 decomposed VCF
- 附带转换/提取子命令：convert-seq、convert-plines、extract-csv、extract-sample
- `pav` 子命令：从 converted VCF 提取 PAV（存在/缺失）矩阵
- call 与序列转换均支持断点续传

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools swave call -i assemblies.tsv -r ref.fa -g graph.gfa -s minigraph
```

最小输入：样本组装 TSV + 参考基因组 FASTA + 泛基因组图 GFA，并指定图的来源。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 泛基因组图 | 把多个样本的基因组「叠」在一起画成一张图，共同部分共用一条路，差异部分分叉 |
| GFA | 图的一种标准文本格式，记录节点和它们怎么连 |
| 结构变异(SV) | 大片段差异，如整段缺失、插入、倒位、重复，区别于单碱基 SNP |
| snarl | 图上结构复杂、需要单独处理的分叉「枢纽」 |
| Decomposed VCF | 预先拆分好的变异结果，cactus/pggb 这类图构建工具的配套产物 |
| PAV 矩阵 | 存在/缺失矩阵，行=变异，列=样本，1=有、0=没有 |

## 输入 | Input { #input }

### 样本组装 TSV

制表符分隔，第一行表头（含 NAME 关键字），之后每行一个样本：第一列样本名，后续列为该样本的组装 FASTA 路径（可为多个）。工具会自动把 TSV 里的相对路径转成绝对路径（基于 TSV 所在目录），避免 Swave 内部切换工作目录后找不到文件。

```text
NAME	assembly
sample1	/path/to/sample1.fa
sample2	/path/to/sample2.fa
```

### 参考基因组 FASTA

参考基因组序列文件，序列名需与 VCF contig 名一致（不一致时工具会自动重命名一份临时 FASTA 适配）。

### GFA 文件与来源

泛基因组图 GFA 文件，`-s/--gfa-source` 指定来源：`minigraph` / `cactus` / `pggb`。cactus/pggb 还需 `--decomposed-vcf`。

## 参数说明 | Parameters { #parameters }

### call 必需输入 | Required (call)

**通俗理解|In plain words:** 四个必填：样本组装表、参考基因组、图文件、图的来源。来源决定检测细节，cactus/pggb 还额外要 decomposed VCF。

相关参数：`-i/--assemblies-tsv`、`-r/--ref-fasta`、`-g/--gfa-file`、`-s/--gfa-source`、`--decomposed-vcf`。

### 路径与输出 | Paths & output

**通俗理解|In plain words:** `--swave-path` 是 Swave 软件本体的安装目录（默认 `~/software/swave/Swave-main`），一般不用改；`-o` 是输出目录（默认 `./swave_output`）。

相关参数：`--swave-path`（默认 `~/software/swave/Swave-main`）、`-o/--output-dir`（默认 `./swave_output`）。

### SV 检测参数 | SV detection

**通俗理解|In plain words:** `--min-sv-size`/`--max-sv-size` 划定「多长算 SV」的范围（默认 50 bp 到 1 Mb），设太小会掺入短片段噪音，设太大漏掉大片段；`--max-sv-comps` 限制单个 SV 最多由多少个图组件拼成，越大越能还原复杂结构但更慢。**一般用默认即可。**

相关参数：`--min-sv-size`（默认 50）、`--max-sv-size`（默认 1000000）、`--max-sv-comps`（默认 5）。

### 处理选项与样本 | Processing options & samples

**通俗理解|In plain words:** 三个布尔开关按需开：`--dup-to-ins` 把重复报告成插入（下游工具只认插入时用）；`--remove-small` 删掉小于 `--min-sv-size` 的图节点；`--force-reverse` 强制检测反向映射的 snarl（更全但更慢）。`--output-mode` 选输出模式（auto 自动判断群体/单样本）；`--spec-samples` 只对指定样本出结果。

相关参数：`--dup-to-ins`、`--remove-small`、`--force-reverse`、`--output-mode`（默认 auto）、`--spec-samples`。

### 性能与工具路径 | Performance & tool paths

**通俗理解|In plain words:** `-t` 线程数默认 12；`--minigraph-path`/`--gfatools-path` 只在系统找不到这两个工具时指定。

相关参数：`-t/--threads`（默认 12）、`--minigraph-path`（默认 minigraph）、`--gfatools-path`（默认 gfatools）。

### 高级定位 | Advanced targeting

**通俗理解|In plain words:** 排查单个位点时，`--spec-snarl`/`--spec-path` 只检测指定的某个 snarl 或某条 path，正常跑全量不要用。

相关参数：`--spec-snarl`、`--spec-path`。

### 其他子命令 | Other subcommands

**通俗理解|In plain words:** `convert-seq` 把 VCF 里图路径编号还原成碱基序列；`convert-plines` 把 VCF 转成 GFA P 行；`extract-csv` 从 VCF 提 INV/DUP 等 CSV；`extract-sample` 只保留指定样本；`pav` 提取存在/缺失矩阵。这些大多把 Swave 原生命令透传，参数即 Swave 的参数。

相关子命令：`convert-seq`（`--vcf-path`、`--gfa-path`、`--ref-path`、`--force-pangenie`）、`convert-plines`（`--gfa-path`、`--vcf-path`、`--ref-vcf-path`、`--force-vg`）、`extract-csv`（`--vcf-path`、`--spec-csv`）、`extract-sample`（`--vcf-path`、`--spec-samples`）、`pav`（`-i/--vcf-file`、`-o/--output-file`、`--min-ac`、`--no-strip-prefix`、`--svtype`）。

## 分析流程 | Pipeline { #pipeline }

```text
call 子命令
    |
    v
调用 Swave.py call(conda 环境 swave_v.1.2)
    |
    v
输出 sample_level.vcf / sample_level.split.vcf
    |
    v
自动 convert_seq(图路径编号 -> 真实碱基序列)
    |
    v
*.converted.vcf
```

## 输出 | Output { #output }

```text
swave_output/
├── swave.sample_level.vcf                    # 主结果：样本水平 SV
├── swave.sample_level.split.vcf              # 拆分后的 SV
├── swave.sample_level.vcf.converted.vcf      # 序列还原后的 VCF
├── swave.sample_level.split.vcf.converted.vcf
├── swave.log                                 # Swave 原生日志
└── 99_logs/
    └── swave_pipeline.log                    # 流程日志
```

（若 VCF 与 FASTA 的染色体名不一致，会临时生成 `ref_renamed.fa` 并自动清理。）

## 结果解读 | Interpreting Results { #interpreting-results }

- `swave.sample_level.vcf`：核心结果，每个 SV 一行，基因型列给出每个样本该 SV 的状态
- `*.converted.vcf`：把图路径编号还原成真实 REF/ALT 序列后的版本，下游凡是要用碱基序列的步骤（如 pangenie）应优先用它
- `swave.log`：Swave 原生输出，若出现 `Skip generating for snarl`，说明上游对部分复杂 snarl 跳过了 dotplot 可视化（Swave 内部 bug），但 VCF 基因型结果仍完整
- SV 数量明显偏少时，检查 `--min-sv-size` 是否过大、`--max-sv-size` 是否过小

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 标准流程：默认参数即可（min 50 bp，max 1 Mb）
- 只关注大片段 SV：调大 `--min-sv-size`（如 500、1000）
- 下游用 pangenie：检测后用 `--force-pangenie` 的 `convert-seq`（或直接依赖自动转换结果）
- 只要某几个样本：用 `--spec-samples` 减少无关计算
- 单一样本：可设 `--output-mode single`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --assemblies-tsv` | 必填 |  | 样本组装TSV文件｜Assemblies TSV file |
| `-r, --ref-fasta` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-g, --gfa-file` | 必填 |  | 泛基因组图GFA文件｜Pangenome graph GFA file |
| `-s, --gfa-source` | 必填 | minigraph/cactus/pggb | GFA文件来源｜GFA file source (minigraph/cactus/pggb) |
| `--swave-path` | `~/software/swave/Swave-main` |  | Swave软件路径｜Swave software path |
| `-o, --output-dir` | `./swave_output` |  | 输出目录｜Output directory |
| `--decomposed-vcf` | — |  | Decomposed VCF文件(cactus/pggb必需)｜Decomposed VCF file (required for cactus/pggb) |
| `--output-mode` | `auto` | auto/population/single | 输出模式｜Output mode |
| `--spec-samples` | — |  | 指定样本｜Specify samples |
| `--min-sv-size` | `50` | int | 最小SV大小｜Minimum SV size |
| `--max-sv-size` | `1000000` | int | 最大SV大小｜Maximum SV size |
| `--max-sv-comps` | `5` | int | 最大SV组件数｜Maximum number of SV components |
| `--dup-to-ins` | — |  | 将duplication报告为insertion｜Report duplications as insertions |
| `--remove-small` | — |  | 移除小于min_sv_size的节点｜Remove nodes smaller than min_sv_size |
| `--force-reverse` | — |  | 强制调用反向映射snarls｜Force call reversed mapping snarls |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minigraph-path` | `minigraph` |  | minigraph工具路径｜minigraph tool path |
| `--gfatools-path` | `gfatools` |  | gfatools工具路径｜gfatools tool path |
| `--spec-snarl` | — |  | 只调用特定snarl｜Only call specific snarl |
| `--spec-path` | — |  | 只调用特定path｜Only call specific path |
| `--vcf-path` | 必填 |  | VCF文件路径｜VCF file path |
| `--gfa-path` | 必填 |  | GFA文件路径｜GFA file path |
| `--ref-path` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--output-path` | — |  | 输出路径｜Output path |
| `--force-pangenie` | — |  | 强制输出满足pangenie要求的序列｜Force output sequences to meet pangenie requirements |
| `--ref-vcf-path` | — |  | 参考VCF文件路径｜Reference VCF file path |
| `--force-vg` | — |  | 强制输出满足vg要求的序列｜Force output sequences to meet vg requirements |
| `--spec-csv` | — | INV/DUP/All | 特定CSV类型｜Specific CSV type (INV/DUP/All) |
| `-i, --vcf-file` | 必填 |  | swave converted VCF文件｜swave converted VCF file |
| `-o, --output-file` | `pav_matrix.tsv` |  | 输出TSV文件｜Output TSV file |
| `--min-ac` | `1` | int | 最小等位基因数｜Minimum allele count |
| `--no-strip-prefix` | — |  | 保留CHROM中的样本前缀｜Keep sample prefix in CHROM |
| `--svtype` | — |  | 仅保留指定SV类型｜Only keep specified SV types (e.g., DUP INS DEL) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --assemblies-tsv` | 必填 |  | 样本组装TSV文件｜Assemblies TSV file |
| `-r, --ref-fasta` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-g, --gfa-file` | 必填 |  | 泛基因组图GFA文件｜Pangenome graph GFA file |
| `-s, --gfa-source` | 必填 | minigraph/cactus/pggb | GFA文件来源｜GFA file source (minigraph/cactus/pggb) |
| `--swave-path` | `~/software/swave/Swave-main` |  | Swave软件路径｜Swave software path |
| `-o, --output-dir` | `./swave_output` |  | 输出目录｜Output directory |
| `--decomposed-vcf` | — |  | Decomposed VCF文件(cactus/pggb必需)｜Decomposed VCF file (required for cactus/pggb) |
| `--output-mode` | `auto` | auto/population/single | 输出模式｜Output mode |
| `--spec-samples` | — |  | 指定样本｜Specify samples |
| `--min-sv-size` | `50` | int | 最小SV大小｜Minimum SV size |
| `--max-sv-size` | `1000000` | int | 最大SV大小｜Maximum SV size |
| `--max-sv-comps` | `5` | int | 最大SV组件数｜Maximum number of SV components |
| `--dup-to-ins` | — | store_true | 将duplication报告为insertion｜Report duplications as insertions |
| `--remove-small` | — | store_true | 移除小于min_sv_size的节点｜Remove nodes smaller than min_sv_size |
| `--force-reverse` | — | store_true | 强制调用反向映射snarls｜Force call reversed mapping snarls |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--minigraph-path` | `minigraph` |  | minigraph工具路径｜minigraph tool path |
| `--gfatools-path` | `gfatools` |  | gfatools工具路径｜gfatools tool path |
| `--spec-snarl` | — |  | 只调用特定snarl｜Only call specific snarl |
| `--spec-path` | — |  | 只调用特定path｜Only call specific path |
| `--vcf-path` | 必填 |  | VCF文件路径｜VCF file path |
| `--gfa-path` | 必填 |  | GFA文件路径｜GFA file path |
| `--ref-path` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--output-path` | — |  | 输出路径｜Output path |
| `--force-pangenie` | — | store_true | 强制输出满足pangenie要求的序列｜Force output sequences to meet pangenie requirements |
| `--ref-vcf-path` | — |  | 参考VCF文件路径｜Reference VCF file path |
| `--force-vg` | — | store_true | 强制输出满足vg要求的序列｜Force output sequences to meet vg requirements |
| `--spec-csv` | — | INV/DUP/All | 特定CSV类型｜Specific CSV type (INV/DUP/All) |
| `-i, --vcf-file` | 必填 |  | swave converted VCF文件｜swave converted VCF file |
| `-o, --output-file` | `pav_matrix.tsv` |  | 输出TSV文件｜Output TSV file |
| `--min-ac` | `1` | int | 最小等位基因数｜Minimum allele count (default: 1) |
| `--no-strip-prefix` | — | store_true | 保留CHROM中的样本前缀｜Keep sample prefix in CHROM |
| `--svtype` | — |  | 仅保留指定SV类型｜Only keep specified SV types (e.g., DUP INS DEL) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Swave 软件本体（默认 `~/software/swave/Swave-main`，须含 `Swave.py`）
- conda 环境 `swave_v.1.2`（提供带 pysam 的 Python；可用环境变量 `SWAVE_CONDA_ENV` 覆盖）
- minigraph（minigraph 来源必需；推荐 conda 环境 `pan`）
- gfatools（可选）
- Python pysam（contig 名检测与 FASTA 重命名用）

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
会（部分）。`call` 主输出 `swave.sample_level.vcf` 已存在时跳过 call；序列转换按 VCF 逐一判断 `*.converted.vcf` 是否存在并跳过。换参数重跑前需先删旧产物。

**Q2：cactus/pggb 报错要求 decomposed-vcf？**
cactus/pggb 图的检测需要预先拆分好的 VCF，用 `--decomposed-vcf` 传入。

**Q3：为什么 swave.log 里出现 Skip generating for snarl？**
Swave 对部分复杂 snarl 生成 dotplot 时触发 numpy 越界/NaN 而跳过（上游 bug），属逐 snarl 优雅降级，VCF 结果完整，仅缺失这些位点的可视化。

**Q4：报「找不到 Swave conda 环境」？**
代码固定找 conda 环境 `swave_v.1.2`，若你的环境名不同，用环境变量 `SWAVE_CONDA_ENV=你的环境名` 覆盖。

**Q5：TSV 里写了相对路径会找不到文件吗？**
不会，工具会先把 TSV 里的相对路径基于 TSV 所在目录转成绝对路径再执行。
