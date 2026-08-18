# WGDI 比较基因组学分析 | WGDI Comparative Genomics Analysis

一句话理解：做「两个物种基因组的同源关系」三件套——先画同源基因点图看整体结构，再做共线性分析找出对应的基因块，最后算 Ka/Ks 看这些同源基因在演化中受了多大的选择压力，共 3 个子命令。

## 功能概述 | Overview { #overview }

- 封装 WGDI 工具的 3 个子命令：`dotplot`（同源基因点图）、`collinearity`（共线性分析）、`calks`（Ka/Ks 计算）
- `dotplot`：根据 BLAST 结果画两物种同源基因点图，直观看出染色体对应关系与结构变异
- `collinearity`：识别两物种间的共线性区块（同源基因的保守排列）
- `calks`：基于共线性结果计算同源基因对的 Ka/Ks，评估选择压力
- 自动生成 WGDI 配置文件（.conf）并调用 wgdi，无需手写配置

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools wgdi dotplot -b blast.txt --gff1 s1.gff --gff2 s2.gff --lens1 s1.lens --lens2 s2.lens
```

最小输入：一个 BLAST 结果 + 两个物种的 GFF + 两个 LENS 文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 同源基因 | 两物种里「同一个祖先基因」的后代，功能往往相似 |
| 点图（dotplot） | 横轴一个物种、纵轴另一个物种，把同源基因画成点；点连成的对角线条带就是共线性区块 |
| 共线性（collinearity） | 基因在染色体上「排队的顺序」在两个物种间保持一致的现象 |
| Ka / Ks | 非同义替换率 / 同义替换率；Ks 近似「中性演化速度」，Ka/Ks 比值衡量选择压力：<1 纯化选择，≈1 中性，>1 正选择 |
| BLAST 结果 | 用 BLAST 把两物种蛋白/序列互相比对得到的「谁和谁像」的表格 |
| GFF | 基因注释文件，记录每个基因在染色体上的位置 |
| LENS | WGDI 约定的染色体长度文件（记录各染色体长度） |

## 输入 | Input { #input }

按子命令不同（文件需按 WGDI 约定自行准备）：

- `dotplot` / `collinearity`：BLAST 结果（`-b`）+ 物种1/2 的 GFF（`--gff1/--gff2`）+ 物种1/2 的 LENS（`--lens1/--lens2`）
- `calks`：上一步的共线性结果（`-c`）+ 两个物种的 FASTA（`--fasta1/--fasta2`）

## 参数说明 | Parameters { #parameters }

### 通用参数 | Common

**通俗理解|In plain words:** `-o` 是输出目录（WGDI 配置文件与结果都放这里），`-t` 是线程数。`--wgdi-path` 是 WGDI 软件路径，只在你的 wgdi 不在默认位置时才需要给，一般不用动。

### dotplot（同源基因点图）| Dotplot

**通俗理解|In plain words:** 必填 BLAST + 两套 GFF/LENS。`--multiple` 是每个基因保留几个最佳同源，默认 1；`--score/--evalue` 是过滤 BLAST 结果的阈值，默认较宽松，一般不用动；`--position` 决定点按「基因顺序」还是「物理坐标」排，默认 order。`--savefig` 是输出图名，默认 `dotplot.png`。

### collinearity（共线性分析）| Collinearity

**通俗理解|In plain words:** 必填与 dotplot 相同。`--comparison` 选「基因组级」还是「染色体级」比较；`--grading`（红/蓝/灰三级打分）与 `--mg`（最大 gap）控制共线性块识别的松紧，默认值通用，一般不用动；`--process` 是并行进程数，核多可调大。

### calks（Ka/Ks 计算）| CalKs

**通俗理解|In plain words:** 必填共线性结果 + 两个物种 FASTA。输出 `ks_results.txt`，核心是每个同源基因对的 Ks 与 Ka/Ks，用于判断选择压力与估计分化时间。

## 分析流程 | Pipeline { #pipeline }

```text
（建议顺序）
  dotplot:      BLAST + GFF + LENS -> dotplot.png（看整体结构）
  collinearity: BLAST + GFF + LENS -> collinearity.txt（找共线性块）
  calks:        collinearity.txt + FASTA -> ks_results.txt（算 Ka/Ks）
```

## 输出 | Output { #output }

```text
output_dir/
├── dotplot.conf          # 自动生成的 WGDI 配置（dotplot 时）
├── collinearity.conf     # 自动生成的 WGDI 配置（collinearity 时）
├── calks.conf            # 自动生成的 WGDI 配置（calks 时）
├── dotplot.png           # 同源基因点图（dotplot，默认名）
├── collinearity.txt      # 共线性结果表（collinearity，默认名）
├── ks_results.txt        # Ka/Ks 结果表（calks，默认名）
└── wgdi.log              # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **dotplot.png**：横竖两轴各一个物种，每个点是一对同源基因。点连成一条条对角线 = 共线性区块（两物种染色体间的大段对应）；对角线断裂/错位 = 可能的重排或结构变异
- **collinearity.txt**：共线性区块清单，一行记录一个区块的基因对应关系，是 calks 的输入
- **ks_results.txt**：每个同源基因对的 Ks 与 Ka/Ks。Ks 的峰值分布常用来估计全基因组复制（WGD）事件与物种分化时间；Ka/Ks 判断选择压力（<1 纯化、≈1 中性、>1 正选择）
- **.conf 文件**：程序按你给的参数自动写出的 WGDI 配置，可打开核对实际生效的参数

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 先跑 `dotplot` 快速目检两物种结构，再决定要不要做共线性与 Ka/Ks
- 亲缘远、重复多时，`collinearity` 可调小 `--multiple`（取更保守的唯一同源）
- 关注分化时间 / 全基因组复制时，重点看 `ks_results.txt` 的 Ks 峰值
- 文件输入（BLAST/GFF/LENS）需符合 WGDI 约定，一般由上游自行准备

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --blast` | — |  | BLAST结果文件｜BLAST result file |
| `--gff1` | — |  | 物种1 GFF文件｜Species 1 GFF file |
| `--gff2` | — |  | 物种2 GFF文件｜Species 2 GFF file |
| `--lens1` | — |  | 物种1 LENS文件｜Species 1 LENS file |
| `--lens2` | — |  | 物种2 LENS文件｜Species 2 LENS file |
| `-c, --collinearity` | — |  | 共线性分析结果文件｜Collinearity analysis result file |
| `--fasta1` | — |  | 物种1 FASTA文件｜Species 1 FASTA file |
| `--fasta2` | — |  | 物种2 FASTA文件｜Species 2 FASTA file |
| `-o, --output-dir` | `./wgdi_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `8` | int | 线程数｜Number of threads |
| `--wgdi-path` | — |  | WGDI软件路径｜WGDI software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-b, --blast` | 必填 |  | BLAST结果文件｜BLAST result file |
| `--gff1` | 必填 |  | 物种1 GFF文件｜Species 1 GFF file |
| `--gff2` | 必填 |  | 物种2 GFF文件｜Species 2 GFF file |
| `--lens1` | 必填 |  | 物种1 LENS文件｜Species 1 LENS file |
| `--lens2` | 必填 |  | 物种2 LENS文件｜Species 2 LENS file |
| `--genome1-name` | `Genome1` |  | 物种1名称｜Species 1 name (default: Genome1) |
| `--genome2-name` | `Genome2` |  | 物种2名称｜Species 2 name (default: Genome2) |
| `--blast-reverse` | — | store_true | BLAST结果前两列互换｜Swap first two columns of BLAST result |
| `--multiple` | `1` | int | 最佳同源基因数｜Number of best homologous genes (default: 1) |
| `--score` | `100` | int | BLAST分数阈值｜BLAST score threshold (default: 100) |
| `--evalue` | `1e-05` | float | BLAST e-value阈值｜BLAST e-value threshold (default: 1e-5) |
| `--repeat-number` | `10` | int | 重复基因最大数量｜Max number of repetitive genes (default: 10) |
| `--position` | `order` | order/start/end | 基因位置类型｜Gene position type (default: order) |
| `--ancestor-left` | — |  | 左侧物种祖先染色体区域｜Ancestral chromosome region for left species |
| `--ancestor-top` | — |  | 顶部物种祖先染色体区域｜Ancestral chromosome region for top species |
| `--markersize` | `0.5` | float | 点大小｜Marker size (default: 0.5) |
| `--figsize` | `10,10` |  | 图像尺寸｜Figure size (default: 10,10) |
| `--savefig` | `dotplot.png` |  | 输出图像文件｜Output image file (default: dotplot.png) |
| `-o, --output-dir` | `./wgdi_output` |  | 输出目录｜Output directory (default: ./wgdi_output) |
| `-t, --threads` | `8` | int | 线程数｜Number of threads (default: 8) |
| `--comparison` | `genomes` | genomes/chromosomes | 比较模式｜Comparison mode (default: genomes) |
| `--process` | `8` | int | 进程数｜Number of processes (default: 8) |
| `--grading` | `50,40,25` |  | 评分参数(红,蓝,灰)｜Scoring parameters for red,blue,gray (default: 50,40,25) |
| `--mg` | `40,40` |  | 最大gap值｜Max gap values (default: 40,40) |
| `--pvalue` | `1.0` | float | P值阈值｜P-value threshold (default: 1.0) |
| `--savefile` | `collinearity.txt` |  | 输出文件｜Output file (default: collinearity.txt) |
| `-c, --collinearity` | 必填 |  | 共线性分析结果文件｜Collinearity analysis result file |
| `--fasta1` | 必填 |  | 物种1 FASTA文件｜Species 1 FASTA file |
| `--fasta2` | 必填 |  | 物种2 FASTA文件｜Species 2 FASTA file |
| `--wgdi-path` | — |  | WGDI软件路径｜WGDI software path |
| `-v, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- WGDI 软件（自动解析 phylo 域环境并经 conda run 调用，可用 --wgdi-path 或环境变量 WGDI_PATH 覆盖；域环境缺失时回退旧默认路径 `~/miniforge3/envs/phylo/bin/wgdi`）

## 常见问题 | FAQ { #faq }

**Q1：输入文件（BLAST/GFF/LENS）怎么准备？**
本模块只做 WGDI 的封装，输入需符合 WGDI 约定格式自行准备（BLAST 比对结果、两物种 GFF 注释、两物种 LENS 染色体长度文件）。

**Q2：断点续传支持吗？**
不支持。WGDI 三个子命令都是单次调用，输出已存在时不会自动跳过，重跑会覆盖；需要保留旧结果请先移走或改名。

**Q3：.conf 文件是什么？**
程序根据你的命令行参数自动生成的 WGDI 配置文件，再交给 `wgdi -d/-icl/-ks` 执行。可打开核对实际传给 WGDI 的参数。

**Q4：怎么知道点图里哪些是共线性？**
点图中沿对角线方向连成一条条「带」的点群就是共线性区块；散乱无规则的点通常是无同源关系的噪声或重复序列。

**Q5：calks 需要什么前置？**
需要先有 `collinearity` 的共线性结果文件（`collinearity.txt`）加上两个物种的 FASTA 序列，三者都必填。
