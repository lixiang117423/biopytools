# VG 变异图分析 | VG Variation Graph Analysis

一句话理解：一套围绕「变异图」的操作工具——用参考基因组加变异表建图、给图建索引、把测序 reads 比到图上、再从图反向导出变异表，共 4 个子命令。

## 功能概述 | Overview { #overview }

- 封装 vg（variation graph）工具包的 4 个核心子命令：`construct`（建图）、`index`（建索引）、`giraffe`（reads 比对）、`deconstruct`（导出 VCF）
- `construct`：用参考 FASTA + VCF 构建变异图
- `index`：给图建 XG / GCSA / GBWT / Giraffe 索引，比对前需要
- `giraffe`：把测序 reads 比对到图（支持单端/双端，输出 GAM 或 GAF）
- `deconstruct`：从变异图导出 VCF
- 每个子命令若输出文件已存在会自动跳过（断点续传）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg
```

最小输入：参考 FASTA + 一个 VCF。其余子命令见 `biopytools vg --help` 与下方各小节。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 变异图（variation graph） | 用「图」表示参考 + 变异：主线是参考，变异处开叉，比线性参考更能代表群体 |
| 索引 | 给图预先算好的「加速结构」，像书的目录，比对/查询靠它才能快 |
| XG / GCSA / GBWT / Giraffe | 四种不同用途的索引：XG 管坐标与序列，GCSA 管 k-mer 查找，GBWT 管单倍型，Giraffe 是比对专用的一整套（.min/.dist/.xg） |
| reads | 测序仪读出的短序列片段，比对就是把它们「挂」回图上的位置 |
| GAM / GAF | 比对结果的两种格式：GAM 是 vg 二进制格式，GAF 是文本表格格式 |
| 参考路径（reference path） | 图里代表参考基因组的那条主线路径，deconstruct 以它为基准还原变异 |

## 输入 | Input { #input }

按子命令不同：

- `construct`：参考基因组 FASTA（`-r`）+ VCF（`-v`）。VCF 建议用 bgzip 压缩并建索引（`.tbi`），否则会给出警告
- `index`：一个图文件（`-i`，如 `.vg`）
- `giraffe`：一个「已建好索引」的图前缀（`-g`，需已有 `前缀.xg/.min/.dist`）+ reads 文件（`-f`，双端再给 `-f2`）
- `deconstruct`：一个图文件（`-i`）+ 参考路径名（`-r`）

## 参数说明 | Parameters { #parameters }

### 全局参数 | Global

**通俗理解|In plain words:** `--vg-env` 是 vg 所在 conda 环境名，只在环境不叫 `vg_v.1.7.0` 时才改；`--log-level` 管日志详略，默认 INFO 即可，一般不用动。

### construct（建图）| Build graph

**通俗理解|In plain words:** 必填三件套：参考 FASTA、VCF、输出图路径。`-R` 可只取某条染色体/某段区域建图（省时省内存）。`--alt-paths` 额外保留 alt 等位基因路径，需要更完整的等位信息时才开。线程数默认 12 一般够用。

### index（建索引）| Create indexes

**通俗理解|In plain words:** 用 `--xg/--gcsa/--gbwt/--giraffe` 指定要建哪几种索引，至少一种，否则报错。只做比对就 `--giraffe`（会自动顺带建 XG）。`-k` 是 GCSA 的 k-mer 大小，只有建 GCSA 才用到，默认 16 一般不用动。

### giraffe（reads 比对）| Read alignment

**通俗理解|In plain words:** 必须先有 `index --giraffe` 生成的索引前缀。`-f` 给 reads，双端加 `-f2`；`-l/-s` 是片段长度及其标准差，默认 0 自动检测；`--min-identity` 是最低相似度，默认 0 即不过滤；`--format` 选 GAM（默认）还是 GAF 文本。

### deconstruct（导出 VCF）| Export VCF

**通俗理解|In plain words:** 把图还原成 VCF 变异表，必须给参考路径名 `-r`（图里那条参考主线的名字）。`-s` 可只导出指定样本的变异（可重复多次），不指定则全部样本。

## 分析流程 | Pipeline { #pipeline }

```text
典型用法（建图 -> 建索引 -> 比对）：
  construct: 参考 FASTA + VCF -> graph.vg
  index --giraffe: graph.vg -> graph.xg/.min/.dist
  giraffe: reads + graph 前缀 -> alignments.gam
（可选）deconstruct: graph.vg + 参考路径 -> output.vcf
```

## 输出 | Output { #output }

每个子命令输出一个主文件，日志写在当前目录 `vg_<子命令>.log`：

```text
construct   -> graph.vg              # 变异图
index       -> 前缀.xg / .gcsa / .gbwt / .min / .dist   # 各索引
giraffe     -> alignments.gam / .gaf # 比对结果
deconstruct -> output.vcf            # 导出的变异表
vg_<子命令>.log                       # 运行日志（当前目录）
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **graph.vg**：变异图本体，是建索引、比对、导出的共同入口
- **前缀.xg**：坐标与序列索引，deconstruct、可视化等都要用
- **前缀.min / 前缀.dist**：Giraffe 比对专用索引，`.min` 是最小化器、`.dist` 是距离索引，两者缺一不可
- **alignments.gam / .gaf**：reads 比对结果。GAM 是二进制（用 `vg view -a` 查看），GAF 是可直接打开的文本表
- **output.vcf**：从图还原的标准 VCF，可与传统变异分析流程衔接

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只想比对：`index --giraffe` 一步到位（自动含 XG）；要做图查询/坐标操作再补 `--xg`
- 大 VCF / 大基因组建图：用 `-R` 限定区域，显著省时省内存
- 双端测序务必给 `-f2`，并让 `-l/-s` 保持 0 自动检测
- 需要文本格式的比对结果供脚本处理时，giraffe 用 `--format GAF`
- `--vg-env`、`--log-level`、`-t`、`-k` 等默认值一般不用动

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vg-env` | `vg_v.1.7.0` |  | VG conda环境名称｜VG conda environment name (default: vg_v.1.7.0) |
| `--log-level` | `INFO` | DEBUG/INFO/WARN/ERROR | 日志级别｜Log level (default: INFO) |
| `-r, --reference` | 必填 |  | 参考基因组FASTA文件｜Reference FASTA file |
| `-v, --vcf` | 必填 |  | VCF文件｜VCF file |
| `-o, --output` | 必填 |  | 输出VG文件｜Output VG file |
| `-R, --region` | — |  | 指定染色体区域｜Specify chromosome region |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--alt-paths` | — | store_true | 保存alt等位基因路径｜Save alt allele paths |
| `--no-progress` | — | store_true | 不显示进度｜Do not show progress |
| `-i, --input` | 必填 |  | 输入图文件｜Input graph file |
| `--xg` | — | store_true | 创建XG索引｜Create XG index |
| `--gcsa` | — | store_true | 创建GCSA索引｜Create GCSA index |
| `--gbwt` | — | store_true | 创建GBWT索引｜Create GBWT index |
| `--giraffe` | — | store_true | 创建GIRAFFE索引｜Create GIRAFFE indexes |
| `-k, --kmer-size` | `16` | int | GCSA k-mer大小｜GCSA k-mer size (default: 16) |
| `-g, --graph` | 必填 |  | 图文件前缀（索引）｜Graph file prefix (indexed) |
| `-f, --reads` | 必填 |  | 输入reads文件｜Input reads file |
| `-f2, --reads2` | — |  | 第二个reads文件（双端测序）｜Second reads file (paired-end) |
| `-l, --fragment-length` | `0` | int | 片段长度（0=自动检测）｜Fragment length (0=auto) |
| `-s, --fragment-std-dev` | `0` | int | 片段长度标准差（0=自动检测）｜Fragment length std dev (0=auto) |
| `--min-identity` | `0.0` | float | 最小相似度｜Min identity (default: 0.0) |
| `--format` | `GAM` | GAM/GAF | 输出格式｜Output format (default: GAM) |
| `--progress` | — | store_true | 显示进度｜Show progress |
| `-r, --reference-path` | 必填 |  | 参考路径名称｜Reference path name |
| `-s, --samples` | — | append | 样本列表（可多次使用）｜Sample list (can be used multiple times) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- vg（conda 环境 `vg_v.1.7.0`，通过 `--vg-env` 可改）

## 常见问题 | FAQ { #faq }

**Q1：giraffe 报找不到 .min/.dist 文件？**
giraffe 需要已建好的 Giraffe 索引（图前缀 + `.xg/.min/.dist`）。先用 `biopytools vg index -i graph.vg -o 前缀 --giraffe` 建索引，再用同一个前缀跑 giraffe。

**Q2：index 报「必须指定至少一种索引类型」？**
`--xg/--gcsa/--gbwt/--giraffe` 至少要给一个，否则程序不知道建什么。

**Q3：construct 提示 VCF 未压缩/未索引？**
程序会给警告但会继续。建议对 VCF 用 bgzip 压缩并 `tabix` 建索引（`.tbi`），速度和可靠性更好。

**Q4：输出文件已存在会怎样？**
每个子命令都会先检查输出文件，已存在则直接跳过并报告文件大小，不重复计算。想重算需先删除旧输出文件。

**Q5：日志文件在哪？**
日志写在**当前工作目录**的 `vg_<子命令>.log`（如 `vg_construct.log`），不是输出目录。
