# Pixy 群体遗传统计 | Pixy Population Genetics Statistics

一句话理解：**一次性算出 pi、dxy、fst、Watterson's theta、Tajima's D 等群体遗传学核心统计量**，并且能正确纳入「不变位点」，避免 fst 被高估。

## 功能概述 | Overview

- 封装 pixy，计算 pi(核苷酸多样性)、dxy(群体间差异)、fst、Watterson's theta、Tajima's D
- 支持滑窗(`-w`)、BED 窗口(`-b`)、特定位点(`-s`)三种模式
- 能纳入不变位点(invariant sites)，是无偏统计的金标准工具
- 自动检测 VCF 是否含不变位点，缺失时自动加 `--bypass_invariant_check`
- 每个统计量分别调用 pixy 计算，互不影响

## 快速开始 | Quick Start

```bash
biopytools pixy -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000
```

`-w 100000` 表示按 100 kb 滑窗计算；pixy 要求必须指定 `-w`、`-b` 或 `-s` 之一。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| pi(核苷酸多样性) | 群体内随机挑两条序列，平均多少个位点不一样；越高=内部越多样 |
| dxy | 两个群体各挑一条序列的平均差异；衡量群体间「隔多远」 |
| fst | 群体间分化占比打分，0=不分，1=彻底分家 |
| Watterson's theta | 按分离位点数估的群体大小，和 pi 一起看能推断群体历史 |
| Tajima's D | 比较 pi 与 theta 的差异；偏负=近期扩张/选择，偏正=群体收缩 |
| 不变位点 | 所有样本都一样的位点；漏掉它们会让 pi/fst 算出来偏大 |
| 滑窗 | 把染色体切成固定大小的一段段，逐段算统计量 |

## 输入 | Input

### VCF 文件

必须用 bgzip 压缩并建立 tabix 索引(即 `.vcf.gz` + `.vcf.gz.tbi` 或 `.csi`)：

```bash
bgzip variants.vcf && tabix -p vcf variants.vcf.gz
```

### 群体文件

两列：`样本ID` `群体名`：

```text
sample1	popA
sample2	popA
sample3	popB
sample4	popB
```

## 输出 | Output

```text
pixy_output/
├── <前缀>_pi.txt          # pi 结果(按窗口/位点)
├── <前缀>_dxy.txt         # dxy 结果
├── <前缀>_fst.txt         # fst 结果
├── <前缀>_watterson_theta.txt  # Watterson's theta 结果
├── <前缀>_tajima_d.txt    # Tajima's D 结果
└── pixy.log               # 运行日志
```

每个统计量一个结果文件(文件名含统计量名，位于输出目录下)，列含染色体、窗口坐标、对应统计量值。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 各统计量回答不同问题：pi/dxy 看「多样性高低」，fst 看「分得多开」，Tajima's D 看「群体近期是扩张还是收缩」。

- **pi**：群体内核苷酸多样性，值越大内部越多样；可按窗口画曲线找高/低多样区
- **dxy**：群体间平均差异，比 fst 更「绝对值」，不受群体内多样性干扰
- **fst**：0-1，越大分化越大；因为纳入了不变位点，比 VCFtools 的 fst 更准确(不会高估)
- **Tajima's D**：显著为负提示近期扩张或选择清除，显著为正提示群体收缩或平衡选择
- 窗口模式下，某一窗口统计量异常高/低，常对应选择信号或结构变异区域

## 参数选择建议 | Parameter Guidance

- **`-w`/`-b`/`-s` 必须三选一**：常规用 `-w 100000`(100 kb 滑窗)；有基因/外显子等自定义区间用 `-b`；只关心少数位点用 `-s`
- **`--stats`**：默认全算(pi,dxy,fst,watterson,tajima)；只算部分时用 `--stats pi,fst`
- **`--min-samples`/`--max-missing`/`--min-maf`**：默认不限制，**一般不用动**
- **不变位点**：程序自动检测 VCF 是否含不变位点并自动处理，**一般不用手动干预**；确信要强制时用 `--bypass-invariant-check`
- **`-c/--chromosomes`**：只想算某几条染色体时用，逗号分隔
<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径（需用bgzip压缩并建立tabix索引）｜VCF file path (must be bgzip-compressed and tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件路径（两列：样本ID 群体名）｜Population file path (two columns: sample_id population_name) |
| `-o, --output-dir` | `./pixy_output` | Path | 输出目录｜Output directory |
| `--stats` | `pi,dxy,fst,watterson,tajima` |  | 要计算的统计量，逗号分隔｜Statistics to calculate, comma-separated |
| `--calc-pi` | — |  | 计算pi（核苷酸多样性）｜Calculate pi (nucleotide diversity) |
| `--calc-dxy` | — |  | 计算dxy（群体间核苷酸差异）｜Calculate dxy (nucleotide divergence) |
| `--calc-fst` | — |  | 计算fst（遗传分化系数）｜Calculate fst (genetic differentiation) |
| `--calc-watterson-theta` | — |  | 计算Watterson's theta｜Calculate Watterson's theta |
| `--calc-tajima-d` | — |  | 计算Tajima's D｜Calculate Tajima's D |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (null for genome-wide) |
| `-b, --bed-file` | — |  | BED文件定义窗口｜BED file defining windows |
| `-s, --sites-file` | — |  | 位点文件（只计算特定位点）｜Sites file (calculate only specific sites) |
| `--min-samples` | `0` | int | 每个群体最小样本数｜Minimum samples per population |
| `--max-missing` | `1.0` | float | 最大缺失率｜Maximum missing rate |
| `--min-maf` | `0.0` | float | 最小等位基因频率｜Minor allele frequency |
| `--zscore-window` | — | int | Z-score过滤窗口大小｜Z-score filtering window size |
| `-c, --chromosomes` | — |  | 指定染色体列表，逗号分隔｜List of chromosomes, comma-separated |
| `--pixy-path` | `pixy` |  | pixy可执行文件路径｜pixy executable path |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` |  | conda环境路径｜conda environment path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--bypass-invariant-check` | — |  | 强制绕过不变位点检查（默认自动检测VCF并自动绕过）｜Force bypass invariant sites check (default: auto-detect VCF and bypass if needed) |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | 输入VCF文件（需用bgzip压缩并建立tabix索引）｜Input VCF file (must be bgzip-compressed and tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件（两列：样本ID 群体名）｜Population file (two columns: sample_id population_name) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--stats` | `pi,dxy,fst,watterson,tajima` |  | 要计算的统计量，逗号分隔（默认: pi,dxy,fst,watterson,tajima）｜Statistics to calculate, comma-separated (default: pi,dxy,fst,watterson,tajima) |
| `--calc-pi` | — | store_true | 计算pi（核苷酸多样性）｜Calculate pi (nucleotide diversity) |
| `--calc-dxy` | — | store_true | 计算dxy（群体间核苷酸差异）｜Calculate dxy (nucleotide divergence) |
| `--calc-fst` | — | store_true | 计算fst（遗传分化系数）｜Calculate fst (genetic differentiation) |
| `--calc-watterson-theta` | — | store_true | 计算Watterson's theta｜Calculate Watterson's theta |
| `--calc-tajima-d` | — | store_true | 计算Tajima's D｜Calculate Tajima's D |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (null for genome-wide) |
| `-b, --bed-file` | — |  | BED文件定义窗口（自定义大小窗口）｜BED file defining windows (custom-sized windows) |
| `-s, --sites-file` | — |  | 位点文件（只计算特定位点）｜Sites file (calculate only specific sites) |
| `--min-samples` | `0` | int | 每个群体最小样本数（默认: 0=不限制）｜Minimum samples per population (default: 0=no limit) |
| `--max-missing` | `1.0` | float | 最大缺失率（默认: 1.0=不限制）｜Maximum missing rate (default: 1.0=no limit) |
| `--min-maf` | `0.0` | float | 最小等位基因频率（默认: 0.0=不限制）｜Minor allele frequency (default: 0.0=no limit) |
| `--zscore-window` | — | int | Z-score过滤窗口大小（不设置则不过滤）｜Z-score filtering window size (null=no filter) |
| `-c, --chromosomes` | — |  | 指定染色体列表，逗号分隔（不设置则全部）｜List of chromosomes, comma-separated (null for all) |
| `--pixy-path` | `pixy` |  | pixy可执行文件路径（默认: pixy）｜pixy executable path (default: pixy) |
| `--conda-env` | `~/miniforge3/envs/pixy_v.2.0.0` |  | conda环境路径｜conda environment path |
| `-t, --threads` | `12` | int | 线程数（默认: 12）｜Number of threads (default: 12) |
| `--bypass-invariant-check` | — | store_true | 强制绕过不变位点检查（默认自动检测VCF并自动绕过）｜Force bypass invariant sites check (default: auto-detect VCF and bypass if needed) |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- pixy(默认通过 conda 环境 `~/miniforge3/envs/pixy_v.2.0.0` 调用，可用 `--pixy-path`/`--conda-env` 调整)
- bgzip/tabix(htslib，用于 VCF 压缩与索引)

## 常见问题 | FAQ

**Q1：为什么报「必须指定窗口」？**
pixy 要求 `-w`(窗口大小)、`-b`(BED 文件)、`-s`(位点文件)三者至少指定一个。加 `-w 100000` 即可。

**Q2：VCF 必须压缩和索引吗？**
是。pixy 依赖 tabix 索引做随机访问，输入必须是 bgzip 压缩的 `.vcf.gz` 并带 `.tbi`/`.csi` 索引，否则报错。

**Q3：不变位点是什么，为什么要自动加 bypass？**
不变位点(所有样本基因型相同的位点)若不纳入，fst/pi 会系统性偏高。程序会先检测 VCF 是否含不变位点；不含则自动加 `--bypass_invariant_check` 让 pixy 跳过该检查继续运行。

**Q4：换参数重跑要删旧文件吗？**
pixy 模块无断点续传，每次运行重算。换 `-w`、`--stats` 等重跑同一 `-o` 会覆盖同名输出；想保留多组结果请换输出目录。