# LD衰减分析 | Linkage Disequilibrium Decay Analysis (PopLDdecay)

一句话理解：**计算「连锁不平衡(LD)」随距离衰减的速度**——基因上挨得近的位点往往「绑定遗传」，离得远就各走各的；这条下降曲线既能反映群体历史，也能告诉我们做 GWAS 需要多密的标记。

## 功能概述 | Overview

- 封装 PopLDdecay2（输出与经典 PopLDdecay 逐字节一致的加速重写版），计算 r²(或 D')随物理距离的衰减曲线，支持多线程
- 支持 VCF 和 Genotype 两种输入，可指定子群体分别计算并汇总
- 按距离分箱(bin)求平均，自动绘制 LD 衰减曲线图
- 用 Hill & Weir 模型拟合衰减曲线，自动推荐 LD 阈值、背景 r² 与 GWAS 建议窗口
- 支持 MeanBin / MedianBin / PercentileBin / HW 等多种绘图统计方法
- 无断点续传：每次运行直接重算并覆盖同名输出（换参数重跑需换前缀，见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools poplddecay -i variants.vcf -o output
```

`-i` 输入 VCF，`-o` 是输出文件前缀(不是目录)，会生成 `output.stat.gz`、`output.png` 等文件。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 连锁不平衡(LD) | 相邻位点「绑定遗传」的程度：本应随机组合的两个位点，却总是一起出现 |
| r² | 两个位点关联强度的打分，0=完全独立，1=完全绑定 |
| LD 衰减 | 距离越远 r² 越低的下降趋势，像「亲戚关系」，住得越远越不沾亲 |
| 衰减距离 | r² 降到某个阈值时的距离；衰减越快=历史重组越多(群体越古老) |
| bin(分箱) | 把距离相近的位点对归成一格再求平均，让曲线光滑 |
| 背景 r² | 远距离 r² 的「地板」，代表随机噪声水平，正常应接近 0 |
| GWAS 窗口 | 做关联分析时，围绕显著位点向外看多远的区间，由衰减距离决定 |

## 输入 | Input

### VCF / Genotype 文件

标准 VCF(支持 .vcf/.vcf.gz，需含 GT 基因型)或 Genotype 格式，用 `-t/--type` 指定。

### 子群体文件(可选，`-s/--subpop`)

两列制表符分隔：`样本ID` 与 `群体名`。程序按群体名分组，先算全体样本、再逐个群体算、最后合并对比：

```text
sample1	popA
sample2	popA
sample3	popB
sample4	popB
```

每个样本名必须与 VCF 头部的样本名完全一致，否则该样本被忽略；程序会先验证文件再计算。

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入 VCF 或 Genotype 文件 |
| `-o, --output` | 输出文件前缀（不是目录） |

### 位点过滤 | Site filtering

**通俗理解|In plain words:** 决定「什么样的位点信得过、参与计算」。MAF 太低=少数派等位基因太少、区分度低；杂合率太高=数据可能有问题；缺失率太高=依据不足。**绝大多数项目用默认值即可，几乎不需要动。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-m, --min-maf` | 0.005 | 最小等位基因频率，低于此值的位点剔除 |
| `--max-het` | 0.88 | 最大杂合率，高于此值的位点剔除 |
| `--max-miss` | 0.25 | 最大缺失率，高于此值的位点剔除 |
| `-d, --max-dist` | 300 | 只计算多少 kb 以内的位点对 |

### 子群体与输出类型 | Subpopulation & output type

**通俗理解|In plain words:** `-s` 用来「把样本分成几组分别算、再画到一张图上对比」；`--out-type` 决定统计表里输出哪些 LD 度量(r² 还是再加 D')。不指定 `-s` 就对全体样本算一次。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-s, --subpop` | 无 | 子群体文件（样本 ID + 群体名，两列） |
| `-t, --type` | vcf | 输入类型：vcf / genotype |
| `--out-type` | 1 | 输出类型：1=r² / 2=r²+D' / 3=Pairwise LD |

### 绘图参数 | Plotting

**通俗理解|In plain words:** 这一组只影响图的平滑度和呈现方式，**不影响 `.stat.gz` 统计结果**。bin 是把相近距离的位点对归成一格再取代表值；`--method` 决定每格取什么代表值(平均/中位数/百分位)。**一般不用动。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--bin1` | 10 | 短距离 bin 大小(kb) |
| `--bin2` | 100 | 长距离 bin 大小(kb) |
| `--break-point` | 100 | 短/长距离分界点(kb) |
| `--max-x` | 无 | 图的最大 X 坐标(kb) |
| `--measure` | r2 | LD 度量方法：r2 / D / both |
| `--method` | MeanBin | 绘图统计方法：MeanBin / HW / MedianBin / PercentileBin |
| `--percentile` | 0.5 | 百分位数（仅 PercentileBin 用） |
| `--no-plot` | 关 | 不绘制图像 |
| `--no-recommend-threshold` | 关 | 不推荐 LD 阈值 |

## 分析流程 | Pipeline

```text
输入 VCF / Genotype
   |
   v
PopLDdecay 计算位点对 r² → 按距离分箱求平均(.stat.gz)
   |
   ├─ 绘图(.png，默认开，依赖 Perl 绘图脚本)
   ├─ LD 阈值推荐(.tsv，默认开，Hill & Weir 模型拟合)
   └─ 指定子群体时：先算全体(_all)，再逐个群体(_pop)，最后合并(_summary.tsv)
```

## 输出 | Output

```text
output.stat.gz                       # LD 衰减统计表(核心结果，gzip 压缩 TSV)
output.png                           # LD 衰减曲线图(默认生成)
output.log                           # 运行日志
output_threshold_recommendations.tsv # LD 阈值推荐(默认生成)
```

指定子群体(`-s`)时额外生成：

```text
output_all.stat.gz                   # 全体样本的 LD 衰减
output_{群体名}.stat.gz              # 每个子群体的 LD 衰减
output_summary.tsv                   # 合并表(Population/Dist/Mean_r2 列)
```

- `output.stat.gz`：gzip 压缩的 TSV，核心列 `#Dist`(距离 bp)、`Mean_r^2`(平均 r²)，可直接画图
- `output_threshold_recommendations.tsv`：每群体的样本数、背景 r²、推荐阈值、衰减距离、GWAS 建议窗口

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看 `output.png` 那条从左上到右下的曲线：起点越高=短距离内位点绑定越紧；降得越快=历史重组越多(群体越古老、有效群体越大)。

- **衰减距离(Decay_Distance_kb)**：r² 衰减到推荐阈值时对应的距离(kb)。越小=重组越频繁；GWAS 标记间距应小于衰减距离才能捕获信号
- **Recommended_Threshold**：推荐的 r² 阈值，r² 降到该值时位点间基本不再关联，是判断「标记是否够密」的参考
- **完整衰减表**：每个群体所有候选阈值的 `threshold / decay_kb / bg_ratio / recommendation` 会写入运行日志(.out 文件)；高于拟合曲线最大值的候选(如 r²=0.5——H&W 模型在 0 kb 处上限约 0.455，恒达不到)标为「not recommended: above fitted curve max」，不进入推荐
- **Background_r2**：远距离 r² 的「地板」，代表随机噪声水平，正常应接近 0
- **Rho**：Hill & Weir 模型拟合出的重组率参数，描述曲线下降的快慢
- **BG_Ratio**：推荐阈值除以背景 r² 的比值，接近 2 表示阈值选在「明显高于噪声」的合理位置
- **GWAS_Window_kb**：建议的关联分析窗口(lead SNP +/- X kb)，向上取整到 50 kb

## 参数选择建议 | Parameter Guidance

- `-d/--max-dist`：只算到多少 kb 内的位点对。默认 300 kb 覆盖绝大多数物种；衰减很慢的物种(如自交系)可调到 1000
- `-m/--min-maf`、`--max-het`、`--max-miss`：剔除低 MAF、高杂合、高缺失位点，**默认值一般不用动**
- `--bin1/--bin2/--break-point`：曲线分箱方式，**只影响图的平滑度，一般不用动**
- `--method`：绘图统计方法，默认 MeanBin(每箱取平均)；数据噪声大时可试 MedianBin
- `-s/--subpop`：想比较不同群体的衰减速度时用；不指定则对全体样本算一次
- `--no-recommend-threshold`：只想看曲线、不想要阈值推荐时加上
- `--threads`：**通俗理解|In plain words:** 并行计算用的「人手数」。人手按染色体分活——数据里有几条染色体/scaffold，最多就同时用几个人；给的数字超过染色体数也不报错，只是白给。人手多了求和顺序会变，结果可能在第 4 位小数上有 ±0.0001 的浮动（科学上无影响）；同样的 `--threads` 值重跑结果完全一致。**一般不用动**（默认 12），单染色体数据调了也没用

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 | Path | 输入VCF或Genotype文件｜Input VCF or Genotype file |
| `-o, --output` | 必填 | Path | 输出文件前缀｜Output file prefix |
| `-t, --type` | `vcf` | vcf/genotype | 输入文件类型｜Input file type (vcf/genotype) [default: vcf] |
| `-d, --max-dist` | `300` | int | 最大距离(kb)｜Max distance in kb [default: 300] |
| `-m, --min-maf` | `0.005` | float | 最小等位基因频率｜Min minor allele frequency [default: 0.005] |
| `--max-het` | `0.88` | float | 最大杂合率｜Max heterozygous rate [default: 0.88] |
| `--max-miss` | `0.25` | float | 最大缺失率｜Max missing rate [default: 0.25] |
| `-s, --subpop` | — | Path | 子群体样本列表文件｜Subpopulation sample list file |
| `--out-type` | `1` | int | 输出类型(1:R^2, 2:R^2&D', 3:Pairwise LD)｜Output type [default: 1] |
| `--bin1` | `10` | int | 短距离bin大小｜Bin size for short distance [default: 10] |
| `--bin2` | `100` | int | 长距离bin大小｜Bin size for long distance [default: 100] |
| `--break-point` | `100` | int | 短/长距离分界点｜Break point [default: 100] |
| `--max-x` | — | int | 最大X坐标(kb)｜Max X coordinate in kb |
| `--measure` | `r2` | r2/D/both | LD度量方法｜LD measure method [default: r2] |
| `--method` | `MeanBin` | MeanBin/HW/MedianBin/PercentileBin | 绘图方法｜Plotting method [default: MeanBin] |
| `--percentile` | `0.5` | float | 百分位数｜Percentile for PercentileBin [default: 0.5] |
| `--no-plot` | — |  | 不绘制图像｜Do not plot figure |
| `--no-recommend-threshold` | — |  | 不推荐LD阈值｜Do not recommend LD threshold |
| `--threads` | `12` | int | 计算线程数(按染色体分区,超过染色体数无效)｜Threads for LD calculation (partitioned per chromosome) [default: 12] |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- PopLDdecay2（默认 `~/miniforge3/envs/pop/bin/PopLDdecay2`，可用环境变量 `POPLDDECAY2_PATH` 覆盖；pop 环境含其依赖 htslib 1.24，二进制内嵌 RPATH 无需额外设置；经典版 PopLDdecay 仍保留在同目录可随时回退）
- Plot_OnePop.pl / Plot_MultiPop.pl（绘图 Perl 脚本，默认 `~/software/PopLDdecay2/bin/`，与经典版脚本逐字节相同）
- Perl（需 Data::Dumper、Getopt::Long 模块；绘图脚本自身还需 GD 等图形模块）
- Python 依赖 numpy / pandas / scipy（仅阈值推荐步骤）

## 常见问题 | FAQ

**Q1：`-o` 是前缀还是目录？**
是输出文件前缀。写成 `-o output` 生成 `output.stat.gz` 等；若 `-o` 以 `/` 结尾或指向已存在目录，程序会自动从输入文件名推导前缀(避免产出 `.stat.gz` 这类隐藏文件)。

**Q2：换参数重跑，旧结果还在吗？**
PopLDdecay 每次运行直接覆盖同名输出，不判断「已存在即跳过」。换了 `--max-dist`、`--bin1` 等参数重跑同一前缀即可，无需手动删旧文件；想保留多组结果请换不同前缀。

**Q3：子群体文件是什么格式？**
两列制表符分隔：`样本ID` 与 `群体名`（不是每行一个样本）。程序按群体名分组，先算全体(_all)、再逐个群体(_pop)、最后合并成 `_summary.tsv` 供对比。样本名须与 VCF 头部完全一致，否则被忽略。

**Q4：绘图失败怎么办？**
绘图依赖 Perl 的绘图脚本(Plot_OnePop.pl 等)及其 GD 模块；缺模块时程序会跳过绘图并给出警告，但 `.stat.gz` 统计文件仍会正常生成，可自己用 ggplot2 等工具画。

**Q5：阈值推荐是怎么算出来的？**
用 Hill & Weir (1988) 模型拟合衰减曲线，估计远距离的「背景 r²」噪声水平，再在候选阈值里挑一个约为背景 r² 两倍的 r² 作为推荐阈值，并据此给出衰减距离和 GWAS 建议窗口。高于拟合曲线最大值的候选会被标记并剔除——H&W 公式在 0 kb 处上限约 0.455(与 rho 无关)，所以候选 r²=0.5 永远达不到，背景 r² 偏高(>0.2)时它曾因 bg_ratio 最接近 2 被误选，输出「衰减距离 0.0 kb / GWAS 窗口 ±0」的无效结果；全部候选阈值的衰减表写入运行日志，可自行查看每个阈值对应的衰减距离。

**Q6：为什么我的数据加 `--threads` 没变快？**
多线程按染色体/scaffold 分区并行。数据只有 1 条染色体时只有一个「活」，加多少线程都是单线程在跑（属正常，不是故障）；真核基因组等多染色体数据才能吃满线程。另外单染色体数据可先 bgzip + tabix 建索引，PopLDdecay2 的快速通道会边读边算，仍有提速。

**Q7：能把别处编译好的 PopLDdecay2 二进制直接拷来超算用吗？**
不能。C++ 二进制绑定了编译机的 glibc 与动态库（如在 Mac/homebrew 新系统上编译的版本要求 GLIBC 2.33+，超算 RHEL8 只有 2.28，会报 `version GLIBC_2.33 not found`），且依赖 libhts/libdeflate/libcrypto。超算上的正确部署是本仓库 envs/pop.yml 的方式：pop 环境装 htslib 后用 `HTSLIB_ROOT=~/miniforge3/envs/pop bash make.sh` 本机编译，再用 patchelf `--force-rpath --set-rpath '$ORIGIN/../lib'` 嵌 RPATH（务必 DT_RPATH 而非默认 RUNPATH——登录环境的 LD_LIBRARY_PATH 可能指向其它软件带的旧 libhts，RUNPATH 会被它压过导致符号版本冲突）。

**Q8：子群体最少要几个样本？**
PopLDdecay 工具自身要求子群体 ≥3 个样本，少于 3 个会报 `sub Group Population szie is too small`（新旧版本行为一致）；模块会记 WARNING 并跳过该群体，其余群体与全体(_all)结果照常产出。