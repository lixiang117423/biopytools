# ADMIXTURE 群体结构分析 | ADMIXTURE Population Structure Analysis

一句话理解：**估计每个个体来自 K 个「祖先群体」的比例**，从而看清样本的群体结构、混合(admixture)情况和可能的谱系来源。

## 功能概述 | Overview

- 封装 ADMIXTURE 或 ADAMIXTURE 两种方法(最大似然估计祖先成分)
- 自动跑 K=min_k 到 max_k 每个值，用交叉验证(CV)误差选出最优 K
- 自动 VCF 预处理：双等位+SNP 过滤、质控(MAF/位点缺失,可选HWE)、LD 剪枝；**只做位点层面过滤,保留全部样品**(不删个体)
- 输出每个个体的祖先成分(Q)、各祖先等位基因频率(P)、CV 误差、最优 K
- 顺带生成 GWAS 协变量文件、可视化堆叠图与总结报告
- 支持断点续传与 `--force` 强制重算

## 快速开始 | Quick Start

```bash
biopytools admixture -i input.vcf -o admixture_results
```

默认跑 K=2 到 10，用交叉验证选最优 K。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| K 值 | 假设的「祖先群体」个数；K=3 就是假设样本由 3 个祖先群体混合而来 |
| 祖先成分(Q) | 每个个体的基因组按比例「拆」到各祖先群体的占比，加起来=1 |
| 交叉验证(CV) | 把数据切几份轮流留一份验证，误差越小说明该 K 拟合越好 |
| 群体结构 | 样本按遗传背景自然分层/分群，常对应地理或谱系 |
| 混合(admixture) | 个体的祖先是来自不同群体的「混血」，一个个体可同时属多个成分 |
| LD 剪枝 | 去掉「绑定遗传」的冗余位点，避免它们重复计票(ADMIXTURE 假设位点独立) |
| MAF | 「少数派」等位基因频率；太低=位点无区分度 |

## 输入 | Input

标准 VCF(支持 .vcf/.vcf.gz)，含基因型(GT)。建议至少含十几个样本、数千个 SNP。

## 分析流程 | Pipeline

```text
VCF
  │ 步骤1: 预处理(bcftools 只留双等位 SNP → biallelic.vcf.gz)
  ▼
步骤2: 转 PLINK 格式(.bed/.bim/.fam)
  │
  ▼
步骤3: 质控(MAF/位点缺失;HWE默认关闭,可 --skip-preprocessing 跳过)
  │
  ▼
LD 剪枝(默认开) → 修复染色体编号
  │
  ▼
步骤4: ADMIXTURE 逐 K 计算(交叉验证 → cv_results.csv)
  │
  ▼
步骤5: 选最优 K → 祖先成分表(admixture_proportions.csv)
  │
  ├─ 步骤6: 生成 GWAS 协变量文件
  └─ 可视化堆叠图 + 总结报告
```

## 输出 | Output

```text
admixture_results/
├── 00_pipeline_info/software_versions.yml   # 软件版本记录
├── 01_preprocessing/
│   ├── biallelic.vcf.gz                     # 双等位 SNP 过滤后的 VCF
│   └── chromosome_mapping.txt               # 染色体重命名对应表
├── 02_plink/                                # PLINK 中间文件(.bed/.bim/.fam)
├── 03_admixture/
│   ├── cv_results.csv                       # 各 K 的交叉验证误差(选最优 K 用)
│   ├── loglikelihood_results.csv            # 各 K 的对数似然(ADAMIXTURE)
│   └── <前缀>.<K>.Q / <前缀>.<K>.P          # 各 K 的祖先成分 / 等位基因频率
├── 04_results/
│   ├── admixture_proportions.csv            # 每个个体的祖先成分(核心结果)
│   ├── admixture_statistics.txt             # 混合程度统计
│   ├── gwas_covariates.txt                  # GWAS 协变量文件
│   ├── analysis_summary.txt                 # 分析总结报告
│   └── *.pdf                                # 个体堆叠图/可视化
└── 99_logs/admixture_analysis.log           # 运行日志
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看 `cv_results.csv` 找 CV 误差最小的 K(那就是最优 K)；再看该 K 的 `admixture_proportions.csv`，每行一个个体、各祖先成分占比加起来=1。

- **cv_results.csv**：每行一个 K 值对应一个 CV 误差，误差越小拟合越好；取误差最小(或开始「平台」)的 K 作为最优
- **admixture_proportions.csv**：列为各祖先成分(K1..Kn)，某个体某成分接近 1 说明它基本「纯」，多个成分都有值说明它是「混血」
- **admixture_statistics.txt**：高度混合个体(max 成分 < 0.7)、纯合个体(max 成分 > 0.9)的数量，辅助判断混合程度
- **堆叠图(*.pdf)**：每个个体一根竖条，用颜色段表示各祖先成分比例，一眼看出分群
- **gwas_covariates.txt**：把祖先成分转成协变量，可直接喂给 GWAS 软件控制群体分层

## 参数选择建议 | Parameter Guidance

- **K 范围**：`-k`/`-K` 决定假设祖先数范围；先跑宽范围(如 2-10)，看 CV 曲线再收窄
- **`--method`**：默认 `admixture`；数据量大/想更快可试 `adamixture`(随机梯度下降版，参数 `--adamixture-*` 一般不用动)
- **质控阈值**：`--maf 0.05`/`--missing 0.1`(位点缺失率,`--geno`) **一般不用动**；样品少位点丢得多可放宽 `--missing 0.2`。**HWE 默认关闭**——混合群体偏离 HWE 属预期,过滤反而会误删群体分化位点;确需启用传 `-H 1e-6`。个体缺失率(`--mind`)已移除,保证不删样品；要跳过全部质控用 `--skip-preprocessing`
- **LD 剪枝**：默认开启(ADMIXTURE 假设位点独立)，`--ld-window 50`/`--ld-step 10`/`--ld-r2 0.1`(经典保守组合) **一般不用动**；过滤后 SNP 偏多(>5万)可收紧 r²,SNP 偏少(<1000)可放宽到 `100 10 0.2`
- **`--keep-intermediate`**：加它保留中间 PLINK/VCF 文件，排查问题或复用


<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -i` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `admixture_results` |  | 输出目录｜Output directory |
| `--method` | `admixture` | admixture/adamixture | 分析方法｜Analysis method (admixture or adamixture) |
| `--min-k, -k` | `2` | int | 最小K值｜Minimum K value |
| `--max-k, -K` | `10` | int | 最大K值｜Maximum K value |
| `--cv-folds, -c` | `5` | int | 交叉验证折数｜Cross-validation folds (仅ADMIXTURE｜ADMIXTURE only) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--maf, -m` | `0.05` | float | 最小等位基因频率阈值｜MAF threshold |
| `--missing, -M` | `0.1` | float | 位点缺失率阈值(--geno,不删个体)｜Site-level missing threshold (--geno, never drops samples) |
| `--hwe, -H` | `1.0` | float | HWE p值阈值,1=关闭(混合群体默认不过滤)｜HWE p-value threshold, 1=disabled (default for admixed panels) |
| `--ld-prune/--no-ld-prune` | `True` |  | LD剪枝(默认开启,ADMIXTURE假设位点独立)｜LD pruning (on by default; ADMIXTURE assumes unlinked sites) |
| `--ld-window` | `50` |  | LD剪枝窗口(kb或SNP数)｜LD pruning window (kb or SNP count) |
| `--ld-step` | `10` | int | LD剪枝步长｜LD pruning step size |
| `--ld-r2` | `0.1` | float | LD剪枝r2阈值｜LD pruning r2 threshold |
| `--skip-preprocessing, -s` | — |  | 跳过VCF预处理｜Skip VCF preprocessing |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--verbose, -v` | — |  | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--force, -f` | — |  | 强制覆盖已存在的文件｜Force overwrite existing files |
| `--dry-run` | — |  | 试运行模式(不实际执行)｜Dry run without execution |
| `--adamixture-path` | `~/miniforge3/envs/adamixture_v.1.0.2/bin/adamixture` |  | ADAMIXTURE可执行文件路径｜ADAMIXTURE executable path |
| `--adamixture-lr` | `0.005` | float | ADAMIXTURE学习率｜ADAMIXTURE learning rate |
| `--adamixture-beta1` | `0.8` | float | ADAMIXTURE beta1参数｜ADAMIXTURE beta1 parameter |
| `--adamixture-beta2` | `0.88` | float | ADAMIXTURE beta2参数｜ADAMIXTURE beta2 parameter |
| `--adamixture-max-iter` | `1500` | int | ADAMIXTURE最大迭代次数｜ADAMIXTURE maximum iterations |
| `--adamixture-seed` | `42` | int | ADAMIXTURE随机种子｜ADAMIXTURE random seed |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | `admixture_results` |  | 输出目录｜Output directory |
| `--method` | `admixture` | admixture/adamixture | 分析方法｜Analysis method (admixture or adamixture) |
| `-k, --min-k` | `2` | int | 最小K值｜Minimum K value |
| `-K, --max-k` | `10` | int | 最大K值｜Maximum K value |
| `-c, --cv-folds` | `5` | int | 交叉验证折数｜Cross-validation folds |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--adamixture-path` | `~/miniforge3/envs/adamixture_v.1.0.2/bin/adamixture` |  | ADAMIXTURE可执行文件路径｜ADAMIXTURE executable path |
| `--adamixture-lr` | `0.005` | float | ADAMIXTURE学习率｜ADAMIXTURE learning rate |
| `--adamixture-beta1` | `0.8` | float | ADAMIXTURE beta1参数｜ADAMIXTURE beta1 parameter |
| `--adamixture-beta2` | `0.88` | float | ADAMIXTURE beta2参数｜ADAMIXTURE beta2 parameter |
| `--adamixture-max-iter` | `1500` | int | ADAMIXTURE最大迭代次数｜ADAMIXTURE maximum iterations |
| `--adamixture-seed` | `42` | int | ADAMIXTURE随机种子｜ADAMIXTURE random seed |
| `-m, --maf` | `0.05` | float | MAF阈值｜MAF threshold |
| `-M, --missing` | `0.1` | float | 位点缺失率阈值(--geno,不删个体)｜Site-level missing threshold (--geno, never drops samples) |
| `-H, --hwe` | `1.0` | float | HWE p值阈值,1=关闭(混合群体默认不过滤)｜HWE p-value threshold, 1=disabled (default for admixed panels) |
| `--no-ld-prune` | — | store_true | 关闭LD剪枝(默认开启)｜Disable LD pruning (on by default) |
| `--ld-window` | `50` |  | LD剪枝窗口(kb或SNP数)｜LD pruning window (kb or SNP count) |
| `--ld-step` | `10` | int | LD剪枝步长｜LD pruning step size |
| `--ld-r2` | `0.1` | float | LD剪枝r2阈值｜LD pruning r2 threshold |
| `-s, --skip-preprocessing` | — | store_true | 跳过VCF预处理和质控｜Skip VCF preprocessing and QC |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level (default: INFO) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run without execution |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- ADMIXTURE(默认 `~/miniforge3/envs/pop/bin/admixture`；ADAMIXTURE 默认 `~/miniforge3/envs/adamixture_v.1.0.2/bin/adamixture`)
- PLINK(默认 `~/miniforge3/envs/pop/bin/plink`，用于转格式/质控/LD 剪枝)
- bcftools(默认 `~/miniforge3/envs/align/bin/bcftools`，用于双等位过滤)
- R(用于生成可视化堆叠图)

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
本模块有断点续传：按输出文件是否已存在判断步骤是否完成。换质控参数(如 `--maf`)或 LD 参数重跑旧目录前，先删对应旧产物(如 `01_preprocessing/biallelic.vcf.gz`、`02_plink/*.bed`、`03_admixture/cv_results.csv`)，或直接加 `--force` 强制重算。

**Q2：最优 K 怎么选？**
看 `cv_results.csv`：CV 误差最小(或不再明显下降、进入「平台」)的 K 即最优。但也要结合生物学背景——有时 CV 最小的 K 不一定最有解释力。

**Q3：日志提示「未确定最优 K 值」？**
说明各 K 的交叉验证没有成功产出结果(常见于位点过少、样本过少或 ADMIXTURE 运行失败)。检查 `03_admixture/` 下各 K 的日志(`log_k.out`)，并核对输入 VCF 的 SNP 数量与样本数。

**Q4：ADMIXTURE 和 ADAMIXTURE 有什么区别？**
ADMIXTURE 是经典的最大似然实现，结果稳定；ADAMIXTURE 用随机梯度下降，速度快、适合超大 SNP 集，代价是需要调学习率等参数(默认值一般可用)。两者结果应大致一致。

**Q5：K 值必须大于样本数吗？**
K 不能超过样本数，否则 ADMIXTURE 会失败。设置 `-K` 时别超过 VCF 里的样本数量。