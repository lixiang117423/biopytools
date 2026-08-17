# K-mer GWAS 分析 | K-mer-based GWAS

一句话理解：**不用传统 SNP、而是把基因组切成一个个固定长度的小片段(k-mer)来做关联分析，专门解决「没有参考基因组」和「多倍体」物种找不到性状关联位点的难题**。输入一堆重测序 FASTQ + 一个表型表，输出「哪些 k-mer 与性状显著相关」及它们落在基因组的什么位置。

## 功能概述 | Overview

- 完整 k-mer GWAS 流程：`count`(计数)→ `kctm`(建矩阵)→ `filter`(过滤)→ `m2b`(转 BIMBAM)→ `asso`(关联分析)，可一键 `pipeline` 跑完，也可分步跑
- 基于 KMERIA 软件，对无参考基因组、多倍体物种友好
- Post-GWAS：把显著 k-mer 用 bwa/blast 比对回参考基因组定位候选基因(需参考基因组)
- 断点续传：产出存在即跳过(注意 `--force` 只在 pipeline/count 有，见 FAQ)
- 自动生成质控报告和曼哈顿图/QQ图数据

## 快速开始 | Quick Start

```bash
biopytools kmeria pipeline -i fastq_dir --samples samples.txt -d depth.txt -p pheno.txt -o results
```

最小输入：FASTQ 目录 + 样本列表 + 测序深度文件 + 表型文件。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| k-mer | 把 DNA 序列切成固定长度 k 的「小碎片」(如 k=31 就是 31 个碱基一段)，像把长句子切成固定字数的词组 |
| k-mer GWAS | 不做 SNP 分型、直接拿「每个样本里有哪些 k-mer、有多少」做关联，适合无参考/多倍体 |
| 丰度(abundance) | 一个 k-mer 在样本里出现多少次；太低=测序噪声，太高(>1000)=重复序列污染 |
| 缺失率(missing ratio) | 一个 k-mer 在多少样本里「没出现」；太高=不稳定，过滤掉 |
| 倍性(ploidy) | 基因组有几套染色体；二倍体=2、四倍体=4，影响丰度阈值判断 |
| 测序深度(depth) | 每个位置平均被测到几遍，用于把丰度折算成「相对拷贝数」 |
| BIMBAM | 关联分析软件的经典输入格式(基因型矩阵)，m2b 步骤就是把 k-mer 矩阵转成它 |
| Post-GWAS | 把显著 k-mer「翻译」回基因组坐标并找附近基因 |

## 输入 | Input

### 1. FASTQ 目录

双端测序命名自动配对：`{sample}_R1.fq.gz` + `{sample}_R2.fq.gz`，或 `{sample}_1/`_2`。支持 `.fq.gz`/`.fastq.gz`/`.fq`/`.fastq`。

### 2. 样本列表(`--samples`)

纯文本，每行一个样本名，无表头，忽略空行和 `#` 开头的行：

```text
sample1
sample2
sample3
```

### 3. 测序深度文件(`-d/--depth-file`)

Tab 分隔，无表头：`样本名  深度`(第 3 列可选，给「混合倍性」样本单独指定倍性)：

```text
sample1    45.2
sample2    52.8
sample3    38.9
```

### 4. 表型文件(`-p/--pheno-file`)

空白或 Tab 分隔，无表头：`样本名  表型值  [协变量…]`。默认用第 2 列(第 1 列是样本名)作表型，用 `--pheno-col` 改(1-based)：

```text
sample1    1.5
sample2    2.3
sample3    1.8
```

## 参数说明 | Parameters

### pipeline - 一键流程 | Pipeline

**通俗理解|In plain words:** 一个命令跑完五步+可选 Post-GWAS。`-k` 是 k-mer 长度(默认 31 基本不动)；`--max-abund` 砍掉重复序列污染；`--missing-ratio` 过滤「在很多样本里缺席」的 k-mer；`--ploidy` 按物种倍性设。**注意：`biopytools kmeria` 入口的默认值是 missing-ratio 0.05、ploidy 2、threads 12(与 `python -m biopytools.kmeria` 直调的 0.8/4/24 不同)。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i, --fastq-dir` | 必填 | FASTQ 目录 |
| `--samples` | 必填 | 样本列表文件 |
| `-d, --depth-file` | 必填 | 测序深度文件 |
| `-p, --pheno-file` | 必填 | 表型文件 |
| `-k, --kmer-size` | `31` | k-mer 大小 |
| `--max-abund` | `1000` | 最大丰度(去重复污染) |
| `--missing-ratio` | `0.05` | 缺失率阈值 |
| `--ploidy` | `2` | 倍性 |
| `--step` | 全跑 | 从指定步骤开始(count/kctm/filter/m2b/asso) |
| `-f, --force` | 关 | 强制重跑所有步骤 |

### Post-GWAS 注释(可选) | Post-GWAS annotation

**通俗理解|In plain words:** 拿到显著 k-mer 后想「翻译成基因」，就给参考基因组。`--genome-file` + `--gff-file` 开启注释，`--alignment-tool` 选 bwa(快)或 blast(慢但更准)。**不给基因组文件就跳过这步，不影响关联分析本身。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--genome-file` | 无 | 参考基因组(注释用) |
| `--gff-file` | 无 | GFF 注释文件 |
| `--alignment-tool` | `bwa` | bwa/blast |
| `--window-size` | `200000` | 基因查找窗口(bp) |

### 分步子命令 | Step subcommands

**通俗理解|In plain words:** 五步可以拆开单独跑(每步的输入是上一步的输出)。各步参数与 pipeline 里同名参数一致，一般只在「只想重跑某一步」时用。

| 子命令 | 作用 | 输入 → 输出 |
|--------|------|-------------|
| `count` | k-mer 计数 | FASTQ → 01_kmer_counts |
| `kctm` | 建 k-mer 矩阵 | 计数目录 → 02_kmer_matrices |
| `filter` | 过滤(丰度/缺失/倍性) | 矩阵 → 03_filtered_matrices |
| `m2b` | 转 BIMBAM 格式 | 过滤目录 → 04_bimbam |
| `asso` | 关联分析 | BIMBAM + 表型 → 05_association |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先把每个样本的 k-mer 数一遍 → 汇总成一张「样本 × k-mer」大表 → 按丰度/缺失率过滤掉噪声 → 转成关联软件认识的格式 → 逐 k-mer 算与性状的相关性 → (可选)把显著 k-mer 比对回基因组找基因。

```text
FASTQ 目录
    │ count:  每个样本数 k-mer → {sample}_k31.bin
    ▼
kctm:  合并成 k-mer 矩阵 → kmer_matrix.*.bin
    ▼
filter: 按 max-abund / missing-ratio / ploidy 过滤
    ▼
m2b:  转成 BIMBAM 格式 → *.bimbam.gz
    ▼
asso:  逐 k-mer 关联分析 → *.ps(结果)
    ▼
(可选) Post-GWAS: 显著 k-mer 比对回基因组 → 定位候选基因
```

## 输出 | Output

```text
results/
├── 01_kmer_counts/            # 每样本 {sample}_k31.bin
├── 02_kmer_matrices/          # kmer_files.txt + kmer_matrix.*.bin
├── 03_filtered_matrices/      # filtered_*.bin / .txt
├── 04_bimbam/                 # *.bimbam.gz
├── 05_association/            # *.ps(关联结果) + samples.txt
├── 06_qc_reports/             # qc_report.json
├── 07_post_gwas_bwa/          # Post-GWAS 结果(bwa 比对, 显著 k-mer 定位)
│   └── (或 07_post_gwas_blast/)
└── kmeria.log                 # 日志
```

## 结果解读 | Interpreting Results

### 关联结果(`05_association/*.ps`)

**通俗理解|In plain words:** 核心结果，一行一个 k-mer，含 P 值。P 越小越显著；画成曼哈顿图找「冒尖的塔」。显著 k-mer(如 P < 1e-5)就是候选的「性状相关片段」。

### Post-GWAS 定位(`07_post_gwas_*/`)

**通俗理解|In plain words:** 把显著 k-mer 比对回参考基因组，看它们落在哪些基因附近——这是「从 k-mer 到候选基因」的最后一步。结果里的基因组坐标+附近基因就是候选名单。

### 质控报告(`06_qc_reports/qc_report.json`)

**通俗理解|In plain words:** 各步骤的统计(每个样本多少 k-mer、过滤掉了多少)。过滤比例异常高说明数据质量或阈值有问题。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规分析 | 全部默认，`pipeline` 一键跑 |
| 二倍体 | `--ploidy 2`(默认) |
| 四倍体/六倍体 | `--ploidy 4` / `--ploidy 6` |
| 重复序列多 | 调小 `--max-abund`(如 500) |
| 想做基因注释 | 加 `--genome-file` + `--gff-file`，`--alignment-tool bwa` |
| 无参考基因组 | 不给 `--genome-file`，只看关联结果本身 |
| 某步失败续跑 | 用 `--step` 从失败步骤开始 |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --fastq-dir` | 必填 |  | FASTQ文件目录｜FASTQ files directory |
| `--sample` | 必填 |  | 样本列表文件｜Sample list file |
| `-d, --depth-file` | 必填 |  | 测序深度文件｜Sequencing depth file |
| `-p, --pheno-file` | 必填 |  | 表型文件｜Phenotype file |
| `-o, --output-dir` | `./kmeria_results` |  | 输出目录｜Output directory |
| `-f, --force` | `False` |  | 强制重新运行所有步骤｜Force re-run all steps even if output exists |
| `-k, --kmer-size` | `31` |  | K-mer大小｜K-mer size (default: 31) |
| `--max-abund` | `1000` |  | 最大丰度｜Maximum abundance (default: 1000) |
| `--missing-ratio` | `0.05` |  | 缺失率｜Missing ratio (default: 0.05) |
| `--ploidy` | `2` |  | 倍性｜Ploidy (default: 2) |
| `--step` | — | count/kctm/filter/m2b/asso | 从指定步骤开始｜Start from specified step |
| `-t, --threads` | `12` |  | 线程数｜Thread count (default: 12) |
| `--batch-size` | `4` |  | 批处理大小｜Batch size (default: 4) |
| `--pheno-col` | `1` |  | 表型列｜Phenotype column (default: 1) |
| `--kinship-file` | — |  | 亲缘关系矩阵｜Kinship matrix file |
| `--covar-file` | — |  | 协变量文件｜Covariate file |
| `--enable-qc` | `True` |  | 启用质控｜Enable QC (default: True) |
| `--enable-visualization` | `True` |  | 启用可视化｜Enable visualization (default: True) |
| `--enable-annotation` | — |  | 启用k-mer注释｜Enable k-mer annotation |
| `--genome-file` | — |  | 参考基因组｜Reference genome (for annotation) |
| `--gff-file` | — |  | GFF注释文件｜GFF annotation file (for annotation) |
| `--sample-ratio` | `0.1` | float | 高p值位点抽样比例，用于减少绘图点数 (默认: 0.1 = 10%)｜Sampling ratio for high p-value loci to reduce plot points (default: 0.1 = 10%) |
| `--window-size` | `200000` | int | 基因查找窗口大小，单位bp (默认: 200000 = 200kb)｜Gene search window size in bp (default: 200000 = 200kb) |
| `--alignment-tool` | `bwa` | bwa/blast | Post-GWAS比对工具选择 (默认: bwa)｜Alignment tool for Post-GWAS (default: bwa) |
| `--bwa-k` | `9` | int | BWA mem -k 参数，最小种子长度 (默认: 9)｜BWA mem -k parameter, minimum seed length (default: 9) |
| `--bwa-t-min-score` | `10` | int | BWA mem -T 参数，最小输出分数 (默认: 10)｜BWA mem -T parameter, minimum score to output (default: 10) |
| `--as-ratio` | `0.95` | float | BWA AS过滤阈值，保留AS >= 最高AS * as_ratio的所有比对 (默认: 0.95)｜BWA AS filtering threshold, keep alignments with AS >= max_AS * ratio (default: 0.95) |
| `--log-file` | — |  | 日志文件｜Log file |
| `-b, --batch-size` | `4` |  | 批处理大小｜Batch size (default: 4) |
| `-C, --count-separate-strands` | — |  | 分别计数链｜Count strands separately |
| `-T, --text-output` | — |  | 文本输出｜Text output |
| `-i, --input-dir` | 必填 |  | 输入目录｜Input directory |
| `-c, --max-abund` | `1000` |  | 最大丰度｜Maximum abundance (default: 1000) |
| `-s, --missing-ratio` | `0.05` |  | 缺失率｜Missing ratio (default: 0.05) |
| `-p, --ploidy` | `2` |  | 倍性｜Ploidy (default: 2) |
| `--in, -i` | 必填 |  | 输入目录｜Input directory |
| `--out, -o` | `./04_bimbam` |  | 输出目录｜Output directory |
| `--no-normalize` | — |  | 不归一化｜No normalization |
| `--quantile-norm` | — |  | 分位数归一化｜Quantile normalization |
| `-n, --pheno-col` | `1` |  | 表型列｜Phenotype column (default: 1) |
| `-c, --covar-file` | — |  | 协变量文件｜Covariate file |
| `-k, --kinship-file` | — |  | 亲缘关系矩阵｜Kinship matrix file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --fastq-dir` | 必填 |  | FASTQ文件目录｜FASTQ files directory |
| `--sample, --samples` | 必填 |  | 样本列表文件｜Sample list file |
| `-d, --depth-file` | 必填 |  | 测序深度文件｜Sequencing depth file |
| `-p, --pheno-file` | 必填 |  | 表型文件｜Phenotype file |
| `-o, --output-dir` | `./kmeria_results` |  | 输出目录｜Output directory |
| `-k, --kmer-size` | `31` | int | K-mer大小｜K-mer size (default: 31) |
| `--max-abund` | `1000` | int | 最大丰度｜Maximum abundance (default: 1000) |
| `--missing-ratio` | `0.8` | float | 缺失率（k-mer默认0.8，与GitHub一致）｜Missing ratio (k-mer default 0.8, matches GitHub) |
| `--ploidy` | `4` | int | 倍性（默认4，适配多倍体）｜Ploidy (default 4, for polyploids) |
| `--step` | — | count/kctm/filter/m2b/asso | 从指定步骤开始｜Start from specified step |
| `-t, --threads` | `24` | int | 线程数｜Thread count (default: 24) |
| `--batch-size` | `4` | int | 批处理大小｜Batch size (default: 4) |
| `--pheno-col` | `1` | int | 表型列｜Phenotype column (default: 1) |
| `--kinship-file` | — |  | 亲缘关系矩阵｜Kinship matrix file |
| `--covar-file` | — |  | 协变量文件｜Covariate file |
| `--enable-qc` | `True` | store_true | 启用质控｜Enable QC (default: True) |
| `--genome-file` | — |  | 参考基因组｜Reference genome (for Post-GWAS analysis) |
| `--gff-file` | — |  | GFF注释文件｜GFF annotation file (optional) |
| `--sample-ratio` | `0.1` | float | 高p值位点抽样比例｜Sampling ratio for high p-value loci (default: 0.1) |
| `--window-size` | `200000` | int | 基因查找窗口大小｜Gene search window size (default: 200000 = 200kb) |
| `--alignment-tool` | `bwa` | bwa/blast | Post-GWAS比对工具选择 (默认: bwa)｜Alignment tool for Post-GWAS (default: bwa) |
| `--bwa-k` | `9` | int | BWA mem -k 参数，最小种子长度 (默认: 9)｜BWA mem -k parameter, minimum seed length (default: 9) |
| `--bwa-t-min-score` | `10` | int | BWA mem -T 参数，最小输出分数 (默认: 10)｜BWA mem -T parameter, minimum score to output (default: 10) |
| `--as-ratio` | `0.95` | float | BWA AS过滤阈值，保留AS >= 最高AS * ratio的所有比对 (默认: 0.95)｜BWA AS filtering threshold, keep alignments with AS >= max_AS * ratio (default: 0.95) |
| `--log-file` | — |  | 日志文件｜Log file |
| `-f, --force` | — | store_true | 强制重新运行所有步骤｜Force re-run all steps |
| `-b, --batch-size` | `4` | int | 批处理大小｜Batch size (default: 4) |
| `-C, --count-separate-strands` | — | store_true | 分别计数链｜Count strands separately |
| `-T, --text-output` | — | store_true | 文本输出｜Text output |
| `-i, --input-dir` | 必填 |  | 输入目录｜Input directory |
| `-c, --max-abund` | `1000` | int | 最大丰度｜Maximum abundance (default: 1000) |
| `-s, --missing-ratio` | `0.8` | float | 缺失率（默认0.8，与GitHub一致）｜Missing ratio (default 0.8, matches GitHub) |
| `-p, --ploidy` | `4` | int | 倍性（默认4）｜Ploidy (default 4) |
| `--in, -i` | 必填 |  | 输入目录｜Input directory |
| `--out, -o` | `./04_bimbam` |  | 输出目录｜Output directory |
| `--no-normalize` | — | store_true | 不归一化｜No normalization |
| `--quantile-norm` | — | store_true | 分位数归一化｜Quantile normalization |
| `-n, --pheno-col` | `1` | int | 表型列｜Phenotype column (default: 1) |
| `-c, --covar-file` | — |  | 协变量文件｜Covariate file |
| `-k, --kinship-file` | — |  | 亲缘关系矩阵｜Kinship matrix file |
| `--use-gemma` | — | store_true | 使用gemma替代bimbamAsso｜Use gemma instead of bimbamAsso |
| `-m, --analysis-method` | — |  | 分析方法（default/lm/lmm）｜Analysis method |
| `--kin-method` | `3` | int | kinship计算方法（1=IBS均值, 2=IBS随机, 3=Balding-Nichols）｜Kinship method (default: 3) |
| `--no-kinship` | — | store_true | 不使用kinship矩阵｜Do not use kinship matrix |
| `--disable-gls` | — | store_true | 禁用GLS，使用OLS｜Disable GLS, use OLS |
| `--write-eigen` | — | store_true | 输出特征值/特征向量｜Write eigenvalues/eigenvectors |
| `--kin-precision` | `10` | int | kinship精度｜Kinship precision (default: 10) |
| `--out-precision` | `5` | int | 输出精度｜Output precision (default: 5) |
| `--maf` | — | float | 次等位基因频率过滤｜Minor allele frequency filter |
| `--miss` | — | float | 缺失阈值｜Missing threshold |
| `--start-marker` | — | int | 起始marker索引｜Start marker index |
| `--end-marker` | — | int | 结束marker索引｜End marker index |
| `--generate-plots` | — | store_true | 生成图表｜Generate plots |
| `--compress` | — | store_true | 压缩输出｜Compress output |
| `--verbose` | — | store_true | 详细输出｜Verbose output |
| `--dry-run` | — | store_true | 仅显示命令不执行｜Show commands without executing |
| `--no-validate` | — | store_true | 跳过输入验证｜Skip input validation |
| `--no-check-deps` | — | store_true | 跳过依赖检查｜Skip dependency check |
| `--no-cleanup` | — | store_true | 保留临时文件｜Keep temporary files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- KMERIA(conda 环境 `kmeriaenv`，二进制 `~/software/kmeria/bin/kmeria`，自动 `conda run` 包装)
- bwa(Post-GWAS 比对，**需在 PATH**，裸调用不经 conda)、blast(makeblastdb/blastn，Post-GWAS 备选)
- Rscript(画图)、samtools(算深度时用)
- gemma 可选(仅 `--use-gemma` 时，默认用 bimbamAsso)

## 常见问题 | FAQ

**Q1：为什么 `biopytools kmeria` 和 `python -m biopytools.kmeria` 的默认值不一样？**
click 包装器会把自己的默认值(threads=12、missing-ratio=0.05、ploidy=2)显式传给底层 argparse，而直调入口用的是代码默认(24 线程、0.8、4)。**用 `biopytools kmeria` 入口时按 12/0.05/2 记**。

**Q2：想强制重跑某一步怎么办？**
`--force` 只在 `pipeline` 和 `count` 暴露；kctm/filter/m2b/asso 没有 `--force`，产物存在就永远跳过。要重跑请先删掉对应步骤的输出目录。

**Q3：换参数(如 --kmer-size)重跑，结果没变？**
断点续传按产物存在性判断。换 `-k` 等参数后需删除旧产物(如 `01_kmer_counts`)再重跑，否则复用旧结果。

**Q4：需要参考基因组吗？**
关联分析本身不需要；只有想「把显著 k-mer 翻译成基因」时才需要 `--genome-file`。

**Q5：k-mer 大小怎么选？**
默认 31 适合多数场景。较小(15-21)更敏感但假阳性高；较大(25-31)更特异但费资源。
