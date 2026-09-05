
# Filter SNP/INDEL - VCF 变异过滤 | SNP/INDEL Variant Filtering

一句话理解：**给已 call 出的 VCF 做一次「体检打分」，按 GATK 最佳实践把 SNP 和 INDEL 分开、各用一套标准过滤掉低质量变异**，得到一份可信的变异清单。

## 功能概述 | Overview { #overview }

- 自动分离 SNP 与 INDEL（bcftools view --types），分别应用不同的质控标准
- 按 GATK 最佳实践过滤：QUAL/DP/MQ/QD/FS/SOR/MQRankSum/ReadPosRankSum/MAF 等
- SNP 额外支持「只保留双等位位点」与 MAF 过滤
- 可选 --repair-vcf 自动修复列数不匹配等损坏的 VCF
- 合并过滤后的 SNP+INDEL，输出规整的最终 VCF 与统计报告

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools filter-snp-indel -i variants.vcf -o filtered_output/
```

最小输入：一个 VCF（.vcf / .vcf.gz）。输出目录默认 ./filtered_vcf。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| SNP | 单个碱基的改变（A 变 G 之类） |
| INDEL | 插入或删除一小段碱基（InDels） |
| QUAL | 变异整体的质量分，越高越可信 |
| DP(深度) | 这个位点被 reads 覆盖了多少层 |
| MQ(比对质量) | reads 贴到这个位置的信心 |
| QD(质量/深度比) | 把质量分摊到每层深度，排除「高深度堆出来的假质量」 |
| FS(链偏倚) | 正负链证据是否一边倒，偏得太厉害可能是假变异 |
| SOR(链比值) | 类似 FS，衡量链偏倚的另一种打分 |
| MQRankSum/ReadPosRankSum | 衡量支持变异的 reads 是否「质量偏低/位置集中」 |
| MAF | 少数等位基因频率，太低=罕见到可疑 |
| 双等位位点 | 只有两种碱基可能性的位点（更可靠） |

## 输入 | Input { #input }

- **VCF**（`-i`）：标准 VCF，支持 .vcf / .vcf.gz，可含 SNP 和 INDEL 混合

## 参数说明 | Parameters { #parameters }

### 必需与基础 | Required & basic

**通俗理解|In plain words:** -i 是输入 VCF，-o 输出目录，-t 线程数（默认 12）。--variant-type 告诉程序输入里是什么（both/snp_only/indel_only），默认 both 自动分离。

### SNP 过滤 | SNP filtering

**通俗理解|In plain words:** 这组是 SNP 的「及格线」：QUAL>=30、DP>=10、MQ>=40、QD>=2、FS<=60、SOR<=3、MQRankSum>=-12.5、ReadPosRankSum>=-8、MAF>=0.05。调严=删更多、更保守；调松=保留更多。**默认是文献常用 GATK 标准，一般不用动。**

### INDEL 过滤 | INDEL filtering

**通俗理解|In plain words:** INDEL 用更宽松的标准（FS<=200、SOR<=10），因为 INDEL 本身链偏倚天然更高。**默认即合理，一般不用动。**

### 双等位位点过滤 | Biallelic filtering

**通俗理解|In plain words:** 默认 --snp-biallelic 开启，只保留「只有两种等位基因」的 SNP（更可靠）；多等位位点想保留时用 --no-snp-biallelic 关闭。

### 执行与日志 | Execution & logging

**通俗理解|In plain words:** --repair-vcf 自动修复列数不匹配等损坏；--force 覆盖已存在结果；--dry-run 只打印命令不执行；-v/--quiet/--log-level 控制日志详细程度。**日常用默认。**

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF
    │
    ▼
步骤0: (可选 --repair-vcf) 检查并修复损坏的 VCF
    │
    ▼
步骤1: 检查依赖 (bcftools)
    │
    ▼
步骤2: 分离 SNP / INDEL (bcftools view --types)
    │
    ▼
步骤3a: SNP 过滤 (QUAL/DP/MQ/QD/FS/SOR/MQRankSum/ReadPosRankSum/MAF)
    │
    ▼
步骤3b: INDEL 过滤 (宽松标准)
    │
    ▼
步骤3c: SNP 双等位位点过滤 (bcftools view -m2 -M2)
    │
    ▼
步骤4: 合并 SNP+INDEL (bcftools concat + sort)
    │
    ▼
步骤5: 统计报告
```

## 输出 | Output { #output }

```text
filtered_vcf/
├── variation_raw_snp.vcf.gz               # 分离出的原始 SNP (+.tbi)
├── variation_raw_indel.vcf.gz             # 分离出的原始 INDEL (+.tbi)
├── variation_filtered_snp.vcf.gz          # 过滤后的 SNP (+.tbi)
├── variation_filtered_indel.vcf.gz        # 过滤后的 INDEL (+.tbi)
├── variation_filtered_snp_biallelic.vcf.gz # 双等位位点过滤后的 SNP (+.tbi)
├── variation_filtered_merged.vcf.gz       # 合并后的最终 VCF (+.tbi)
└── (统计报告)
```

文件名前缀固定为 variation（模块内部 base_name 默认值）。

## 结果解读 | Interpreting Results { #interpreting }

### 1. variation_filtered_merged.vcf.gz（最终结果）

**通俗理解|In plain words:** 这是最终交付的、过滤+合并后的 SNP+INDEL VCF，直接用于下游。

### 2. 分类型文件

**通俗理解|In plain words:** raw_* 是分离后未过滤的，filtered_* 是过滤后的，biallelic_* 是再经过双等位过滤的 SNP。想单独用 SNP 或 INDEL 时取对应文件。

### 3. 统计报告

**通俗理解|In plain words:** 各步骤过滤前后的变异数量，看「过滤删掉了多少」。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规数据**：全部默认即可（GATK 最佳实践标准）。
- **低深度数据**：可下调 --snp-dp/--indel-dp（如 5）。
- **保留多等位位点**：用 --no-snp-biallelic。
- **只处理 SNP 或只处理 INDEL**：用 --variant-type snp_only / indel_only，跳过不必要的步骤。
- **VCF 报「列数不匹配」**：加 --repair-vcf 自动修复。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径(支持.vcf/.vcf.gz)｜Input VCF file path (supports .vcf/.vcf.gz) |
| `--output-dir, -o` | `./filtered_vcf` |  | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--snp-qual` | `30.0` | float | SNP最小质量值｜SNP minimum QUAL |
| `--snp-dp` | `10` | int | SNP最小测序深度｜SNP minimum DP |
| `--snp-mq` | `40.0` | float | SNP最小比对质量｜SNP minimum MQ |
| `--snp-qd` | `2.0` | float | SNP最小质量/深度比｜SNP minimum QD |
| `--snp-fs` | `60.0` | float | SNP最大FisherStrand值｜SNP maximum FS |
| `--snp-sor` | `3.0` | float | SNP最大StrandOddsRatio｜SNP maximum SOR |
| `--snp-mqrs` | `-12.5` | float | SNP最小MappingQualityRankSum｜SNP minimum MQRankSum |
| `--snp-rprs` | `-8.0` | float | SNP最小ReadPosRankSum｜SNP minimum ReadPosRankSum |
| `--snp-maf` | `0.05` | float | SNP最小次等位基因频率｜SNP minimum MAF |
| `--snp-biallelic` | `True` |  | 只保留双等位位点SNP｜Keep only biallelic SNPs |
| `--no-snp-biallelic` | `False` |  | 禁用双等位点过滤｜Disable biallelic filtering |
| `--indel-qual` | `30.0` | float | INDEL最小质量值｜INDEL minimum QUAL |
| `--indel-dp` | `10` | int | INDEL最小测序深度｜INDEL minimum DP |
| `--indel-mq` | `40.0` | float | INDEL最小比对质量｜INDEL minimum MQ |
| `--indel-qd` | `2.0` | float | INDEL最小质量/深度比｜INDEL minimum QD |
| `--indel-fs` | `200.0` | float | INDEL最大FisherStrand值｜INDEL maximum FS |
| `--indel-sor` | `10.0` | float | INDEL最大StrandOddsRatio｜INDEL maximum SOR |
| `--indel-rprs` | `-20.0` | float | INDEL最小ReadPosRankSum｜INDEL minimum ReadPosRankSum |
| `--variant-type` | `both` | both/snp_only/indel_only | 输入VCF文件的变异类型｜Variant type in input VCF |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |
| `--verbose, -v` | — |  | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level |
| `--force, -f` | — |  | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — |  | 模拟运行(不实际执行)｜Dry run without execution |
| `--repair-vcf` | — |  | 自动修复损坏的VCF文件（列数不匹配等问题）｜Auto-repair corrupted VCF files (column mismatch, etc.) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径(支持压缩和未压缩)｜Input VCF file path (supports compressed and uncompressed) |
| `-o, --output-dir` | `./filtered_vcf` |  | 输出目录｜Output directory (default: ./filtered_vcf) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--variant-type` | `both` | both/snp_only/indel_only | 输入VCF文件的变异类型｜Variant type in input VCF (default: both) |
| `--snp-qual` | `30.0` | float | SNP最小质量值｜SNP minimum QUAL (default: 30.0) |
| `--snp-dp` | `10` | int | SNP最小测序深度｜SNP minimum DP (default: 10) |
| `--snp-mq` | `40.0` | float | SNP最小比对质量｜SNP minimum MQ (default: 40.0) |
| `--snp-qd` | `2.0` | float | SNP最小质量/深度比｜SNP minimum QD (default: 2.0) |
| `--snp-fs` | `60.0` | float | SNP最大FisherStrand值｜SNP maximum FS (default: 60.0) |
| `--snp-sor` | `3.0` | float | SNP最大StrandOddsRatio｜SNP maximum SOR (default: 3.0) |
| `--snp-mqrs` | `-12.5` | float | SNP最小MappingQualityRankSum｜SNP minimum MQRankSum (default: -12.5) |
| `--snp-rprs` | `-8.0` | float | SNP最小ReadPosRankSum｜SNP minimum ReadPosRankSum (default: -8.0) |
| `--snp-maf` | `0.05` | float | SNP最小次等位基因频率｜SNP minimum MAF (default: 0.05) |
| `--snp-biallelic` | `True` | store_true | 只保留双等位位点SNP｜Keep only biallelic SNPs (default: True) |
| `--no-snp-biallelic` | — | store_false | 禁用双等位位点过滤｜Disable biallelic filtering |
| `--indel-qual` | `30.0` | float | INDEL最小质量值｜INDEL minimum QUAL (default: 30.0) |
| `--indel-dp` | `10` | int | INDEL最小测序深度｜INDEL minimum DP (default: 10) |
| `--indel-mq` | `40.0` | float | INDEL最小比对质量｜INDEL minimum MQ (default: 40.0) |
| `--indel-qd` | `2.0` | float | INDEL最小质量/深度比｜INDEL minimum QD (default: 2.0) |
| `--indel-fs` | `200.0` | float | INDEL最大FisherStrand值｜INDEL maximum FS (default: 200.0) |
| `--indel-sor` | `10.0` | float | INDEL最大StrandOddsRatio｜INDEL maximum SOR (default: 10.0) |
| `--indel-rprs` | `-20.0` | float | INDEL最小ReadPosRankSum｜INDEL minimum ReadPosRankSum (default: -20.0) |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path (default: bcftools) |
| `-v, --verbose` | `0` | count | 详细输出模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式(只输出ERROR)｜Quiet mode (ERROR only) |
| `--log-level` | — |  | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)｜Log level (default: INFO) |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行(不实际执行)｜Dry run without execution |
| `--repair-vcf` | — | store_true | 自动修复损坏的VCF文件（列数不匹配等问题）｜Auto-repair corrupted VCF files (column mismatch, etc.) |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- bcftools（bcftools，分离/过滤/合并/索引；可用 --bcftools-path 或 BCFTOOLS_PATH 覆盖）
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
不支持。本工具是单次过滤流程，每次运行都会重新执行各步骤；用 --force 覆盖已有输出。

**Q2：为什么 SNP 和 INDEL 的阈值不一样？**
INDEL 的链偏倚(FS/SOR)天然更高，所以用更宽松的标准（FS<=200、SOR<=10 vs SNP 的 FS<=60、SOR<=3），避免误删真实 INDEL。

**Q3：输入 VCF 报列数不匹配怎么办？**
加 --repair-vcf，程序会先检查并自动修复损坏的 VCF 再继续。

**Q4：过滤条件在哪看？**
运行时日志会逐条打印 SNP/INDEL 的完整过滤表达式；也可加 -v 看更详细的命令。

**Q5：输出文件名为什么都叫 variation？**
模块内部 base_name 固定为 variation，所以输出前缀统一是 variation_*，无法通过 CLI 修改。
