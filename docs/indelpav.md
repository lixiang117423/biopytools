# INDEL PAV 分析 | INDEL Presence/Absence Variation Analysis

一句话理解：**把一个多样本 VCF 里的所有插入/缺失(INDEL)转成一张「0/1 存在缺失」矩阵**，直接告诉你每个样本里有哪些 INDEL 有、哪些没有，解决「快速看清群体里 INDEL 的分布」的问题。

## 功能概述 | Overview { #overview }

- 输入一个 VCF(支持 `.vcf` / `.vcf.gz`)，自动提取所有 INDEL
- 排除 SNP 与复杂变异(多等位、`*` 等位)，按长度/质量/深度/缺失率过滤
- 输出 PAV 矩阵(每样本 0/1，1=存在、0=不存在或缺失) + 摘要报告
- 支持 `--compress` 压缩输出

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools indelpav -v variants.vcf -o pav_results.txt
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| PAV | Presence/Absence Variation，存在/缺失变异，只关心「有没有」，不关心具体多少 |
| INDEL | 小片段插入(insertion)/缺失(deletion) |
| 缺失率(missing rate) | 多少样本在这个 INDEL 上「没检出」，缺失太多说明不可靠 |
| 0/1 矩阵 | 每行一个 INDEL、每列一个样本，1=该样本有、0=没有，像一张勾选表 |
| QUAL | VCF 里每个变异的置信度打分，越高越可信 |

## 输入 | Input { #input }

### VCF 文件

多样本 VCF(压缩或未压缩均可)。程序从 `#CHROM` 行读取样本名，然后逐行处理变异。

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2 ...
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-v` 是输入 VCF，`-o` 是输出文件(默认 `./indel_pav.txt`)。只给 `-v` 就能跑。

### 过滤参数 | Filtering

**通俗理解|In plain words:** 这一组决定「什么样的 INDEL 才写进结果」。`--min-length`(默认 1)/`--max-length`(默认不限制)限定 INDEL 长度；`-q/--min-quality`(默认 20.0)是最低 QUAL，低于就丢弃；`-d/--min-depth`(默认 5)是最低深度(从 INFO 的 DP 字段读取，无 DP 则不查)；`--max-missing`(默认 0.8)是「未检出该 INDEL 的样本比例」上限。**默认值较宽松，一般不用动。**

### 输出参数 | Output

**通俗理解|In plain words:** `--compress` 或 `-o` 文件名以 `.gz` 结尾时输出 gzip 压缩文件；`--include-complex` 用于包含复杂变异，但当前版本对复杂变异(多等位/`*`)始终跳过，此参数暂不影响结果。`-t` 线程数(默认 12)为预留参数，当前按单线程顺序处理。`--bcftools-path`(默认 `bcftools`)指定 bcftools 路径，仅在依赖检查时用到。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先读样本名，再逐行判断是不是 INDEL、过不过滤，最后把每个样本的 0/1 拼成一张表输出。

```text
输入 VCF
    |
    ▼
步骤1: 依赖检查(bcftools)
    |
    ▼
步骤2: 解析 VCF 头部(#CHROM 行读样本名)
    |
    ▼
步骤3: 逐行处理变异(判断 INDEL → 长度/质量/深度/缺失率过滤 → 提取序列 → 解析 0/1 基因型)
    |
    ▼
步骤4: 结果过滤(PAVFilter)
    |
    ▼
步骤5: 写出 PAV 矩阵 + 摘要报告
```

## 输出 | Output { #output }

```text
./
├── indel_pav.txt              # PAV 矩阵(核心结果)
├── indel_pav_summary.txt      # 摘要报告
└── indel_pav_analysis.log     # 运行日志
```

### PAV 矩阵格式

```text
Chromosome    Start    End    REF    ALT    INDEL_Type    sample1    sample2 ...
Chr1          100      105    ATG    A      deletion      1          0 ...
Chr1          200      201    A      AT     insertion     0          1 ...
```

- 前 6 列是 INDEL 的坐标与信息，之后每列一个样本，值为 `0`(不存在/缺失)或 `1`(存在)
- `INDEL_Type` 为 `insertion` 或 `deletion`；`ALT` 列是该 INDEL 的插入/缺失序列

## 结果解读 | Interpreting Results { #interpreting }

### 1. PAV 矩阵(核心)

**通俗理解|In plain words:** 一张「谁有谁没有」的大表，可直接拿去做群体分布统计、聚类或画热图。

- `1` = 该样本携带这个 INDEL(至少一个非参考等位基因)；`0` = 不携带或该位点缺失
- 某列几乎全 1 或全 0 的 INDEL 信息量低，可结合下游分析再筛选

### 2. indel_pav_summary.txt(摘要)

**通俗理解|In plain words:** 一份数字汇总，含样本数、INDEL 总数、插入/删除占比、平均长度、长度范围和染色体分布。

## 参数选择建议 | Parameter Guidance { #guidance }

- **快速看整体 PAV**：默认参数直接跑
- **只要可信度高的 INDEL**：调高 `-q`(如 30)和 `-d`(如 10)
- **只要大片段 INDEL**：`--min-length 50`；只要小片段则 `--max-length 100`
- **样本缺失严重、想收紧**：调低 `--max-missing`(如 0.5)
- **结果文件大**：加 `--compress` 压缩输出

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -v` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `./indel_pav.txt` | Path | 输出文件路径｜Output file path |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-length` | `1` | int | 最小INDEL长度(bp)｜Minimum INDEL length (bp) |
| `--max-length` | — | int | 最大INDEL长度(bp)｜Maximum INDEL length (bp) |
| `--min-quality, -q` | `20.0` | float | 最小质量分数｜Minimum quality score |
| `--min-depth, -d` | `5` | int | 最小深度｜Minimum depth |
| `--max-missing` | `0.8` | float | 最大缺失率(0-1)｜Maximum missing rate (0-1) |
| `--include-complex` | — |  | 包含复杂变异｜Include complex variants |
| `--compress` | — |  | 压缩输出文件｜Compress output file |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-v, --vcf` | 必填 |  | 输入VCF文件路径 (支持压缩和未压缩)｜Input VCF file path (supports compressed and uncompressed) |
| `-o, --output` | `./indel_pav.txt` |  | 输出文件路径｜Output file path |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--min-length` | `1` | int | 最小INDEL长度(bp)｜Minimum INDEL length (bp) |
| `--max-length` | — | int | 最大INDEL长度(bp)｜Maximum INDEL length (bp) |
| `-q, --min-quality` | `20.0` | float | 最小质量分数｜Minimum quality score |
| `-d, --min-depth` | `5` | int | 最小深度｜Minimum depth |
| `--max-missing` | `0.8` | float | 最大缺失率 (0-1)｜Maximum missing rate (0-1) |
| `--include-complex` | — | store_true | 包含复杂变异｜Include complex variants |
| `--compress` | — | store_true | 压缩输出文件｜Compress output file |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- bcftools(依赖检查用 `--version` 探测；自动解析 align 域环境并经 conda run 调用，可用 --bcftools-path 或环境变量 BCFTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用)

## 常见问题 | FAQ { #faq }

**Q1：`--include-complex` 加了没效果？**
当前版本的 INDEL 判定始终排除多等位(ALT 含逗号)和 `*` 等位的复杂变异，`--include-complex` 参数虽被接受但暂未改变这一行为。

**Q2：`-t` 线程数为什么感觉没加速？**
当前版本按单线程顺序逐行处理 VCF，线程数为预留参数，暂不影响并行度。

**Q3：`--max-missing` 到底算的是什么？**
是「基因型为 0(不含该 INDEL)的样本比例」，缺失和真实不存在都算在内；超过该比例就丢弃该 INDEL。

**Q4：输出想压缩怎么处理？**
加 `--compress`，或让 `-o` 文件名以 `.gz` 结尾，两者都会输出 gzip 压缩文件。

**Q5：没有输出结果？**
常见原因：VCF 里没有合格的 INDEL(全是 SNP 或被过滤光)。检查日志里的「INDEL数量」和「过滤掉」计数。
