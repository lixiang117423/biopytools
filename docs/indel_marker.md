# 抗病/感病 INDEL 共显性标记开发 | R/S INDEL Codominant Marker Development

一句话理解：**从一群抗病个体和一群感病个体的测序数据里，自动找出「抗病个体都携带、感病个体都没有(或反过来)」的插入缺失(INDEL)位点，并给每个候选位点设计好 PCR 引物**，解决「找到能区分抗病/感病表型的分子标记」的问题。

## 功能概述 | Overview { #overview }

- 输入多样本合并 VCF + 样品分组表(samplesheet) + 参考基因组 FASTA
- 提取 INDEL → 群体基因型判定(共显性) → 覆盖度验证(含 deletion 骤降) → 提取侧翼序列 → 设计 PCR 引物
- 输出候选主表、BED 区间、侧翼 FASTA、覆盖度矩阵和摘要报告
- 断点续传：VCF 提取步骤已产出矩阵时自动跳过
- 依赖 bcftools / samtools / primer3-py 三个工具

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools indel-marker -v v.vcf.gz -s samplesheet.tsv -g ref.fa -o out/
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| INDEL | 基因组上的小片段「插入(insertion)/缺失(deletion)」，像一句话里多写或少写几个字 |
| 共显性标记 | 能区分「纯合变异/杂合/纯合参考」三种状态的标记，信息量比只区分有无更高 |
| 抗病/感病(resistant/susceptible) | 两个表型分组，程序要找的就是只在其中一组出现的 INDEL |
| 基因型(genotype) | 一个位点上两条同源染色体的构成：hom_ref=纯合参考、het=杂合、hom_alt=纯合变异 |
| 覆盖度(depth) | 一个位点被测到的 reads 数，像「这个位置被读了多少遍」 |
| deletion 骤降 | 携带缺失的样本该区域覆盖度应明显低于不携带的样本，用于验证缺失真实性 |
| 侧翼序列(flank) | INDEL 两侧的一段序列，引物就设计在这里 |
| PCR 引物 | 一小段短 DNA，像「夹住目标片段两端的夹子」，用于实验室放大目标片段 |

## 输入 | Input { #input }

### VCF 文件

多样本合并的 VCF(支持 `.vcf` / `.vcf.gz`)。用 bcftools 提取 INDEL 并导出每个样本的 GT 基因型矩阵。

### samplesheet 样品分组表

三列：`sample_name  group  bam_path`。无需表头(首行第二列是 `group` 时自动识别并跳过)；默认 TAB 分隔，无 TAB 时回退为空白分隔；任意两种不同的分组标签都会自动映射为抗病/感病。

```text
sample1    resistant    /data/bam/sample1.sorted.bam
sample2    resistant    /data/bam/sample2.sorted.bam
sample3    susceptible  /data/bam/sample3.sorted.bam
sample4    susceptible  /data/bam/sample4.sorted.bam
```

### 参考基因组 FASTA

用于提取侧翼序列(`-g` 指定)。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-v` 是多样本合并 VCF，`-s` 是样品分组表，`-g` 是参考基因组，`-o` 是输出目录。samplesheet 里的样本名必须都能在 VCF 里找到，且抗病、感病两组都至少有一个样本。

### INDEL 过滤参数 | INDEL filtering

**通俗理解|In plain words:** 这一组决定「什么样的 INDEL 才进候选」。`--min-indel-size`(默认 10)/`--max-indel-size`(默认 100)限定 INDEL 长度，太短难设计引物、太长难扩增；`--min-quality`(默认 20.0)是最低 QUAL 值，缺失 QUAL 的保留；`--max-candidates`(默认 0=不限)是候选数上限，候选太多时用来截断。**默认值适合常规 PCR 标记，一般不用动。**

### 群体判定参数 | Population calling

**通俗理解|In plain words:** `--min-group-consistency`(默认 0.9)是「组内有多少比例的样本要一致」的阈值，1.0=最严格(要求整组完全一致)，调低能容忍个别离群样本但假阳性变多；`--min-samples-per-group`(默认 1)是每组最少样本数。**默认 0.9 是常用折中，一般不用动。**

### 覆盖度参数 | Coverage

**通俗理解|In plain words:** `--min-depth`(默认 10)是最低覆盖度，所有样本都要达到才算通过质控；`--deletion-depth-ratio`(默认 0.3)是 deletion 骤降阈值，携带缺失组的平均覆盖度低于不携带组的 0.3 倍才算「确实缺失」。**默认值经实践验证，一般不用动。**

### 引物参数 | Primer

**通俗理解|In plain words:** `--flank-length`(默认 300)是 INDEL 两侧各取多长的序列用于设计引物，取短了引物没地方放，取长了设计慢。**默认 300 即可，一般不用动。**

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 是线程数(默认 12)。覆盖度计算是逐个样本串行跑的，线程数目前主要影响外部工具调用，一般不用特别调。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先筛出合格的 INDEL，再看它们在抗病/感病两组里是否「一组有、一组没有」，然后验证覆盖度、取侧翼序列、设计引物，最后打包输出。

```text
输入 VCF + samplesheet + 参考 FASTA
    |
    ▼
步骤1: 依赖检查(bcftools / samtools / primer3-py)
    |
    ▼
步骤2: 解析 samplesheet + 校验分组
    |
    ▼
步骤3: 提取 INDEL 矩阵(bcftools view -v indels | query → 01_vcf_extract/)
    |
    ▼
步骤4: 群体基因型判定(一组 present + 另一组 absent 才保留)
    |
    ▼
步骤5: 覆盖度计算(samtools depth) + 质控 + deletion 骤降验证
    |
    ▼
步骤6: 提取侧翼序列(samtools faidx)
    |
    ▼
步骤7: 设计 PCR 引物(primer3-py)
    |
    ▼
步骤8: 写出候选主表/BED/FASTA/覆盖度矩阵/摘要
```

## 输出 | Output { #output }

```text
out/
├── 00_pipeline_info/
│   ├── software_versions.yml    # 软件版本与运行参数
│   └── pipeline_params.yaml     # 关键运行参数
├── 01_vcf_extract/
│   └── indels_gt_matrix.tsv     # 提取的 INDEL 基因型矩阵
├── 02_genotype_call/
├── 03_coverage/
│   └── indels_coverage.tsv      # 每候选每样本的平均覆盖度矩阵
├── 04_sequence/
│   └── indels_flank.fa          # 每候选的侧翼序列(FASTA)
├── 05_primer/
├── 06_results/
│   ├── indel_marker_candidates.tsv  # 候选主表(核心结果)
│   ├── indel_marker_candidates.bed  # 候选区间 BED
│   └── indel_marker_summary.txt     # 筛选摘要
└── 99_logs/
    └── indel_marker.log         # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. indel_marker_candidates.tsv(候选主表，核心)

**通俗理解|In plain words:** 一张「哪些 INDEL 能当标记」的最终清单，一行一个候选，含方向、两组基因型比例、覆盖度、引物序列。

关键列：

- `candidate_id`：候选编号，格式 `{染色体}:{起点}-{终点}:{DEL|INS}:{R|S}_spec`，末尾 `R`=抗病特异、`S`=感病特异
- `indel_type` / `indel_size`：插入还是缺失、多长
- `R_hom_alt_rate` / `S_hom_alt_rate`：抗病/感病组里「纯合变异」的比例，高的一方就是「携带组」
- `passes_coverage_qc`：是否通过最低覆盖度质控；`passes_deletion_drop`：是否通过 deletion 骤降验证(insertion 为 NA)
- `left_primer` / `right_primer` / `product_size`：设计出的左右引物和产物长度；`primer_status` 为 `ok` 表示设计成功

### 2. indel_marker_candidates.bed

**通俗理解|In plain words:** 候选区间转成 BED 格式(0-based 起点)，可直接丢给 IGV 等可视化工具看。

### 3. indels_flank.fa

**通俗理解|In plain words:** 每个候选的侧翼序列，FASTA 头的名字就是 candidate_id，用于后续引物设计或人工核对。

### 4. indels_coverage.tsv

**通俗理解|In plain words:** 每个候选在每个样本里的平均覆盖度矩阵，可核对「携带组覆盖度确实更低」的预期。

### 5. indel_marker_summary.txt

**通俗理解|In plain words:** 一份数字汇总：合格 INDEL 总数、群体判定候选数、通过覆盖度质控数、成功设计引物数。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规 PCR 标记开发**：全部用默认值即可
- **要严格排除个别离群样本**：`--min-group-consistency 1.0`(要求整组完全一致)
- **群体有少量杂株/污染**：适当调低 `--min-group-consistency`(如 0.8)，容忍离群
- **只想快速预览前 N 个候选**：`--max-candidates 50`
- **缺失片段较大**：`--max-indel-size` 调大；若只想看小片段，`--min-indel-size` 调大

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-v, --vcf` | 必填 |  | 多样本合并VCF｜Multi-sample VCF |
| `-s, --samplesheet` | 必填 |  | 样品分组TSV(sample_name/group/bam_path)｜Samplesheet TSV |
| `-g, --genome-fasta` | 必填 |  | 参考基因组｜Reference FASTA |
| `-o, --output-dir` | `./indel_marker_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--min-indel-size` | `10` |  | 最小INDEL长度｜Min INDEL size |
| `--max-indel-size` | `100` |  | 最大INDEL长度｜Max INDEL size |
| `--min-quality` | `20.0` |  | 最低QUAL过滤(缺失QUAL保留)｜Min QUAL filter (missing QUAL kept) |
| `--max-candidates` | `0` |  | 候选数上限(0=不限)｜Candidate cap (0=no limit) |
| `--min-group-consistency` | `0.9` |  | 组内纯合一致比例(1.0=严格)｜Within-group consistency |
| `--min-samples-per-group` | `1` |  | 每组最少样品数(默认1)｜Min samples per group (default 1) |
| `--min-depth` | `10` |  | 最低覆盖度｜Min depth |
| `--deletion-depth-ratio` | `0.3` |  | deletion骤降阈值｜deletion drop threshold |
| `--flank-length` | `300` |  | 侧翼长度｜Flank length |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-v, --vcf` | 必填 |  | 多样本合并VCF｜Multi-sample VCF |
| `-s, --samplesheet` | 必填 |  | 样品分组TSV(sample_name/group/bam_path)｜Samplesheet TSV |
| `-g, --genome-fasta` | 必填 |  | 参考基因组｜Reference FASTA |
| `-o, --output-dir` | `./indel_marker_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--min-indel-size` | `10` | int | 最小INDEL长度｜Min INDEL size |
| `--max-indel-size` | `100` | int | 最大INDEL长度｜Max INDEL size |
| `--min-quality` | `20.0` | float | 最低QUAL过滤(缺失QUAL保留)｜Min QUAL filter (missing QUAL kept) |
| `--max-candidates` | `0` | int | 候选数上限(0=不限)｜Candidate cap (0=no limit) |
| `--min-group-consistency` | `0.9` | float | 组内纯合一致比例阈值｜Min within-group consistency (1.0=strict) |
| `--min-samples-per-group` | `1` | int | 每组最少样品数(默认1)｜Min samples per group (default 1) |
| `--min-depth` | `10` | int | 最低覆盖度｜Min depth |
| `--deletion-depth-ratio` | `0.3` | float | deletion骤降阈值｜deletion drop threshold |
| `--flank-length` | `300` | int | 侧翼长度｜Flank length |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- bcftools(提取 INDEL 与基因型矩阵)
- samtools(覆盖度计算 + 侧翼序列提取)
- primer3-py(Python 包，设计 PCR 引物；`pip install primer3-py`)

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持 VCF 提取这一最耗时步骤：`01_vcf_extract/indels_gt_matrix.tsv` 已存在则跳过。其余步骤(覆盖度/引物)每次重跑。换过滤参数重跑前需先删除该矩阵文件。

**Q2：引物那列全是空的？**
多半是 primer3-py 没装，或本地仓库里的 `biopytools/primer3/` 子包覆盖了 `primer3` 包名。请确认 `pip install primer3-py` 且运行目录不含同名子包。

**Q3：samplesheet 一定要写 resistant/susceptible 吗？**
不用。任意两种不同的分组标签都可以，会自动映射为抗病/感病(标准别名如 R/S、抗病/感病会直接识别，其它标签按出现顺序：第 1 组=抗病、第 2 组=感病)。

**Q4：deletion 骤降为什么有的候选是 NA？**
只有长度达到 30bp 的 deletion 才做骤降验证；insertion 和较短 deletion 该字段为 NA，属正常。

**Q5：为什么没有候选输出？**
常见原因：两组都携带(数据矛盾)或两组都不携带(不分离)的 INDEL 会被跳过；也可能 QUAL 太低或长度超出范围。检查日志里的「群体判定候选」计数。
