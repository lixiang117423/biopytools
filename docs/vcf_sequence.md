# vcf-sequence - VCF序列提取 | VCF Sequence Extraction

一句话理解：**给定参考基因组 FASTA 和 VCF，把指定染色体区间里「每个样本自己的那套」DNA 序列重建出来**（把变异套到参考序列上），连同参考序列一起导出成表格/FASTA/CSV，方便逐碱基比较不同样本。

## 功能概述 | Overview { #overview }

- 从参考基因组提取指定区间的参考序列
- 读取该区间内的 VCF 变异
- 把变异「套」到参考序列上，重建每个样本的序列
- 支持 tab / fasta / csv 三种输出格式
- 生成碱基组成统计报告和总结报告
- 支持指定/排除样本、质量过滤、选第一或第二等位基因

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-sequence -v variants.vcf -g genome.fa -c chr1 -s 1000000 -e 1001000 -o seq_output
```

提取 chr1 上 1000000-1001000 区间，重建每个样本的序列并导出到 `seq_output/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 参考序列 | 群体共同参照的那条「标准 DNA 序列」，来自参考基因组 |
| 变异 | 某个样本在某个位置与参考不同的地方（如 A 变成了 G） |
| 等位基因 | 一个位置上两条染色体各自的版本；第一/第二等位基因即其中一条 |
| 坐标 | 染色体上的位置编号，本工具按 VCF 惯例用 1 起算 |
| FASTA | 最通用的序列存储格式，`>` 开头是名字，下面一行是碱基序列 |

## 输入 | Input { #input }

- **VCF 文件**：包含区间内变异的 VCF（需已建立索引 .tbi/.csi，见 FAQ）。
- **基因组 FASTA 文件**：参考基因组，需能被 pysam 读取（未建 .fai 索引时 pysam 会自动建立）。
- **区间坐标**：`-c` 染色体名（须与 VCF/FASTA 里一致）、`-s` 起始位置、`-e` 结束位置（1 起算，`start < end`，`start >= 1`）。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 五个参数（VCF、基因组、染色体、起止位置）缺一不可，共同确定「从哪拿参考、把谁的变异套上去」。

### 输出 | Output

**通俗理解|In plain words:** `-o` 指定输出目录，`--format` 选择导出格式（tab 表 / fasta 序列 / csv 表）。**默认 tab 最适合人眼看和 Excel 打开；要喂给比对软件就选 fasta。**

### 序列构建 | Sequence building

**通俗理解|In plain words:** `--second-allele` 让每个位点用「第二个等位基因」而不是默认的第一个（看另一条染色体的版本）；`--no-reference` 让结果里不包含参考序列、只留样本。**一般用默认即可。**

### 过滤 | Filtering

**通俗理解|In plain words:** `--min-qual` 按位点质量分过滤低质量变异；`--samples` 只处理指定样本、`--exclude-samples` 排除某些样本（都支持逗号分隔的列表或一个文件，每行一个样本名）。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF + 基因组 FASTA + 区间
    │
    ▼
步骤1: 打开基因组和 VCF(pysam)
    │
    ▼
步骤2: 从 FASTA 取参考序列
    │
    ▼
步骤3: 读取样本名(应用指定/排除样本过滤)
    │
    ▼
步骤4: 读取区间内变异(应用质量过滤)
    │
    ▼
步骤5: 逐样本把变异套到参考序列上重建序列(只处理单碱基替换)
    │
    ▼
步骤6-8: 导出序列 + 生成统计报告 + 总结报告
```

## 输出 | Output { #output }

以 `-c chr1 -s 1000000 -e 1001000 --format tab` 为例：

```text
seq_output/
├── chr1_1000000_1001000_sequences.txt   # 序列文件(tab/fasta/csv 由 --format 决定)
├── chr1_1000000_1001000_statistics.txt  # 碱基组成统计(A/T/C/G/N + GC含量)
├── extraction_summary.txt               # 总结报告(输入/处理/过滤信息)
└── sequence_extraction.log              # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **序列文件**（`*_sequences.txt` 等）：一行/一条记录对应一个样本，第一行（默认）是 `Reference` 参考序列，其后是按样本名排序的样本序列。逐列/逐碱基对比即可看到各样本与参考的差异。
- **统计报告**（`*_statistics.txt`）：参考序列和每个样本的 A/T/C/G/N 计数与 GC 含量，缺失位点记作 N。
- **总结报告**（`extraction_summary.txt`）：记录本次提取的区间、样本数、变异数、使用的等位基因与过滤参数，方便归档。
- 注意：缺失基因型会写成 `N`；插入/删除(indel)变异不会套用到序列上（只处理单碱基替换）。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **区间长度**：单次提取区间控制在 10kb 以内（代码校验，超过约 10kb 会报错终止），更长区域请分段提取。
- **`--format`**：给人看/进 Excel 选 `tab` 或 `csv`；要拿序列去比对或建树选 `fasta`。
- **`--second-allele`**：想看每个样本「另一条染色体」的序列时使用（默认取第一个等位基因）。
- **`--samples/--exclude-samples`**：样本很多、只关心其中一部分时用它缩小结果。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--vcf, -v` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--genome, -g` | 必填 | Path | 基因组FASTA文件路径｜Genome FASTA file path |
| `--chrom, -c` | 必填 | str | 染色体名称｜Chromosome name |
| `--start, -s` | 必填 | int | 起始位置｜Start position |
| `--end, -e` | 必填 | int | 结束位置｜End position |
| `--output-dir, -o` | `./sequence_output` | Path | 输出目录｜Output directory |
| `--format` | `tab` | tab/fasta/csv | 输出格式｜Output format |
| `--second-allele` | — |  | 使用第二个等位基因｜Use second allele |
| `--no-reference` | — |  | 不包含参考序列｜Do not include reference sequence |
| `--min-qual` | — | int | 最小质量过滤｜Minimum quality filter |
| `--samples` | — | str | 指定样本｜Specify samples |
| `--exclude-samples` | — | str | 排除样本｜Exclude samples |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-v, --vcf` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-g, --genome` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `-c, --chrom` | 必填 |  | 染色体名称｜Chromosome name |
| `-s, --start` | 必填 | int | 起始位置｜Start position |
| `-e, --end` | 必填 | int | 结束位置｜End position |
| `-o, --output-dir` | `./sequence_output` |  | 输出目录｜Output directory |
| `--format` | `tab` | tab/fasta/csv | 输出格式｜Output format |
| `--second-allele` | — | store_true | 使用第二个等位基因｜Use second allele |
| `--no-reference` | — | store_true | 不包含参考序列｜Do not include reference sequence |
| `--min-qual` | — | int | 最小质量值过滤｜Minimum quality filter |
| `--samples` | — |  | 指定样本列表｜Specify samples |
| `--exclude-samples` | — |  | 排除样本列表｜Exclude samples |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- pysam（必需，读取 FASTA 与 VCF；未安装会直接报错退出）

## 常见问题 | FAQ { #faq }

**Q1：运行报「区间长度过大」？**
代码限制单次区间超过约 10kb 会报错。把大区间拆成多段分别提取。

**Q2：VCF 或基因组没索引会怎样？**
VCF 的区间提取依赖索引（.tbi/.csi），缺失时会报错，先用 `bcftools index` 或 `tabix -p vcf` 建立；基因组 FASTA 的 .fai 索引通常由 pysam 自动建立。

**Q3：indel（插入/删除）会体现在序列里吗？**
不会。本工具只把单碱基替换(SNP)套到参考序列上，indel 会被跳过，相关位置保持参考序列。

**Q4：缺失基因型怎么表示？**
该位置写成 `N`，表示「不知道是什么碱基」。

**Q5：有断点续传吗？**
没有。每次运行都重新提取并覆盖输出目录里的结果文件。
