# VCF文件按染色体拆分工具|Split VCF File by Chromosome Tool

## 功能简介|Function Overview

将单个VCF文件按染色体拆分成多个文件，用于GEC基因组范围多重检验校正分析。

Split single VCF file into multiple files by chromosome for GEC genome-wide error correction analysis.

---

## 为什么需要拆分？|Why Split?

**GEC工具要求的VCF文件格式：**
- 文件名必须包含 `CHROM_` 占位符
- 例如：`chr_CHROM_.vcf.gz`
- GEC会自动将 `CHROM_` 替换为实际染色体编号（1, 2, 3...）

**如果你的VCF是单个文件**（如 `input.vcf.gz`），**必须先拆分**才能使用GEC！

---

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
# 拆分VCF文件
biopytools gwas-split-vcf \
  -i input.vcf.gz \
  -o vcf_split

# 拆分后使用GEC分析
biopytools gwas-gec \
  -i gwas_results.txt \
  -r vcf_split/chr_CHROM_.vcf.gz
```

### 完整示例|Complete Example

```bash
# 1. 拆分VCF文件
biopytools gwas-split-vcf \
  -i input.vcf.gz \
  -o vcf_by_chrom \
  -p chr_CHROM_

# 输出示例|Output example:
# vcf_by_chrom/chr_Chr01.vcf.gz
# vcf_by_chrom/chr_Chr02.vcf.gz
# vcf_by_chrom/chr_Chr03.vcf.gz
# ...

# 2. 使用拆分后的VCF进行GEC分析
biopytools gwas-gec \
  -i gwas_domain_p_file.txt \
  -r vcf_by_chrom/chr_CHROM_.vcf.gz \
  -t 64 \
  -m 200g
```

---

## 参数说明|Parameter Description

| 参数|Parameter | 说明|Description | 默认值|Default |
|------------|-------------|----------------------------------------|----------------|
| `-i, --input` | 输入VCF文件路径|Input VCF file path | 必需|required |
| `-o, --output-dir` | 输出目录|Output directory | 必需|required |
| `-p, --prefix` | 输出文件前缀|Output file prefix | `chr_CHROM_` |

### 关于prefix参数|About prefix parameter

- `CHROM_` 是占位符，会被替换为实际染色体编号
- `CHROM_` is a placeholder, will be replaced with actual chromosome ID

**示例|Examples:**
- `-p chr_CHROM_` → 生成 `chr_Chr01.vcf.gz`, `chr_Chr02.vcf.gz` ...
- `-p CHROM_` → 生成 `Chr01.vcf.gz`, `Chr02.vcf.gz` ...
- `-p chromosome_CHROM_` → 生成 `chromosome_Chr01.vcf.gz` ...

**重要|Important:** 必须确保前缀格式与P值文件的染色体格式一致！

---

## 输出文件|Output Files

### 文件命名规则|File Naming Rules

假设输入VCF包含染色体：`Chr01`, `Chr02`, `Chr03`

使用 `-p chr_CHROM_` 时，输出文件为：
```
vcf_split/
├── chr_Chr01.vcf.gz
├── chr_Chr02.vcf.gz
└── chr_Chr03.vcf.gz
```

### 文件内容|File Content

每个文件包含：
1. 完整的VCF表头（来自原始文件）
2. 该染色体的所有变异记录

Each file contains:
1. Complete VCF header (from original file)
2. All variant records for that chromosome

---

## 程序输出|Program Output

### 运行时信息|Runtime Information

```
======================================================================
VCF文件拆分工具|VCF File Splitter
======================================================================
输入文件|Input file: input.vcf.gz
输出目录|Output directory: vcf_split
文件前缀格式|File prefix format: chr_CHROM_
======================================================================

检测染色体格式|Detecting chromosome format...
  检测到的格式|Detected format: Chr_prefix
  示例染色体|Sample chromosomes: Chr01, Chr02, Chr03

读取VCF文件|Reading VCF file...
  已处理|Processed: 100,000 行|variants
  已处理|Processed: 200,000 行|variants
总共读取|Total read: 1,234,567 行|variants
发现染色体|Found chromosomes: 12 个

写入染色体文件|Writing chromosome files...
  Chr01: 98,765 行|variants -> chr_Chr01.vcf.gz
  Chr02: 87,654 行|variants -> chr_Chr02.vcf.gz
  ...

======================================================================
拆分完成|Splitting completed!
======================================================================

总变异数|Total variants: 1,234,567
染色体数量|Number of chromosomes: 12
输出文件数|Number of output files: 12

各染色体统计|Chromosome-wise statistics:
染色体|Chrom      变异数|Variants
------------------------------
Chr01          98,765
Chr02          87,654
...

======================================================================
使用GEC分析命令|GEC analysis command:
biopytools gwas-gec -i <your_pvalue_file.txt> -r vcf_split/chr_CHROM_.vcf.gz
======================================================================
```

---

## 实际案例|Real Case

### 你的问题|Your Issue

**输入文件：**
- P值文件：`gwas_domain_p_file.txt` (Chr01格式)
- VCF文件：`input.vcf.gz` (单个文件，Chr01格式)

**错误原因：**
```
GEC期望：多个VCF文件，文件名包含CHROM_占位符
你提供的：单个VCF文件 input.vcf.gz
结果：解析到0个LD块，无法计算阈值
```

**解决步骤：**

```bash
# Step 1: 拆分VCF文件
biopytools gwas-split-vcf \
  -i input.vcf.gz \
  -o vcf_by_chrom \
  -p chr_CHROM_

# Step 2: 使用拆分后的VCF运行GEC
biopytools gwas-gec \
  -i gwas_domain_p_file.txt \
  -r vcf_by_chrom/chr_CHROM_.vcf.gz \
  -t 64 \
  -m 200g
```

**预期结果：**
- 成功解析所有染色体
- 计算出正确的有效检验数
- 得到校正后的显著性阈值

---

## 常见问题|FAQ

### Q1: 拆分后的文件很大怎么办？

**A:** 可以使用 `--chrom` 参数只分析部分染色体：
```bash
# 只分析前3个染色体
biopytools gwas-gec \
  -i gwas.txt \
  -r vcf_split/chr_CHROM_.vcf.gz \
  --chrom Chr01-Chr03
```

### Q2: 我的染色体格式是 `1`, `2` 而不是 `Chr01`？

**A:** 工具会自动检测并适配，无需担心。输出文件名会自动匹配。

### Q3: 拆分需要多长时间？

**A:** 取决于VCF文件大小：
- 小文件（<1GB）：几分钟
- 大文件（>10GB）：可能需要1-2小时

### Q4: 能否使用bgzip压缩？

**A:** 拆分后的文件已经是gzip压缩格式（.gz），可以直接用于GEC分析。

---

## 技术细节|Technical Details

### 内存使用|Memory Usage

- 程序会将所有变异记录加载到内存
- 内存需求 ≈ VCF文件大小 × 2
- 对于超大VCF（>50GB），建议先按染色体提取

### 染色体排序|Chromosome Sorting

输出按以下顺序排列：
1. 数字染色体：1, 2, 3, ... 22
2. 性染色体：X, Y
3. 线粒体：M, MT
4. 其他：按字母顺序

---

## 更新日志|Changelog

| 版本|Version | 日期|Date | 更新内容|Updates |
|---------|------------|----------|----------|
| 1.0.0 | 2026-01-14 | 初始版本|Initial release |
