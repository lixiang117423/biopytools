# 覆盖度过滤工具|Coverage Filter Tool

## 功能描述|Function Description

基于BAM文件的覆盖度对基因组序列进行质量分级和过滤，识别高质量序列和可疑的污染/低质量序列。

Sequence quality classification and filtering based on BAM coverage, identifying high-quality sequences and suspicious contaminated/low-quality sequences.

## 安装依赖|Dependencies

- samtools
- seqtk
- seqkit

## 使用方法|Usage

### 基本用法|Basic Usage

```bash
biopytools coverage-filter -i sample.bam -f genome.fa -o filtered
biopytools coverage-filter -i sample.bam -f genome.fa -o filtered -d output_dir
```

### 参数说明|Parameters

| 参数|Parameter | 说明|Description |
|------------|------|------|
| `-i`, `--bam-file` | BAM文件路径|BAM file path |
| `-f`, `--fasta-file` | 基因组FASTA文件|Genome FASTA file |
| `-o`, `--output-prefix` | 输出文件前缀|Output file prefix |
| `-d`, `--output-dir` | 输出目录(默认: 当前目录)|Output directory (default: current directory) |
| `-t`, `--threads` | 线程数(默认: 12)|Number of threads |
| `--high-cov` | 高质量覆盖度阈值(默认: 90.0)|High quality coverage threshold |
| `--medium-cov-min` | 中等质量最小覆盖度(默认: 50.0)|Medium quality minimum coverage |

### 输出文件|Output Files

所有文件保存在输出目录中（`-d/--output-dir`指定，默认为当前目录）

| 文件|File | 说明|Description |
|----------|------|------|
| `{prefix}_high_quality.fa` | 高质量序列，推荐用于后续分析|High quality sequences, recommended |
| `{prefix}_medium_quality.fa` | 中等质量序列|Medium quality sequences |
| `{prefix}_low_quality.fa` | 低质量/可疑序列|Low quality/suspicious sequences |
| `{prefix}_coverage.txt` | 详细覆盖度数据|Detailed coverage data |
| `{prefix}_high_quality.list` | 高质量序列ID列表|High quality sequence ID list |
| `{prefix}_medium_quality.list` | 中等质量序列ID列表|Medium quality sequence ID list |
| `{prefix}_low_quality.list` | 低质量序列ID列表|Low quality sequence ID list |
| `{prefix}_report.txt` | 完整统计报告|Complete statistics report |

### 质量分级标准|Quality Classification Standards

#### 高质量序列|High Quality Sequences
- **覆盖度** ≥ `--high-cov` (默认: 90%)

#### 中等质量序列|Medium Quality Sequences
- **覆盖度**: `--medium-cov-min` 到 `--high-cov` (默认: 50% - 90%)

#### 低质量/可疑序列|Low Quality/Suspicious Sequences
- **覆盖度** < `--medium-cov-min` (默认: < 50%)
- 可能为污染序列或组装错误|Possible contamination or assembly errors

## 使用示例|Examples

### 示例1: 基本用法|Example 1: Basic Usage

```bash
biopytools coverage-filter \
  -i sample.bam \
  -f genome.fa \
  -o filtered
```

### 示例2: 自定义阈值和输出目录|Example 2: Custom Thresholds and Output Directory

```bash
biopytools coverage-filter \
  -i sample.bam \
  -f genome.fa \
  -o filtered \
  -d results \
  --high-cov 95.0 \
  --medium-cov-min 60.0
```

## 工作流程|Workflow

1. **计算覆盖度|Calculate Coverage**: 使用 `samtools coverage` 计算每个序列的覆盖度
2. **分类序列|Classify Sequences**: 根据覆盖度阈值对序列进行分级
3. **提取序列|Extract Sequences**: 使用 `seqtk subseq` 提取各质量等级的序列
4. **生成报告|Generate Report**: 生成详细的统计报告

## 注意事项|Notes

- 确保BAM文件已排序并建立索引|Ensure BAM file is sorted and indexed
- FASTA文件中的序列ID需要与BAM文件中的参考序列一致|Sequence IDs in FASTA must match BAM reference sequences
- 低质量序列可能需要进一步BLAST检查|Low quality sequences may need further BLAST inspection

## 版本|Version

v1.0.0 (2026-01-26)
