# INDEL分子标记开发 | INDEL Marker Development

**从抗病/感病群体VCF开发共显性INDEL分子标记, 含群体判定+覆盖度质控+引物设计 | Develop codominant INDEL markers from R/S bulks: genotype call + coverage QC + primer design**

## 功能概述 | Overview

indel_marker 模块针对抗病(R)/感病(S)群体, 从多样本合并 VCF 中筛选适合育种应用的共显性 INDEL 标记。完整流程:

1. **VCF 提取 INDEL** — 按长度/质量过滤, 生成基因型矩阵
2. **群体共显性判定** — R/S 两组在位点上的基因型一致性, 选共显性候选
3. **覆盖度质控** — 用 BAM 检查候选位点在每样品的 depth, 剔除 deletion 骤降等假阳性
4. **侧翼序列提取** — 从参考基因组提取候选位点侧翼
5. **引物设计** — primer3 设计引物
6. **输出** — 候选主表 + BED + FASTA + 覆盖度表 + summary + 流程元数据

支持断点续传(各步骤输出存在则跳过)。

## 快速开始 | Quick Start

```bash
# 基础用法
biopytools indel-marker -v variants.vcf.gz -s samplesheet.tsv -g ref.fa -o out/

# samplesheet.tsv 格式: sample_name<TAB>group(R或S)<TAB>bam_path
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-v, --vcf` | 多样本合并 VCF |
| `-s, --samplesheet` | 样品分组 TSV(sample_name/group/bam_path) |
| `-g, --genome-fasta` | 参考基因组 FASTA |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--min-indel-size` | `10` | 最小 INDEL 长度 |
| `--max-indel-size` | `100` | 最大 INDEL 长度 |
| `--min-quality` | `20.0` | 最低 QUAL 过滤(缺失 QUAL 保留) |
| `--max-candidates` | `0` | 候选数上限(0=不限) |
| `--min-group-consistency` | `0.9` | 组内纯合一致比例阈值(1.0=严格) |
| `--min-samples-per-group` | `1` | 每组最少样品数 |
| `--min-depth` | `10` | 最低覆盖度 |
| `--deletion-depth-ratio` | `0.3` | deletion 骤降阈值 |
| `--flank-length` | `300` | 侧翼长度 |

(运行 `biopytools indel-marker -h` 查看完整参数列表)

## 输出 | Output

- 候选标记主表 TSV(含引物)
- 候选位点 BED
- 侧翼序列 FASTA
- 覆盖度 TSV(每样品 × 候选位点)
- summary 汇总 + `00_pipeline_info/` 流程元数据(软件版本/参数)

## 断点续传 | Checkpoint Resume

各步骤输出文件已存在时跳过(如 INDEL 提取矩阵存在则跳过提取)。

## 依赖 | Dependencies

- **bcftools**: VCF 提取与基因型
- **samtools**: BAM 覆盖度统计
- **primer3-py**: 引物设计

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
