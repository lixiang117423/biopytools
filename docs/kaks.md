# Ka/Ks 计算 | Ka/Ks Calculation

**计算同源基因对之间的非同义替换率 (Ka) 与同义替换率 (Ks) | Calculate non-synonymous (Ka) and synonymous (Ks) substitution rates between homologous gene pairs.**

## 功能概述 | Overview

`kaks` 模块封装了 `KaKs_Calculator 2.0`，给定两个物种的 CDS FASTA 文件和一对配对列表，自动完成：seqkit 子集抽取 → 序列校验 → 生成 AXT 输入 → KaKs_Calculator 计算 → 结果解析、统计与多格式导出。支持 17 种替换速率计算方法（GMYN、YN、NG、LWL、GY 等），并自动生成选择压力分布（正选择/中性/纯化选择）与生物学解释摘要。

## 快速开始 | Quick Start

```bash
biopytools kaks \
    -1 species1.cds.fa \
    -2 species2.cds.fa \
    -p pairs.tsv \
    -o results/
```

`pairs.tsv` 格式（TAB 或 CSV，三列）：

```
seq1_id<TAB>seq2_id<TAB>pair_name
AT1G01010<TAB>Brara.A01.1<TAB>pair_001
```

更换计算方法 | Switch method:

```bash
biopytools kaks -1 a.fa -2 b.fa -p pairs.tsv -o results/ -m YN --verbose
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-1, --fasta1` | 第一个物种的 CDS FASTA 文件 |
| `-2, --fasta2` | 第二个物种的 CDS FASTA 文件 |
| `-p, --pairs`  | 序列配对文件（TSV/CSV：seq1_id, seq2_id, pair_name） |
| `-o, --output` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --method` | `GMYN` | 计算方法，可选：`GMYN, MYN, YN, NG, LWL, LPB, MLWL, MLPB, GY, MS, MA, GNG, GLWL, GLPB, GMLWL, GMLPB, GYN` |
| `--kaks-path` | `KaKs_Calculator` | KaKs_Calculator 可执行文件路径 |
| `-v, --verbose` | off | 启用详细日志（含命令、调试信息） |
| `--temp-dir` | 自动 | 自定义临时目录 |
| `--keep-temp` | off | 保留临时文件（用于调试） |

（运行 `biopytools kaks -h` 查看完整参数列表）

## 输出 | Output

```
results/
├── kaks_detailed.csv       # 详细结果 CSV
├── kaks_detailed.tsv       # 详细结果 TSV
├── kaks_summary.xlsx       # Excel 汇总（含统计）
├── summary_stats.json      # 统计 JSON：omega 分布、选择压力分类
├── 00_pipeline_info/       # 流程信息
└── 99_logs/                # 运行日志
```

终端会打印 omega (Ka/Ks) 的均值/中位数/范围，以及正选择 (omega>1)、中性 (omega=1)、纯化选择 (omega<1) 的数量与百分比。

## 依赖 | Dependencies

- `KaKs_Calculator 2.0`（可通过 `--kaks-path` 指定）
- `seqkit`（用于按 ID 抽取子集）
- Python: `pandas`, `openpyxl`

## 引用 | Citation

- Zhang Z. et al. KaKs_Calculator 2.0: A toolkit incorporating gamma-series methods and sliding window strategies. Genomics Proteomics Bioinformatics. 2010.
- Wang DP. et al. Gamma-MYN: a new algorithm for estimating Ka and Ks with consideration of variable substitution rates. Biology Direct. 2009. (GMYN)
- Yang Z. & Nielsen R. Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol Biol Evol. 2000. (YN)

## 相关链接 | References

- [KaKs_Calculator](https://sourceforge.net/projects/kakscalculator2/)
- [项目主页](https://github.com/lixiang117423/biopytools)
