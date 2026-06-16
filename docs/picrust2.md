# PICRUSt2微生物功能预测 | PICRUSt2 Functional Prediction

**基于16S rRNA/ASV标记基因预测微生物群落的功能基因丰度与代谢通路 | Predict functional gene abundances and metabolic pathways from marker gene sequences**

## 功能概述 | Overview

封装 PICRUSt2 流程，输入代表序列 FASTA 与对应的序列丰度表，自动识别 BIOM / TSV / Excel / Mothur 等多种丰度表格式并转换为 PICRUSt2 所需格式，完成序列系统发育位置预测、隐状态预测（HSP）、功能丰度（EC / KO / GO / PFAM / BIGG / CAZY）推断以及 MinPath 代谢通路重建。流程支持大样本自动拆分（split）与单参考（single）两种内部策略。

- 自动检测并转换 BIOM / TSV / Excel / Mojur 多种输入表格式
- 默认预测 EC + KO，可扩展到 GO / PFAM / BIGG / CAZY
- 可选分层（stratified）输出，可控制是否做 MinPath 通路重建与 gap filling
- 提供 NSTI、最小比对比例等质控参数

## 快速开始 | Quick Start

```bash
biopytools picrust2 -s study_seqs.fna -i seqabun.biom -o picrust2_out

# 使用 TSV 丰度表并生成分层输出
biopytools picrust2 -s study_seqs.fna -i seqabun.tsv -o picrust2_out --stratified

# 预测 EC、KO、GO 三类功能
biopytools picrust2 -s study_seqs.fna -i seqabun.biom -o picrust2_out --in-traits EC,KO,GO
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-s, --study-fasta` | 代表序列（ASV/OTU）FASTA 文件，未比对 |
| `-i, --input` | 序列丰度表文件（自动识别 BIOM / TSV / Excel / Mothur） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./picrust2_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--max-nsti` | `2.0` | 最大 NSTI 阈值（高于此值的 ASV 将被排除） |
| `--in-traits` | `EC,KO` | 预测的功能数据库：`EC,KO,GO,PFAM,BIGG,CAZY` |
| `--stratified` | False | 生成分层输出（按 ASV 拆分贡献） |
| `--placement-tool` | `epa-ng` | 序列放置工具：`epa-ng` / `sepp` |
| `--hsp-method` | `mp` | HSP 方法：`mp` / `emp_prob` / `pic` / `scp` / `subtree_average` |
| `--edge-exponent` | `0.5` | HSP edge exponent |
| `--pipeline` | `auto` | 流程类型：`auto` / `split` / `single` |
| `--min-align` | `0.8` | 最小比对比例 |
| `--min-reads` | `1` | 每 ASV 最小 reads 数 |
| `--min-samples` | `1` | 每 ASV 最小样本数 |
| `--no-pathways` | False | 跳过通路推断 |
| `--coverage` | False | 计算通路覆盖度 |
| `--skip-minpath` | False | 跳过 MinPath |
| `--no-gap-fill` | False | 跳过 gap filling |
| `--per-sequence-contrib` | False | 逐序列运行 MinPath |
| `--skip-norm` | False | 跳过归一化 |
| `--remove-intermediate` | False | 移除中间文件 |
| `--verbose` | False | 详细输出 |

（运行 `biopytools picrust2 -h` 查看完整参数列表）

## 输出 | Output

输出目录典型内容：

- `sequences.fasta` / `study_sequences.fasta`：标准化后的输入序列
- `ec_metagenome_out/`：EC 功能丰度预测（含分层表）
- `ko_metagenome_out/`：KO 功能丰度预测
- `pathways_out/`：MinPath 代谢通路推断结果（路径丰度、覆盖度）
- `intermediate/`：放置树、HSP 中间文件
- `NSTI` / `marker_nsti.tsv`：每 ASV 的 NSTI 质控值
- `00_pipeline_info/`、`99_logs/`：流程信息与日志

**自动后处理**：流程结束后会自动调用 PICRUSt2 官方 `add_descriptions.py` 给三类功能丰度表（通路 / EC / KO）添加描述列，并附加样本均值列、按均值降序排列。原始 PICRUSt2 输出保留为 `*_raw.tsv` 便于核对。

## 依赖 | Dependencies

- PICRUSt2 主程序（含 `picrust2_pipeline.py`、`epa-ng`、`hsp.py`、`minpath.py` 等）
- HMMER、EPA-ng（或 SEPP）、BIOM 格式工具
- Python 包：`pandas`、`biom-format`、`openpyxl`（读 xlsx）

推荐通过 `conda install -c bioconda picrust2` 安装 PICRUSt2，其会自动拉取所需依赖。

## 引用 | Citation

- Douglas GM, Maffei VJ, Zaneveld JR, et al. PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology*. 2020, 38(6):685-688.
- Langille MGI, et al. Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. *Nature Biotechnology*. 2013, 31(9):814-821.

## 相关链接 | References

- [PICRUSt2 主页](https://github.com/picrust/picrust2)
- [PICRUSt2 文档](https://picrust.github.io/picrust2/)
- [项目主页](https://github.com/lixiang117423/biopytools)
