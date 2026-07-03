# a-liner 共线性可视化 Pipeline

**FASTA → minimap2 → a-liner 一键共线性可视化 | One-command synteny visualization: FASTA → minimap2 → a-liner**

## 功能概述 | Overview

aliner 模块把两条（或配对多条）基因组 FASTA 序列的共线性可视化串成一条 pipeline：先用 **minimap2**（`asm5/asm10/asm20` 预设）做比对生成 PAF，再喂给 [a-liner](https://github.com/thackerrvik/a-liner) 绘制出版级共线性图（PDF）。支持按区段（`chrZ` 或 `chrZ:1-30000000`）配对指定 ref/query 侧序列、断点续传（PAF 已存在则跳过 minimap2）、版本/参数记录。适合近缘组装相互比对、挂载前后对比、特定染色体区段共线性检查。

## 快速开始 | Quick Start

```bash
# 两条染色体的共线性图
biopytools aliner \
    --ref ref.fa --query query.fa \
    --ref-seqs chr1 --query-seq chr1 \
    -o output_dir/

# 指定区段 + 远缘预设 + 配色
biopytools aliner \
    --ref ref.fa --query query.fa \
    --ref-seqs chrZ:1-30000000 --query-seqs chrZ:1-30000000 \
    --preset asm20 --colormap 3 \
    -o output_dir/ --out-prefix sample1
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `--ref` | 参考基因组 FASTA |
| `--query` | 查询基因组 FASTA |
| `--ref-seqs` | ref 侧序列（逗号分隔，如 `chr1,chr2` 或 `chrZ:1-30000000`） |
| `--query-seqs` | query 侧序列（逗号分隔，与 ref 等长，按序配对） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./aliner_output` | 输出目录 |
| `--out-prefix` | `synteny` | 输出文件前缀 |
| `--preset` | `asm5` | minimap2 预设（近缘 `asm5` / 远缘 `asm10,asm20`） |
| `--min-identity` | `70` | a-liner identity 阈值（%） |
| `--min-alignment-len` | `1000` | 最小比对长度（bp） |
| `--colormap` | `5` | 配色（0–5） |
| `--figure-size` | `6 0` | 图尺寸 [宽 高]，高 `0` 自适应 |
| `-t, --threads` | `12` | 线程数 |
| `--extra-args` | 空 | 透传给 a-liner 的额外参数 |

（运行 `biopytools aliner -h` 查看完整参数列表）

### 工具路径配置

- minimap2 / samtools 路径按以下优先级查找：环境变量（`MINIMAP2_PATH`/`SAMTOOLS_PATH`）→ `~/.config/biopytools/config.yml` → 默认 `~/miniforge3/envs/telocomp/bin/`
- a-liner 固定使用 conda 环境 `a-liner`，conda 调用自动加 `--no-capture-output`（§13.2.0）

## 输出 | Output

```
output_dir/
├── 00_pipeline_info/
│   ├── software_versions.yml        # minimap2/samtools/a-liner 版本
│   └── pipeline_params.yaml         # 运行参数
├── 01_alignment/
│   ├── {prefix}.paf                 # minimap2 比对结果
│   └── sequence_config.txt          # ref/query 序列配置
├── 02_aliner/
│   └── {prefix}.pdf                 # 共线性图 ⭐
└── 99_logs/
    └── aliner_pipeline.log          # 运行日志
```

> `--ref-seqs` 与 `--query-seqs` 必须等长并按序配对；区段格式 `name` 或 `name:start-end`。

## 依赖 | Dependencies

- **minimap2**（默认 `~/miniforge3/envs/telocomp/bin/minimap2`）
- **samtools**（默认同环境，用于取序列长度）
- **a-liner**（conda 环境 `a-liner`）
- 推荐通过 conda 环境隔离安装

## 引用 | Citation

- Li, H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 2018, 34(18): 3094–3100.
- Li, H. et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 2009, 25(16): 2078–2079.
- a-liner：见 [项目仓库](https://github.com/thackerrvik/a-liner)

## 相关链接 | References

- [a-liner](https://github.com/thackerrvik/a-liner) | [minimap2](https://github.com/lh3/minimap2)
- [项目主页](https://github.com/lixiang117423/biopytools)
