# PGGB泛基因组图构建 | PGGB Pangenome Graph Builder

**使用 wfmash + seqwish + smoothxg 构建泛基因组变异图 | Build pangenome variation graphs using wfmash + seqwish + smoothxg.**

## 功能概述 | Overview

`pggb` 模块封装了 PGGB (Pangenome Graph Builder) 流程，将一组单倍体/组装基因组以 whole-genome alignment + 变异图的方式构建为泛基因组图（GFA/OG/GBZ），并可输出 VCF。流程内部以 `wfmash` 进行分段全基因组比对、`seqwish` 将比对转化为变异图、`smoothxg` 进行图归一化和排序，再用 `odgi`/`vg` 完成可视化和统计。

PGGB 通过 conda 环境调用，默认环境名为 `pggb_v.0.7.4`，可通过 `--conda-env` 自定义。所有命令均经 `conda run` 包装执行。

## 快速开始 | Quick Start

```bash
# 最简流程（单倍型自动识别）
biopytools pggb -i genomes.fa -o output/

# 指定单倍型数、一致度阈值，生成 VCF 并压缩
biopytools pggb -i genomes.fa -o output/ \
    -n 8 -p 95 -s 3000 \
    --vcf-spec "ref#1" --compress --stats -t 48

# 中断后续跑
biopytools pggb -i genomes.fa -o output/ --resume
```

输入 FASTA 中每条记录被视为一个单倍型/组装，染色体名建议以 `#` 分隔（PGGB 约定）。

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 FASTA 文件（多基因组已合并） |
| `-o, --output` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `24` | 线程数 |
| `--conda-env` | `pggb_v.0.7.4` | PGGB 所在 conda 环境名 |
| `-s, --segment-length` | `5000` | wfmash 比对分段长度 |
| `-p, --map-pct-id` | `90` | wfmash 比对一致度 (%, 整数) |
| `-n, --n-haplotypes` | `0` | 单倍型数；`0` 为自动检测 |
| `--vcf-spec` | - | VCF 输出参考规范，如 `ref#1` |
| `--resume` | off | 断点续传，跳过已有产物 |
| `--compress` | off | 压缩输出 (`.gz`) |
| `--stats` | off | 由 odgi 生成统计信息 |

> 完整 wfmash/smoothxg/odgi 高级参数可通过调用 `python -m biopytools.pggb.main -h` 查看，包括 `-l/--block-length`、`-c/--n-mappings`、`-K/--mash-kmer`、`--no-splits`、`--skip-normalization`、`--input-paf`、`--keep-temp` 等。

## 输出 | Output

```
output/
├── *.gfa(.gz)        # 变异图 GFA 格式
├── *.og(.gz)         # ODGI 二进制图
├── *.gbz             # GBZ 格式（用于 vg）
├── *.paf(.gz)        # wfmash 比对结果
├── *.vcf.gz          # 可选 VCF
├── *.log / *.yml     # 运行日志与配置
└── 99_logs/pggb.log  # 模块自身日志
```

运行结束会在日志中按类别汇总输出文件大小和数量。

## 依赖 | Dependencies

- conda 环境 `pggb_v.0.7.4`（或自定义），内含：
  - `wfmash`, `seqwish`, `smoothxg`, `odgi`, `vg`
- `samtools`（用于 FASTA 索引）

## 引用 | Citation

- Hickey G. et al. Pangenome graph construction from genome alignments with Minigraph-Cactus. Nature Biotechnology. 2023.
- Guarracino A. et al. PGGB: Building pangenome graphs. preprint / Nature Methods. (PGGB)
- Li H. et al. The design and construction of wfmash. (wfmash)

## 相关链接 | References

- [PGGB GitHub](https://github.com/pangenome/pggb)
- [项目主页](https://github.com/lixiang117423/biopytools)
