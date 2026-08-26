# NCBI 基因组批量下载 | ncbi-datasets

输入一个 NCBI taxon 编号(物种分类单元 ID),就能把该物种(含属下全部物种)在 NCBI 上公布的所有基因组序列、注释、蛋白序列一次性批量下载回来,不用自己逐个 accession 手动点。

## 功能概述 | Overview { #overview }

- 基于 NCBI 官方 `datasets` CLI,下载过程与 NCBI 数据库完全一致,不需要自己拼 FTP 链接
- 输入 taxon 编号,自动查询该 taxon 下所有 assembly,先生成清单(数量、accession、组装级别)再下载
- 支持常用筛选:RefSeq/GenBank 来源、组装级别(complete/chromosome/scaffold/contig)、只要参考基因组、只要有注释
- 默认只下载基因组序列(`genome`);GFF3 注释、蛋白序列、CDS、seq-report 按需加参数
- 断点续传:zip 已存在跳过下载、解压目录已存在跳过解压,重跑不重复下载
- `datasets` 工具缺失时自动从 NCBI 官方地址下载安装到 `~/bin/`,无需 conda、无需 root

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools ncbi-datasets -t 67593 -o ~/tmp/ncbi_download/
```

下载 Phytophthora sojae(taxon 67593)所有基因组序列到 `~/tmp/ncbi_download/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| taxon 编号 | NCBI 给每个物种(或属、科等分类单元)分配的身份证号,如 67593 = 大豆疫霉 *Phytophthora sojae*;填属的编号会包含属下所有种 |
| assembly | 一次基因组组装结果,NCBI 上同一物种可能有多个人提交的多份组装 |
| accession | 每份 assembly 的唯一编号,如 `GCF_000149865.2`(GCA 开头=GenBank 源,GCF 开头=RefSeq 源) |
| RefSeq / GenBank | NCBI 两套数据库:RefSeq 是官方整理去重版(质量更整齐),GenBank 是原始提交版(数量更多) |
| 组装级别 | 组装完成度:complete > chromosome > scaffold > contig,级别越高越完整 |

## 输入 | Input { #input }

只需一个 **NCBI taxon 编号**(正整数)。查编号方法:

- NCBI Taxonomy 网站搜物种名,URL 里的 `taxid` 就是编号
- 属级编号(如 *Phytophthora* = 4783)会包含属下所有种的 assembly

下载前建议先 `--dry-run` 看清单,避免误下载超大 taxon。

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** 所有筛选参数都是「先看清单再决定要不要下」——先用 `--dry-run` 试出清单大小,觉得合适再去掉 `--dry-run` 真正下载。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--taxon, -t` | 必填 | int | NCBI taxon 编号｜NCBI taxon ID |
| `--output-dir, -o` | `./output` |  | 输出目录｜Output directory |
| `--assembly-source` | — | refseq/genbank | 只下载 RefSeq 或 GenBank 来源｜Only RefSeq or GenBank assemblies |
| `--assembly-level` | — |  | 组装级别过滤(逗号分隔)｜Assembly level filter (comma-separated): complete, chromosome, scaffold, contig |
| `--reference` | — |  | 只下载参考基因组｜Reference genomes only |
| `--annotated` | — |  | 只下载有注释的基因组｜Annotated genomes only |
| `--include-gff3` | — |  | 额外下载 GFF3 基因注释｜Also download GFF3 gene annotation |
| `--include-protein` | — |  | 额外下载蛋白序列｜Also download protein sequences |
| `--include-cds` | — |  | 额外下载 CDS 序列｜Also download CDS sequences |
| `--include-seq-report` | — |  | 额外下载 seq-report 汇总｜Also download seq-report summary |
| `--dry-run` | — |  | 只查询 assembly 清单,不下载｜Query manifest only, no download |
| `--organize/--no-organize` | `True` |  | 下载后整理到 02_organized｜Organize into 02_organized after download |
| `--datasets-path` | — |  | datasets 工具路径(默认走环境变量/配置/~bin)｜datasets tool path (default: env var / config / ~/bin) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--taxon, -t` | 必填 |  | NCBI taxon 编号｜NCBI taxon ID |
| `--output-dir, -o` | `./output` |  | 输出目录｜Output directory |
| `--assembly-source` | — | refseq/genbank | 只下载 RefSeq 或 GenBank 来源｜Only RefSeq or GenBank assemblies |
| `--assembly-level` | — |  | 组装级别过滤(逗号分隔)｜Assembly level filter (comma-separated): complete, chromosome, scaffold, contig |
| `--reference` | — | store_true | 只下载参考基因组｜Reference genomes only |
| `--annotated` | — | store_true | 只下载有注释的基因组｜Annotated genomes only |
| `--include-gff3` | — | store_true | 额外下载 GFF3 基因注释｜Also download GFF3 gene annotation |
| `--include-protein` | — | store_true | 额外下载蛋白序列｜Also download protein sequences |
| `--include-cds` | — | store_true | 额外下载 CDS 序列｜Also download CDS sequences |
| `--include-seq-report` | — | store_true | 额外下载 seq-report 汇总｜Also download seq-report summary |
| `--dry-run` | — | store_true | 只查询 assembly 清单,不下载｜Query manifest only, no download |
| `--no-organize` | — | store_true | 下载后不整理到 02_organized｜Do not organize into 02_organized after download |
| `--datasets-path` | — |  | datasets 工具路径(默认走环境变量/配置/~bin)｜datasets tool path (default: env var / config / ~/bin) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

## 输出 | Output { #output }

```
output/
├── 00_pipeline_info/
│   ├── 67593.assemblies.tsv      # assembly 清单(下载前先看这个)
│   └── software_versions.yml     # datasets 与模块版本
├── 01_download/                  # 原始下载(NCBI 打包结构,不动)
│   ├── 67593.genomes.zip         # 官方打包的下载文件
│   └── 67593.ncbi_dataset/
│       └── data/
│           ├── GCF_000149865.2/
│           │   ├── GCF_000149865.2_genomic.fna
│           │   ├── GCF_000149865.2_genomic.gff
│           │   └── protein.faa (加 --include-protein 时)
│           └── assembly_data_report.jsonl
├── 02_organized/                 # 整理产物(软链,统一命名,下游直接用)
│   ├── genomes/                  #   GCF_000149865.2.fa
│   ├── gff/                      #   GCF_000149865.2.gff (加 --include-gff3 时)
│   ├── protein/                  #   GCF_000149865.2.faa (加 --include-protein 时)
│   ├── cds/                      #   GCF_000149865.2.cds.fa (加 --include-cds 时)
│   └── files.tsv                 # 索引:accession/type/绝对路径,可直接喂下游 fof
└── 99_logs/
    └── 67593.ncbi_datasets.log   # 全量日志(+ .out.log / .err.log 分离)
```

`02_organized/` 下所有序列文件统一后缀:基因组 `.fa`(不保留 NCBI 的 `_genomic.fna` 名)、蛋白 `.faa`、CDS `.cds.fa`;均为软链,不占额外空间。`files.tsv` 三列 `accession/type/path`,绝对路径,下游(fastani、busco、fof 清单)拿来即用。不想要整理产物用 `--no-organize`。

## 结果解读 | Interpreting Results { #interpreting-results }

- **67593.assemblies.tsv**:每行一个 assembly,列 = accession、物种名、组装级别、状态。`assembly_status` 为 `current` 表示当前有效;下载前先数行数,行数=将下载的基因组数
- **zip 体积判断**:基因组数量 × 单基因组大小粗估;加了 `--include-gff3/--include-protein` 会明显变大
- **质量判断**:优先看清单里 `assembly_level` 是否 `Complete Genome`/`Chromosome`;`GCF_` 前缀(RefSeq)比 `GCA_`(GenBank)通常更权威
- **数据完整性**:`01_download/67593.ncbi_dataset/data/` 下 accession 子目录数应等于清单行数;对不上说明下载被中断,重跑同命令会自动续传
- **`02_organized/files.tsv`**:整理步骤的产物索引——`genomes/` 下 `.fa` 文件数应等于清单行数;缺注释/蛋白的 accession 不会出现在对应类型目录(已选 `--include-gff3` 等却一个都没有时会有 WARNING)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 推荐参数<br>Recommended |
|------|------|
| 泛基因组/比较基因组(只要序列) | 默认,不加 include 参数 |
| 下游要基因结构/注释 | 加 `--include-gff3 --include-protein` |
| 只要高质量参考基因组 | `--assembly-source refseq --assembly-level complete,chromosome` |
| 数据库全量备份 | `--include-gff3 --include-protein --include-cds --include-seq-report` |
| 不确定有多少数据 | 先 `--dry-run` 看清单,再决定 |

**通俗理解|In plain words:** 先小后大——`--dry-run` 白嫖清单,筛选条件让 NCBI 少下几十 G 完全靠这几行参数。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| `datasets` CLI | 必需,未安装时模块自动下载到 `~/bin/datasets`(Linux x86_64) |
| `curl` | 自动安装 datasets 时使用 |
| `unzip` | 解压下载包时使用 |
| 网络 | 需能访问 `ftp.ncbi.nlm.nih.gov`(超算登录节点通常可出网) |

## 常见问题 | FAQ { #faq }

- **下载提示没有 assembly?** 筛选条件过严(如 taxon 下没有 RefSeq complete 基因组),去掉筛选再试;也可能 taxon 编号本身不存在,先用 NCBI Taxonomy 确认
- **超大 taxon 下载太慢/超时?** 大下载不要在登录节点跑,用作业调度系统提交(`sub` 等);模块本身支持断点续传——下载被中断留下的损坏 zip 会被自动检测、删除并重新下载,重跑同命令即可
- **`datasets` 装到哪里了?** 默认 `~/bin/datasets`,也可用 `DATASETS_PATH` 环境变量或 `~/.config/biopytools/config.yml` 的 `tools.datasets` 指定
- **zip 里没有蛋白文件?** 默认只下 `genome`,需要加 `--include-protein --include-gff3` 等参数
- **和 `ena_downloader` 的区别?** `ena_downloader` 从 ENA 下载测序原始数据(FASTQ);本模块从 NCBI 下载已组装基因组(FASTA/注释)
