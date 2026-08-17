# 蛋白功能注释 | Protein Functional Annotation (func-anno)

一句话理解：**给一串蛋白序列自动「查户口」**——用 InterProScan 找结构域、用 eggNOG-mapper 查 GO/KEGG 注释，最后整理成两张标准表，直接喂给下游的富集分析（R）。

## 功能概述 | Overview { #overview }

- 端到端三阶段：InterProScan（结构域，可选）→ eggNOG-mapper（GO/KEGG 源，必需）→ 标准 GO/KEGG 表
- 输出严格 4 列 TSV：GO 表（gene/go_id/go_term/go_ontology）、KEGG 表（gene/kegg_id/kegg_term/kegg_category）
- 输入自动识别：单文件→by-step（不嵌套），目录→多样本 by-sample（每文件一子目录）
- 内置 KEGG 过滤：默认排除植物无关的人类/动物通路（癌症、激素、免疫等误挂项）
- 断点续传：IPS 与 eggNOG 按输出文件存在性跳过，建表阶段秒级、每次重跑重建
- 复用 `interproscan`/`eggnog_mapper` 模块，不重复实现

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools func-anno -i proteins.fa -o out/ -t 24
```

多样本（目录，每文件一个样本）：`biopytools func-anno -i proteins_dir/ -o out/ -t 24`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 结构域 | 蛋白上一段有特定形状/功能的小模块，像乐高积木里反复出现的标准件 |
| InterProScan | 专门扫描蛋白里有哪些「标准件」的软件 |
| GO 注释 | 用统一词汇描述「这蛋白干什么/在细胞哪里/参与什么过程」，方便跨物种比较 |
| KEGG 通路 | 一张张「代谢/信号路线图」，注释能告诉你蛋白参与哪条路 |
| eggNOG | 一个按进化关系组织的大数据库，用来给新蛋白「找亲戚」再借亲戚的注释 |
| 富集分析 | 看一群蛋白里哪些 GO/KEGG 条目比随机多，从而判断它们共同的功能方向 |
| ontology | GO 的三大分类：BP（生物过程）/MF（分子功能）/CC（细胞定位） |

## 输入 | Input { #input }

蛋白序列 FASTA。单文件即单样本；目录则每个蛋白文件当一个样本（识别扩展名 `.fa`/`.faa`/`.pep`/`.fasta`）。

```text
>gene1
MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFP
>gene2
MREIVHIQAGQCGNQIGAKFWEVISDEHGIDPTGSYH
```

## 参数说明 | Parameters { #parameters }

### 输入与布局 | Input & layout

**通俗理解|In plain words:** `-i` 给一个文件就是单样本，给一个目录就是多样本（每文件一个样本，各自建子目录）。`--by-sample` 强制单文件也建子目录，好处是往同一个 `-o` 反复跑不同样本不会互相覆盖。`-s` 只对单文件生效，用来自定义样本名/输出前缀（默认取文件名）。

相关参数：`-i/--input`（必需）、`-o/--output-dir`（必需）、`-s/--sample-name`、`--by-sample`。

### 阶段开关与复用 | Phase toggles & reuse

**通俗理解|In plain words:** IPS 只出结构域、不影响 GO/KEGG 表，赶时间可 `--skip-ips` 关掉；eggnog 是 GO/KEGG 的来源，一般别关。已经跑过一遍想只重做建表，可用 `--ips-result`/`--eggnog-result` 直接复用旧结果，跳过昂贵的比对步骤。

相关参数：`--skip-ips`、`--skip-eggnog`、`--ips-result`、`--eggnog-result`。

### KEGG 过滤 | KEGG filtering

**通俗理解|In plain words:** eggNOG 有时会把植物蛋白误挂到人类/动物的通路名（癌症、激素、免疫等）。`--kegg-exclude-keywords` 默认就用一套内置的「植物无关词」黑名单把这些洗掉；想不过滤就传空字符串。`--kegg-exclude-categories` 按大类排除（默认空，因为整块排除容易误伤植物的同源通路如 P450）。

相关参数：`--kegg-map`、`--kegg-exclude-keywords`（默认内置植物无关词）、`--kegg-exclude-categories`（默认空）。

### eggNOG 透传 | eggNOG passthrough

**通俗理解|In plain words:** `-m` 选搜索模式，mmseqs 最快是默认，diamond 次之，hmmer 最敏感最慢；`--data-dir` 指向 eggNOG 数据库目录，`--emapper-path` 覆盖 emapper.py 路径。**默认一般不用动，除非数据库装在非标准位置。**

相关参数：`-m/--mode`（默认 mmseqs）、`--data-dir`、`--emapper-path`。

## 分析流程 | Pipeline { #pipeline }

```text
蛋白 FASTA(单文件/目录)
    │
    ▼
阶段1: InterProScan 结构域(可选,可跳过/复用)
    │
    ▼
阶段2: eggNOG-mapper 注释(GO/KEGG 来源,必需)
    │
    ▼
阶段3: 建 GO/KEGG 标准表(秒级,每次重跑)
    │
    ▼
{prefix}.go.tsv + {prefix}.kegg.tsv
```

## 输出 | Output { #output }

单样本（by-step）结构：

```text
out/
├── 01_interproscan/
│   └── sample.tsv            # InterProScan 结果(TSV;另有 XML 与整理报告)
├── 02_eggnog/
│   └── 01_emapper/
│       └── sample.emapper.annotations   # eggNOG 原始注释
├── 03_tables/
│   ├── sample.go.tsv         # GO 表(gene/go_id/go_term/go_ontology)
│   └── sample.kegg.tsv       # KEGG 表(gene/kegg_id/kegg_term/kegg_category)
└── 99_logs/
    └── func_anno.log         # 运行日志
```

多样本（by-sample）时，上述结构整体嵌套在 `out/{sample}/` 下。

## 结果解读 | Interpreting Results { #interpreting-results }

- **GO 表**：一行一个「基因-条目」对。`go_ontology` 是 BP/MF/CC；同一基因可有多行。`go_term` 缺失（空）说明该 GO id 不在内置映射里，一般极少
- **KEGG 表**：一行一个「基因-通路」对，只保留 `ko` 编号（同一通路冗余的 map 编号被丢弃）。`kegg_category` 是 BRITE 大类 A
- 日志末尾会打印两表的行数和缺失统计：行数越多注释覆盖越好；KEGG 被排除的行数可用来判断过滤是否合理
- 这两张表可直接作为 clusterProfiler 等 R 富集工具的输入

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只要 GO/KEGG 富集、不关心结构域：`--skip-ips` 大幅省时
- 已跑过 eggNOG、只想调整 KEGG 过滤：`--eggnog-result` 复用旧注释，改 `--kegg-exclude-*` 重跑建表
- 大批量蛋白：保持默认 `-m mmseqs`；对敏感度要求极高才换 `-m hmmer`
- 数据库在非标准位置：`--data-dir` + `--emapper-path` 显式指定

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 蛋白序列 FASTA(单文件→by-step) 或目录(多样本→by-sample)｜Protein FASTA (single→by-step) or dir (multi→by-sample) |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `-s, --sample-name` | — |  | 样本名/前缀(仅单文件模式生效)｜Sample name (single-file mode only) |
| `--by-sample` | — |  | 强制 by-sample(单文件也建 sample 子目录, 往同一 -o 多次跑不覆盖)｜Force by-sample layout |
| `--ips-result` | — |  | 复用已有 IPS 结果目录(跳过 IPS)｜Reuse existing IPS dir |
| `--eggnog-result` | — |  | 复用已有 .emapper.annotations(跳过 eggnog)｜Reuse existing annotations |
| `--skip-ips` | — |  | 跳过 IPS(只要 GO/KEGG)｜Skip IPS |
| `--skip-eggnog` | — |  | 跳过 eggnog｜Skip eggnog |
| `--kegg-map` | — |  | 外部 KEGG 映射 TSV(补 category)｜External KEGG map (fill category) |
| `--kegg-exclude-keywords` | — |  | KEGG 通路 name 黑名单(逗号分隔子串, None=内置植物无关词 cancer/estrogen 等)｜KEGG name blacklist (None=built-in plant-irrelevant) |
| `--kegg-exclude-categories` | `` |  | 排除的 KEGG 分类(逗号分隔, 匹配 category A/B 子串, 如 Human Diseases)｜Exclude KEGG categories (substring match) |
| `-m, --mode` | `mmseqs` | mmseqs/diamond/hmmer | eggnog 搜索模式｜eggnog search mode |
| `--data-dir` | — |  | eggnog DB 目录｜eggnog DB dir |
| `--emapper-path` | — |  | emapper.py 路径覆盖｜emapper.py path override |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 蛋白序列 FASTA(单文件→by-step) 或目录(多样本→by-sample)｜Protein FASTA (single→by-step) or dir (multi→by-sample) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `-t, --threads` | `12` | int | 线程数｜Threads (default 12) |
| `-s, --sample-name` | — |  | 样本名/输出前缀(默认输入文件名, 仅单文件模式生效)｜Sample name (single-file mode only) |
| `--by-sample` | — | store_true | 强制 by-sample(单文件也建 sample 子目录, 往同一 -o 多次跑不覆盖)｜Force by-sample layout |
| `--ips-result` | — |  | 复用已有 IPS 结果目录(跳过 IPS)｜Reuse existing IPS dir |
| `--eggnog-result` | — |  | 复用已有 .emapper.annotations(跳过 eggnog)｜Reuse existing annotations |
| `--skip-ips` | — | store_true | 跳过 IPS(只要 GO/KEGG)｜Skip IPS |
| `--skip-eggnog` | — | store_true | 跳过 eggnog(仅整理已有结果)｜Skip eggnog |
| `--kegg-map` | — |  | 外部 KEGG 映射 TSV(ko_id\tname\tcategory, 补 category)｜External KEGG map to fill category |
| `--kegg-exclude-keywords` | — |  | KEGG 通路 name 黑名单(逗号分隔子串, None=内置植物无关词 cancer/estrogen 等)｜KEGG name blacklist (None=built-in) |
| `-m, --mode` | `mmseqs` | mmseqs/diamond/hmmer | 搜索模式｜Search mode |
| `--data-dir` | — |  | eggnog DB 目录｜eggnog DB dir |
| `--emapper-path` | — |  | emapper.py 路径覆盖｜emapper.py path override |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- InterProScan（`interproscan` 模块调用，默认 `~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh`）
- eggNOG-mapper（`eggnog_mapper` 模块调用，emapper.py 默认 `~/miniforge3/envs/annot/bin/emapper.py`，conda 环境 `annot`）
- eggNOG 数据库（默认 `~/database/eggnog`，环境变量 `EGGNOG_DATA_DIR` 可覆盖）
- 内置 GO 映射（`interproscan/go_data.py`，无需额外下载）

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
会。IPS 按 `01_interproscan/{sample}.tsv`、eggnog 按 `02_eggnog/01_emapper/{sample}.emapper.annotations` 存在性跳过；建表阶段每次重跑重建（保证过滤参数生效）。

**Q2：`--eggnog-result` 报「结果不存在」？**
该参数要求路径必须存在，不存在会直接报错（不是警告后重跑）。传一个真实的 `.emapper.annotations` 文件路径。

**Q3：KEGG 表里为什么没有某些预期通路？**
先看是不是被内置黑名单误杀了（癌症/激素/免疫等词），需要保留就用 `--kegg-exclude-keywords` 传空字符串关闭过滤；也可能是 eggNOG 没给出 `ko` 编号（非 `ko` 开头的被丢弃）。

**Q4：IPS 失败会不会影响 GO/KEGG 表？**
不会。IPS 只是结构域注释，失败/无 TSV 只警告不阻断，GO/KEGG 表照常生成。
