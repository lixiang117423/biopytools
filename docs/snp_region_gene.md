# SNP 区域基因提取 | SNP Region Gene Extractor

一句话理解：**给你一批 SNP 位点，自动找出每个位点附近（含启动子）的基因，并把它们的 CDS 和蛋白序列批量抽出来**，省去逐个手动翻注释找基因的麻烦。

## 功能概述 | Overview { #overview }

- 输入 SNP 位置列表 + GFF3 注释 + 参考基因组，输出相关基因的 CDS 和蛋白序列
- 支持 SNP 上游、下游距离和启动子区域三个可调窗口
- 自动判断 SNP 落在基因的哪个部位：启动子 / 外显子 / 内含子 / 基因间区
- 输出配套基因列表，逐 SNP 记录命中的基因、转录本、链方向、部位、距离
- 重复 SNP 自动去重，临时文件默认清理

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools snp-region-gene -i snp.txt -g annotation.gff3 -G genome.fa -l 100000 -r 100000 -o output
```

最小输入：SNP 文件 + GFF3 注释 + 基因组 FASTA 三个文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| SNP | 基因组上的单碱基差异位点，像马路上的一个「点位」 |
| GFF3 | 基因注释格式，记录每个基因/外显子在哪条染色体、从几号到几号 |
| CDS | 编码序列，真正翻译成蛋白的那段 DNA |
| 启动子 | 基因前面的「开关区」，调控该基因表达，通常取基因上游约 2000 bp |
| 上游/下游 | 相对基因方向而言；对正链基因，上游=坐标更小的一侧 |
| 链方向 | 基因在 DNA 双链中的哪一条上（+ 或 -），决定「上游」实际在坐标的哪一边 |

## 输入 | Input { #input }

### SNP 位置文件

纯文本，一行一个位点，格式 `染色体:位置`：

```text
Chr01:24770
Chr01:128900
Chr02:55555
```

- 空行、`#` 开头的注释行、含 SNP/POSITION/COORDINATE 关键词的表头行会自动跳过
- 重复位点自动去重

### GFF3 注释文件

标准 GFF3 格式，需包含 gene、mRNA（或 transcript）、exon 三类特征，属性用 `ID=`/`Parent=` 关联：

```text
Chr01	source	gene	1000	5000	.	+	.	ID=Gene001
Chr01	source	mRNA	1000	5000	.	+	.	ID=Gene001.1;Parent=Gene001
Chr01	source	exon	1000	1500	.	+	.	Parent=Gene001.1
```

### 基因组 FASTA

参考基因组序列文件，染色体名需与 GFF3 一致（gffread 提取序列时按染色体名匹配）。

## 参数说明 | Parameters { #parameters }

### 必需输入 | Required

**通俗理解|In plain words:** 三个输入缺一不可：SNP 列表告诉工具「找哪几个点」，GFF3 告诉它「基因在哪」，基因组 FASTA 供它「把序列抽出来」。

相关参数：`-i/--snp`、`-g/--gff`、`-G/--genome`。

### 查询窗口 | Search window

**通俗理解|In plain words:** 决定「离 SNP 多远算有关」：`-l`/`-r` 分别是在 SNP 左右各扩多大范围去找基因（默认都是 0，即只在位点本身所在基因里找），调大能抓到更远的基因，但可能带进不相关的基因；`-p` 是启动子长度（默认 2000 bp），它同时决定「SNP 是否算落在某基因的启动子上」以及基因的有效范围是否向外扩出启动子区。**一般 `-p` 用默认 2000 即可；`-l`/`-r` 按研究需要（如找 GWAS 显著位点上下游 100 kb 内的基因就设 100000）。**

相关参数：`-l/--left`（默认 0）、`-r/--right`（默认 0）、`-p/--promoter`（默认 2000）。

### 输出与工具路径 | Output & tool paths

**通俗理解|In plain words:** `-o` 是输出文件前缀（不是目录，三个结果文件都挂在这个前缀后面）；`--gffread-path`/`--seqkit-path` 只在功能域环境缺失、系统 PATH 里也找不到这两个软件、或想指定特定版本时才需要改；`--keep-temp` 用于排查问题时保留中间临时文件，平时不用开。

相关参数：`-o/--output`（默认 `./snp_region_output`）、`--gffread-path`（默认 `gffread`）、`--seqkit-path`（默认 `seqkit`）、`--keep-temp`。

## 分析流程 | Pipeline { #pipeline }

```text
SNP 文件 + GFF3 + 基因组 FASTA
    |
    v
解析 SNP 列表(去重) + 建立基因索引
    |
    v
gffread 一次性提取全基因组 CDS + 蛋白(临时文件)
    |
    v
逐个 SNP 查找相关基因(含启动子范围)
    |
    v
seqkit 按命中转录本 ID 过滤出目标序列
    |
    v
写 CDS FASTA + 蛋白 FASTA + 基因列表,清理临时文件
```

## 输出 | Output { #output }

```text
output_cds.fasta          # 命中基因的 CDS 序列
output_protein.fasta      # 命中基因的蛋白序列
output_gene_list.txt      # SNP-基因对应关系表
output.log                # 运行日志
```

（`output` 替换为你的 `-o` 前缀；临时文件 `output_temp_all_cds.fasta`、`output_temp_all_protein.fasta` 默认删除，`--keep-temp` 时保留。）

`output_gene_list.txt` 列：

| 列名 | 含义 |
|---|---|
| SNP_Position | SNP 标识（如 Chr01:24770） |
| Chromosome / SNP_Pos | SNP 所在染色体与坐标 |
| Gene_ID / mRNA_ID | 命中的基因与转录本 ID |
| Strand | 基因链方向（+/-） |
| Feature | SNP 落在基因的部位（promoter/exon/intron/intergenic） |
| Distance | SNP 到基因的距离（bp） |

## 结果解读 | Interpreting Results { #interpreting-results }

- `_cds.fasta` / `_protein.fasta`：命中的基因序列，FASTA 头为转录本 ID，可直接做下游分析
- `_gene_list.txt`：Feature 列是最关键的判据——`exon`=落在外显子，`intron`=落在内含子，`promoter`=落在启动子，`intergenic`=基因间区
- 若某个 SNP 没有任何命中，说明其 `-l`/`-r` 窗口内没有基因（含启动子范围也算）
- 命中数偏少时，先检查 SNP 文件的染色体名是否与 GFF3 完全一致（大小写、前缀都要一致）

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只关心位点落在哪个基因：`-l 0 -r 0`（默认），再结合 `-p` 判断是否在启动子
- 找 GWAS 显著位点附近的候选基因：`-l 100000 -r 100000`（上下游 100 kb）或更大
- 只关心基因本体不关心启动子：`-p 0`
- 排查「提取序列为空」：加 `--keep-temp` 看中间临时文件，并用日志确认 gffread/seqkit 是否成功

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--snp, -i` | 必填 |  | SNP位置文件（格式：Chr01:24770）｜SNP position file (format: Chr01:24770) |
| `--gff, -g` | 必填 |  | GFF3注释文件｜GFF3 annotation file |
| `--genome, -G` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--left, -l` | `0` | int | SNP上游距离（bp）｜Upstream distance from SNP (bp) |
| `--right, -r` | `0` | int | SNP下游距离（bp）｜Downstream distance from SNP (bp) |
| `--promoter, -p` | `2000` | int | 启动子区域距离（bp）｜Promoter region distance (bp) |
| `--output, -o` | `./snp_region_output` |  | 输出文件前缀｜Output file prefix |
| `--gffread-path` | `gffread` |  | gffread程序路径｜gffread program path |
| `--seqkit-path` | `seqkit` |  | seqkit程序路径｜seqkit program path |
| `--keep-temp` | — |  | 保留临时文件｜Keep temporary files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --snp` | 必填 |  | SNP位置文件（格式：Chr01:24770）｜SNP position file (format: Chr01:24770) |
| `-g, --gff` | 必填 |  | GFF3注释文件｜GFF3 annotation file |
| `-G, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-l, --left` | `0` | int | SNP上游距离（bp）｜Upstream distance from SNP (bp) |
| `-r, --right` | `0` | int | SNP下游距离（bp）｜Downstream distance from SNP (bp) |
| `-p, --promoter` | `2000` | int | 启动子区域距离（bp）｜Promoter region distance (bp) |
| `-o, --output` | `./snp_region_output` |  | 输出文件前缀｜Output file prefix |
| `--gffread-path` | — |  | gffread可执行文件路径(默认域环境自动解析)｜gffread executable path (default: auto domain env) |
| `--seqkit-path` | — |  | seqkit可执行文件路径(默认域环境自动解析)｜seqkit executable path (default: auto domain env) |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep temporary files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- gffread（自动解析 rna 域环境并经 conda run 调用，可用 `--gffread-path` 或环境变量 GFFREAD_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- seqkit（自动解析 misc 域环境并经 conda run 调用，可用 `--seqkit-path` 或环境变量 SEQKIT_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- Python 3 标准库

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
不会。每次运行都重新用 gffread 提取全基因组序列，没有「已完成则跳过」的判断，换参数重跑直接重新运行即可。

**Q2：为什么有的 SNP 在基因列表里没有记录？**
该 SNP 的 `-l`/`-r` 窗口与任何基因的有效范围（含启动子）都不重叠，或染色体名与 GFF3 不一致。

**Q3：GFF 用了 GTF 格式（`gene_id "xxx"`）能识别吗？**
不能。本工具只认 `=`（GFF3 风格），GTF 的空格/引号风格会导致 ID/Parent 解析失败、基因关联不上。

**Q4：提取出来的序列为什么可能对不上？**
gffread 按转录本 ID 从全基因组提取序列，再用 seqkit 按同一 ID 过滤；若 GFF3 里的 mRNA ID 与 gffread 输出 FASTA 头里的转录本 ID 不一致（常见于注释不规范），会过滤不到或过滤错。
