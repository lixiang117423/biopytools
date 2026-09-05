# 从测序数据直接建树 | reads2tree

输入一个装着测序原始数据(fastq)的文件夹,自动识别每对双端文件,不需要组装基因组、不需要做序列比对,直接得到物种进化树——适合大批量样本快速建树,基因组测序覆盖度 1.5X 就够用。

## 功能概述 | Overview { #overview }

- 基于 ASTER 套件的 **WASTER**(Without-Alignment/Assembly Species Tree EstimatoR),从原始 reads 直接推断物种树(先 call SNP,再 CASTER 建树)
- 输入文件夹自动识别双端测序(`_R1/_R2`、`_1/_2`、`.R1.`、`read1/read2`、`_R1_001` 等命名),同一样本多 lane 自动合并
- 支持单端 fastq 与双端混用;`.gz` 压缩自动解压
- 双端默认 `cat` 拼接(WASTER 官方推荐:不重叠 reads 效果更好);重叠双端可用 `--merge` 走 BBMerge 合并
- 物种级支持度(local bootstrap)直接标注在树上;可选追加枝长计算(`--branch-length`)
- 可配个体→物种映射(`--samples-map`),同物种多个体自动合并信号
- 断点续传:合并/建树/枝长步骤完成即跳过,重跑不重复计算

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools reads2tree -i ~/tmp/reads_dir/ -o ~/tmp/reads2tree_out/
```

`reads_dir/` 下每个样本一对双端文件(如 `A_R1.fq.gz`、`A_R2.fq.gz`),输出物种树到 `reads2tree_out/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| reads | 测序仪吐出来的短序列碎片(FASTQ 格式),带质量值;是建树的最原始材料 |
| 双端测序 | 一条 DNA 片段两头都测,得到 R1/R2 两个文件;拼接后信息更完整 |
| WASTER | ASTER 套件里的"免组装免比对"建树工具:直接从 reads 里数 k-mer、call SNP,跳过拼基因组和比对两大步 |
| local bootstrap | 树上每个分支的置信度(>95.0 表示很稳),和常见 bootstrap 类似但更快 |
| 外群 | 和所有研究样本亲缘最远的一个物种,用来给树定根(方向) |

## 输入 | Input { #input }

一个 fastq 目录,自动识别(**支持多级子目录**,如 `raw/`、`clean/`):

- **双端**:`A_R1.fq.gz` + `A_R2.fq.gz` → 样本 A(支持 `_R1/_R2`、`_1/_2`、`.R1.`、`read1/read2`、`_R1_001` 等;多 lane 自动归组)
- **质控后命名**:`SRR123_1_clean.fq.gz` + `SRR123_2_clean.fq.gz` → 样本 SRR123(fastp 等清洗产物直接可用)
- **单端**:`C.fq`、`C.fastq.gz` 等,直接作为样本 C
- 缺一侧的双端(只有 R1 没有 R2)→ 按单端处理
- 非 fastq 文件(如脚本、txt)→ 忽略并 WARNING
- **同一样本出现在多个子目录**(如 `raw/` 与 `clean/` 都有 SRR123)→ 报错并提示,避免原始+质控数据混拼;此时把 `-i` 指向具体数据目录(如 `clean/`)

> 注意:WASTER 要求输入未压缩,模块会自动解压/拼接;重叠双端 reads(插入片段短于读长)建议 `--merge` 先 BBMerge 合并,否则建树精度下降。

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** 绝大多数参数和 genome2tree 一致——`-t` 线程数、`--root` 外群定根、`--branch-length` 要不要枝长。新面孔只有 `--merge`:你的双端 reads 如果重叠(测序片段两头都覆盖了同一段),加上它让 BBMerge 先合并,树会更准。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq 目录(自动识别双端,可.gz)｜FASTQ dir (auto-detect paired-end, .gz ok) |
| `-o, --output-dir` | `./reads2tree_output` |  | 输出目录(默认./reads2tree_output)｜Output directory |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--root` | `` |  | 外群物种名(出有根树)｜Outgroup species name (rooted tree) |
| `--branch-length` | `False` |  | 追加 waster_branchlength 枝长计算｜Also compute branch lengths |
| `--samples-map` | `` |  | 个体→物种映射文件(个体名<TAB>物种名)｜individual-to-species map |
| `--merge` | `False` |  | 重叠双端 reads 用 BBMerge 合并(默认 cat 拼接)｜BBMerge overlapping paired reads (default: cat) |
| `--waster-path` | — |  | waster 路径(默认~/software/ASTER/bin/waster)｜waster binary path |
| `--bbmerge-path` | — |  | bbmerge.sh 路径(--merge 时用)｜bbmerge.sh path (for --merge) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq 目录(自动识别双端 _R1/_R2、_1/_2 等,可.gz)｜FASTQ dir (auto-detect paired-end _R1/_R2, _1/_2, .gz ok) |
| `-o, --output-dir` | `./reads2tree_output` |  | 输出目录(默认./reads2tree_output)｜Output directory |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--root` | `` |  | 外群物种名(出有根树)｜Outgroup species name |
| `--branch-length` | — | store_true | 追加 waster_branchlength 枝长计算｜Also compute branch lengths on the fixed topology |
| `--samples-map` | `` |  | 个体→物种映射文件(个体名<TAB>物种名)｜individual-to-species map |
| `--merge` | — | store_true | 重叠双端 reads 用 BBMerge 合并(默认 cat 拼接)｜BBMerge overlapping paired reads (default: cat) |
| `--waster-path` | — |  | waster 路径(默认~/software/ASTER/bin/waster)｜waster path |
| `--bbmerge-path` | — |  | bbmerge.sh 路径(--merge 时用)｜bbmerge.sh path (for --merge) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline { #pipeline }

```
reads目录
  └─ 自动识别双端(_R1/_R2 等)→ 配对分组
01_input/uncompressed/{sample}.fq   ← 解压 + R1+R2 拼接(cat / BBMerge)
01_input/input.tsv                  ← 样本<TAB>文件路径(WASTER 输入清单)
02_waster/waster_species_tree.nw   ← WASTER 物种树(local bootstrap 支持度)
03_branch_length/                   ← (可选)固定拓扑打分出枝长
```

## 输出 | Output { #output }

```
output/
├── 00_pipeline_info/
│   └── software_versions.yml       # waster 版本/路径/参数/运行时间
├── 01_input/
│   ├── input.tsv                   # 样本→合并后 fastq 路径(可人工检查)
│   ├── samples_map.tsv             # (加 --samples-map 时)
│   └── uncompressed/
│       ├── A.fq                    # A_R1.fq.gz + A_R2.fq.gz 拼接结果
│       └── B.fq
├── 02_waster/
│   ├── waster_species_tree.nw     # 物种树(Newick,含 local bootstrap)
│   └── waster.log                  # waster 运行日志
├── 03_branch_length/               # (加 --branch-length 时)
│   ├── waster_branchlength_species_tree.nw
│   └── waster_branchlength.log
└── 99_logs/
    └── reads2tree.log              # 全量日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **`waster_species_tree.nw`**:Newick 格式物种树。内部节点冒号后是 local bootstrap 支持率,**>95.0 表示该分支很可靠**,50 以下的支要谨慎解读
- **`01_input/input.tsv`**:每行 `样本名<TAB>合并后fastq路径`。先看这文件确认每个样本的 reads 都合并对了(文件大小与预期一致)
- **`02_waster/waster.log`**:waster 运行详情,报错时从这里看
- **质量判断**:样本数 ≥4 才有四聚体支持度统计;物种太少时日志会 WARNING。reads 质量差(接头污染、低质量)会拉低支持度

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 推荐参数<br>Recommended |
|------|------|
| 常规批量建树(不重叠双端) | 默认即可,`-t 64`(ASTER 16 核以上超线性加速) |
| 重叠双端 reads | 加 `--merge` |
| 需要定根 | `--root 外群样本名` |
| 需要枝长信息 | 加 `--branch-length`(同一拓扑上补枝长) |
| 同物种多个体 | `--samples-map` 映射文件 |

**通俗理解|In plain words:** ASTER 并行效率很高,核多给足(`-t 64`);别的参数按需加,默认配置对大多数情况已经够用。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| ASTER `waster` | 必需,默认 `~/software/ASTER/bin/waster`(可用 `--waster-path`/`WASTER_PATH` 指定) |
| BBMerge(`bbmerge.sh`) | 仅 `--merge` 时需要,默认 `~/miniforge3/envs/bbmap_v.39.81/bin/bbmerge.sh` |
| 内存 | **WASTER 建树需 >32GB 内存(通常 <64GB),务必提交到计算节点运行** |

## 常见问题 | FAQ { #faq }

- **登录节点能跑吗?** 不能。WASTER 硬编码 K=9,内存需求 32–64GB,必须 `sub` 提交到计算节点
- **双端 reads 需要先合并吗?** 不重叠:不需要,模块自动 `cat` 拼接(官方推荐);重叠:加 `--merge` 用 BBMerge 合并
- **和 genome2tree 的区别?** genome2tree 输入是已组装的基因组 fasta;本模块输入是测序原始 reads(fastq),跳过组装和比对两步
- **样本名重复怎么办?** 双端样本名与单端文件名撞车会报错;同物种多个体用 `--samples-map` 归组,不要改名硬凑
- **树上有支持度很低的支?** 常见原因:reads 质量差、覆盖度太低、样本太少;优先检查 `input.tsv` 和 fastq 质控
