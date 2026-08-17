# HiTE - 转座子检测与注释 | HiTE Transposable Element Detection and Annotation

一句话理解：**用 singularity 容器跑 HiTE，快速把基因组里的转座子（尤其是全长 LTR）从头识别出来**，输出一份高质量的「可信 TE 库」，可选做全基因组注释和 TE 蛋白结构域预测。

## 功能概述 | Overview { #overview }

- 封装 HiTE 3.3.3（通过 singularity 容器直接挂载调用，无需本地装 HiTE 依赖）
- 全类型 TE 检测：LTR / TIR / Helitron / non-LTR（LINE、SINE），可按 `--te-type` 只跑某类
- 默认植物模式（`--plant`），可关闭用于动物/其它基因组
- 可选全基因组注释（`--annotate` 产出 `.gff` / `.out` / `.tbl`）
- 可选 TE 蛋白结构域预测（`--domain`），用于精细分类
- 模块级断点续传（主结果存在即跳过）+ HiTE 自身 `--recover` 续跑

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools hite -i genome.fa -o hite_results -t 12
```

最小输入：一个基因组 FASTA + 已装好的 singularity 与 HiTE SIF 镜像。结果写入 `<输出目录>/01_hite/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 转座子/TE | 基因组里会「复制并跳到别处」的 DNA 片段，像会自动繁殖的段落 |
| LTR 逆转座子 | 最大的一类转座子，两端有重复的「长末端重复」作为识别标记；HiTE 最擅长找它的全长拷贝 |
| TIR / Helitron / non-LTR | 其它几类转座子；non-LTR 又分 LINE 和 SINE |
| 全长/完整 TE | 结构齐全（首尾都在）的转座子，比残缺碎片更有研究价值 |
| 嵌套 TE | 一个转座子插在另一个转座子「肚子里」，层层套娃；默认会移除，避免重复计数 |
| 结构域(domain) | 蛋白质上的功能模块（如逆转录酶域），是判断 TE 类型与活性的关键证据 |
| singularity / SIF | 一种「把软件连同环境打包成一个文件」的容器技术，SIF 就是那个镜像文件 |

## 输入 | Input { #input }

### 基因组 FASTA

必需，标准 FASTA 格式（可含多序列）：

```text
>chr1
ATGCGTACGTACGTAGCTA...
>chr2
TTGCAAGCTAGCATCGATC...
```

### 可选输入

- `--curated-lib`：可信的人工筛选 TE 库 FASTA，用于预屏蔽已知重复、降低假阳性

## 参数说明 | Parameters { #parameters }

### 容器配置 | Container

**通俗理解|In plain words:** 这两个参数告诉程序「singularity 在哪、HiTE 镜像在哪」。默认值已指向标准安装位置，**装好后一般不用动**；只有 singularity 或 SIF 装在别处时才需改。

本组参数：`--singularity-path`（默认 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`）、`--sif-file`（默认 `~/software/singularity/hite_3.3.3.sif`）。

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 线程数（越大越快）；`-o` 输出目录；`--recover` 是「从上次中断处续跑」（换参数重跑前建议先删旧产物，否则可能复用旧结果）。

本组参数：`-t/--threads`（默认 12）、`-o/--output-dir`（默认 ./hite_output）、`--recover/--no-recover`（默认关）。

### TE 检测范围 | Detection scope

**通俗理解|In plain words:** `--plant` 植物模式（植物基因组开）；`--te-type` 只跑某一类转座子（`all`=全跑，只想找 LTR 就选 `ltr`，更快）；`--min-te-len` 是最短 TE 长度，调小=纳入更短的碎片（更全但更杂）。**默认全类型、植物模式，一般不动。**

本组参数：`--plant/--no-plant`（默认开）、`--te-type`（默认 all，可选 ltr/tir/helitron/non-ltr）、`--min-te-len`（默认 80 bp）。

### 注释与高级 | Annotation & advanced

**通俗理解|In plain words:** `--annotate` 产出全基因组注释(GFF/OUT/TBL)；`--domain` 额外预测 TE 蛋白结构域（更精细但更慢）；`--remove-nested` 移除嵌套 TE（默认开，关掉会保留套娃层）；`--chunk-size` 是基因组分块大小(MB)，超大基因组可调小避免内存爆；`--miu` 是中性突变率用于估算 TE 年龄。**除了 `--annotate`/`--domain` 按需开，其余一般不用动。**

本组参数：`--annotate/--no-annotate`（默认关）、`--domain/--no-domain`（默认关）、`--remove-nested/--keep-nested`（默认开）、`--chunk-size`（默认 400 MB）、`--miu`（默认 1.3e-8）、`--curated-lib`、`--debug/--no-debug`（默认关）。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 程序在主机上组装好参数，用 singularity 把基因组目录和输出目录「挂进」容器，让容器里的 HiTE 直接在主机目录读写，省去来回拷贝。

```text
基因组 FASTA
    |
    v
singularity exec --bind <基因组目录> --bind <输出目录> hite_3.3.3.sif python /HiTE/main.py
    |
    v
HiTE: 分块检测各类型 TE -> 去冗余/移除嵌套 -> 汇总可信 TE 库
    |
    v
写入 <输出目录>/01_hite/（confident_TE.cons.fa 等）
```

## 输出 | Output { #output }

```text
<输出目录>/
├── 00_pipeline_info/
│   ├── software_versions.yml       # 软件版本与关键参数
│   └── hite_summary.json           # 结果汇总(输出文件清单)
├── 01_hite/
│   ├── confident_TE.cons.fa        # 可信 TE 库(核心产物, 断点续传判据)
│   ├── confident_ltr_cut.fa.cons   # 可信 LTR 库
│   ├── confident_tir*.fa           # 可信 TIR
│   ├── confident_helitron*.fa      # 可信 Helitron
│   ├── confident_non_ltr*.fa       # 可信 non-LTR(LINE/SINE)
│   ├── all_TE.fa                   # 全部 TE(含低置信)
│   ├── low_confident_TE.cons.fa    # 低置信 TE
│   ├── HiTE.out / HiTE.gff / HiTE.tbl   # 全基因组注释(--annotate 时)
│   ├── confident_TE.cons.fa.domain      # 结构域预测(--domain 时)
│   └── work/                       # 中间文件
└── 99_logs/
    └── hite_<时间戳>.log           # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 可信 TE 库（`confident_TE.cons.fa`）

**通俗理解|In plain words:** HiTE 筛掉假阳性和碎片后的「可信转座子名单」，是这份分析最核心的产出。

- 序列条数与平均长度反映 TE 多样性；与 `all_TE.fa`（未过滤）对比，可信库应明显更「干净」
- 各类型分库文件（ltr/tir/helitron/non_ltr）方便按类统计

### 2. 全基因组注释（`--annotate` 时）

- `HiTE.gff`：每个 TE 的坐标+类型；`HiTE.out`：RepeatMasker 风格注释；`HiTE.tbl`：统计表
- 用 GFF 统计 TE 总占比（多数植物 30%-80%）；占比异常低提示可能漏检，异常高提示可能假阳性

### 3. 结构域文件（`--domain` 时）

给每条 TE 标注了含哪些蛋白结构域，可用于判断 TE 是否「活着」（含完整酶结构域通常仍有活性）。

## 参数选择建议 | Parameter Guidance { #guidance }

- **`--annotate`**：需要 TE 坐标做下游分析（如统计、屏蔽、进化）时开启；只想要 TE 库可不开，更快
- **`--domain`**：研究 TE 活性/进化时开启；纯屏蔽用途可不开
- **`--te-type`**：只关心 LTR 就选 `ltr`，显著省时；大基因组可先按类型分跑
- **`--chunk-size`**：超大基因组（>2 Gb）或内存紧张时调小（如 200）避免爆内存；小基因组保持默认
- **`--recover`**：HiTE 中途被杀，用 `--recover 1` 从断点续跑；但换参数重跑应先删 `01_hite/` 旧产物
- **`--plant`**：动物/真菌等非植物基因组加 `--no-plant`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./hite_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--singularity-path` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity可执行文件｜Singularity executable |
| `--sif-file` | `~/software/singularity/hite_3.3.3.sif` |  | HiTE SIF镜像｜HiTE SIF image |
| `--plant/--no-plant` | `True` |  | 植物基因组｜Plant genome |
| `--annotate/--no-annotate` | `False` |  | 注释基因组(产出.gff/.out/.tbl)｜Annotate genome |
| `--recover/--no-recover` | `False` |  | HiTE断点续跑｜HiTE recovery mode |
| `--domain/--no-domain` | `False` |  | 预测TE蛋白结构域｜Predict TE domains |
| `--te-type` | `all` | ltr/tir/helitron/non-ltr/all | TE类型｜TE type |
| `--chunk-size` | `400` | int | 基因组分块MB｜Genome chunk size MB |
| `--miu` | `1.3e-08` | float | 中性突变率｜Neutral mutation rate |
| `--min-te-len` | `80` | int | 最小TE长度bp｜Min TE length bp |
| `--remove-nested/--keep-nested` | `True` |  | 移除嵌套TE｜Remove nested TE |
| `--curated-lib` | — |  | 可信curated TE库(预屏蔽)｜Trusted curated TE library |
| `--debug/--no-debug` | `False` |  | HiTE debug模式(保留临时文件)｜Debug mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | `./hite_output` |  | 输出目录｜Output directory (default: ./hite_output) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--singularity-path` | — |  | Singularity可执行文件路径｜Singularity executable path |
| `--sif-file` | — |  | HiTE SIF镜像路径｜HiTE SIF image path |
| `--plant` | `1` | 0/1 | 是否植物基因组(1/0)｜Plant genome (1/0, default: 1) |
| `--annotate` | `0` | 0/1 | 是否注释基因组(1/0)｜Annotate genome (1/0, default: 0) |
| `--recover` | `0` | 0/1 | HiTE断点续跑(1/0)｜HiTE recovery mode (1/0, default: 0) |
| `--domain` | `0` | 0/1 | 预测TE蛋白结构域(1/0)｜Predict TE domains (1/0, default: 0) |
| `--te-type` | `all` | ltr/tir/helitron/non-ltr/all | TE类型｜TE type (default: all) |
| `--chunk-size` | `400` | int | 基因组分块MB｜Genome chunk size MB (default: 400) |
| `--miu` | `1.3e-08` | float | 中性突变率｜Neutral mutation rate (default: 1.3e-8) |
| `--min-te-len` | `80` | int | 最小TE长度bp｜Min TE length bp (default: 80) |
| `--remove-nested` | `1` | 0/1 | 移除嵌套TE(1/0)｜Remove nested TE (1/0, default: 1) |
| `--curated-lib` | — |  | 可信curated TE库｜Trusted curated TE library |
| `--debug` | `0` | 0/1 | HiTE debug模式(保留临时文件)｜Debug mode (1/0) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- singularity（默认 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`，可用环境变量 `SINGULARITY_PATH` 覆盖）
- HiTE SIF 镜像 `hite_3.3.3.sif`（默认 `~/software/singularity/hite_3.3.3.sif`，可用环境变量 `HITE_SIF` 覆盖）
- HiTE 版本 3.3.3（镜像内置，无需单独安装）

## 常见问题 | FAQ { #faq }

**Q1：为什么结果直接写到 `01_hite/`，而不是容器里？**
采用 singularity 直接挂载模式（`--bind 主机目录:主机目录`），容器内 HiTE 直接读写主机目录，省去来回拷贝；输出目录无需二次同步。

**Q2：重跑为什么很快结束？**
模块级断点续传：`01_hite/confident_TE.cons.fa` 已存在就跳过整个检测。换参数重跑前先删该文件（或整个 `01_hite/`）。

**Q3：`--recover` 和断点续传什么关系？**
`--recover` 是 HiTE 自身的「从上次中断续跑」，针对一次没跑完的进程；模块级断点续传是「结果已生成就跳过」。两者互补：跑挂了用 `--recover`，跑完了重跑会直接跳过。

**Q4：报 singularity 或 SIF 不存在？**
检查默认路径是否与实际安装位置一致，或用 `--singularity-path` / `--sif-file` 指定，或设置环境变量 `SINGULARITY_PATH` / `HITE_SIF`。

**Q5：`--te-type ltr` 为什么仍可能产出其它类型？**
`--te-type` 只限制 HiTE 主要检测哪类；汇总时仍会保留少量边缘类型文件，属正常现象。
