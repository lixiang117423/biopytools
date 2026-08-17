# vcf-merger - VCF按染色体合并 | VCF Merge by Chromosome

一句话理解：**把按染色体切成多段、散落在一个目录里的 VCF 文件，按染色体自动归组，再用 bcftools 拼回「一整条染色体」的完整 VCF**。解决「逐条染色体跑完变异检测、结果碎成很多小文件，需要合并回完整文件」的问题。

## 功能概述 | Overview { #overview }

- 自动从文件名识别染色体编号并分组（支持 `Chr19`、`chr19`、`19` 等写法）
- 用 `bcftools concat` 把同一染色体的多个 VCF 拼接成一个
- 输出 gzip 压缩格式（`.vcf.gz`），每个染色体一个文件
- 自动为结果创建索引（`.csi`），方便下游按区间快速访问
- 逐染色体报告成功/失败，日志清晰

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-merger -i /path/to/vcf_files -o /path/to/output
```

把输入目录里所有匹配的文件按染色体合并，结果写到输出目录。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| VCF | 记录「基因组上哪里和参考不一样」的标准文件，像一本变异清单 |
| 染色体 | 遗传物质装订成的一本本「册子」，不同物种数量不同 |
| 按染色体切分 | 把一条染色体的变异检测拆成几段并行跑，每段一个文件，跑完要拼回来 |
| bcftools concat | 把多个 VCF「首尾相接」拼成一个的拼接命令 |
| 索引 (.csi) | 给压缩文件做的「目录」，让程序能快速跳到某一段；没有索引只能从头顺序读 |
| 压缩 (.gz) | 文件打包压缩，省磁盘空间，是 VCF 的通用做法 |

## 输入 | Input { #input }

输入是一个**目录**（不是单个文件），目录里放按染色体命名的 VCF 文件，默认只匹配 `*.joint.vcf.gz`。

关键要求：**文件名必须以「染色体编号」开头**，编号后跟下划线或直接到文件名结尾。程序只把第一个下划线之前的部分当作染色体号，其余部分（区间坐标等）忽略。例如：

```text
Chr19_1-20000000.joint.vcf.gz         -> 归到染色体 Chr19
Chr19_20000001-40000000.joint.vcf.gz  -> 归到染色体 Chr19
chr07_1-10000000.joint.vcf.gz         -> 归到染色体 chr07
19_1-20000000.joint.vcf.gz            -> 归到染色体 19
```

注意：同一个染色体编号（如 `Chr19` 和 `chr19`）大小写不同会被当成不同分组，请保持命名一致。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 这两个参数告诉程序「去哪找文件、结果放哪」，每次都必须填。

### 文件匹配 | File pattern

**通俗理解|In plain words:** 决定「哪些文件算要合并的 VCF」。默认 `*.joint.vcf.gz` 只认以 `.joint.vcf.gz` 结尾的文件；如果你的文件是别的命名（比如 `*.vcf.gz`），才需要改。**一般不用动。**

### 并行与索引 | Threads & index

**通俗理解|In plain words:** 线程数只影响 bcftools 合并时的并行度，合并本身不重，默认值即可。索引是给结果文件做的「目录」，**一般保留**（下游工具常用），只有确定用不到时才用 `--no-index` 关掉以省一点磁盘和时间。

### 日志 | Logging

**通俗理解|In plain words:** 控制屏幕输出的啰嗦程度。`-v` 多一次多一层细节，`--quiet` 只报错误。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先找到文件、再按染色体分堆、最后每堆拼一个文件并做索引，全程自动。

```text
输入目录
    │
    ▼
步骤1: 查找匹配 pattern 的 VCF 文件
    │
    ▼
步骤2: 从文件名提取染色体编号，按染色体分组
    │
    ▼
步骤3: 逐染色体用 bcftools concat 合并 -> {chr}.joint.merged.vcf.gz
    │
    ▼
步骤4: 为每个结果创建索引 -> {chr}.joint.merged.vcf.gz.csi
```

## 输出 | Output { #output }

```text
output/
├── Chr19.joint.merged.vcf.gz        # 染色体 Chr19 合并结果
├── Chr19.joint.merged.vcf.gz.csi    # 对应索引
├── chr07.joint.merged.vcf.gz        # 染色体 chr07 合并结果
├── chr07.joint.merged.vcf.gz.csi
└── ...                              # 每条识别出的染色体各一个 .gz + .csi
```

- 合并过程中会在输出目录生成临时的 `{chr}.joint.merged.vcf_filelist.txt`（文件清单），合并成功后自动删除。
- 文件名里的 `joint.merged` 是固定后缀，染色体号来自输入文件名。

## 结果解读 | Interpreting Results { #interpreting-results }

- 每个染色体一个 `.vcf.gz` 文件，就是该染色体所有分段的拼接结果，可直接用于下游分析。
- `.csi` 是索引，不影响内容，删了也能用 `bcftools index` 重建。
- 判断是否成功：看日志末尾的总结，形如「成功合并 N/N 个染色体」；若 N/N 不齐，说明有染色体合并失败（日志里会列出失败原因，通常是文件损坏或格式不一致）。
- 如果某个文件被跳过并出现 WARNING「无法从文件名提取染色体编号」，说明它没被合并——检查它的文件名是否以染色体编号开头。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **`--threads`**：合并是轻量操作，默认 12 通常够用，无需调大。
- **`--pattern`**：只有文件名不以 `.joint.vcf.gz` 结尾时才改，例如统一命名为 `*.vcf.gz` 时改成 `--pattern '*.vcf.gz'`。
- **`--no-index`**：一般不加；保留索引方便下游（如 bcftools view、IGV 等）快速访问。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-dir, -i` | 必填 |  | 输入VCF文件目录｜Input directory containing VCF files |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--pattern` | `*.joint.vcf.gz` |  | VCF文件名模式｜VCF file name pattern |
| `--no-index` | — |  | 不生成索引文件｜Do not generate index files |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose mode |
| `--quiet` | — |  | 静默模式｜Quiet mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 输入VCF文件目录｜Input directory containing VCF files |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `4` | int | 使用的线程数｜Number of threads |
| `--pattern` | `*.joint.vcf.gz` |  | VCF文件名模式｜VCF file name pattern |
| `--no-index` | — | store_true | 不生成索引文件｜Do not generate index files |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | 静默模式｜Quiet mode (only ERROR) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- bcftools（外部命令行工具，通过系统 PATH 直接调用；未安装会报错并提示 `conda install -c bioconda bcftools`）

## 常见问题 | FAQ { #faq }

**Q1：运行提示「未找到匹配的VCF文件」，为什么？**
默认 pattern 是 `*.joint.vcf.gz`，只匹配以 `.joint.vcf.gz` 结尾的文件。若你的文件名不满足，用 `--pattern` 指定；同时确认文件名以染色体编号开头，否则会被分组时跳过。

**Q2：文件没被合并，日志里只有 WARNING？**
程序无法从该文件名提取染色体编号时只警告并跳过。文件名必须像 `Chr19_xxx`、`chr19_xxx` 或 `19_xxx` 这样以编号开头。

**Q3：有断点续传吗？**
没有。每次运行都会重新合并并覆盖输出文件，中途失败后重跑会从头再来。

**Q4：`.csi` 文件能删吗？**
能。它只是索引，删除后可用 `bcftools index your.vcf.gz` 重新生成。
