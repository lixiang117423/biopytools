# vcf-renamer - VCF样品名称重命名 | VCF Sample Name Renamer

一句话理解：**把 VCF 里又长又乱的样本名统一改成简洁的 S1、S2、S3……**，避免下游软件因样本名过长或含特殊字符而截断/报错，同时保留一张「旧名对应新名」的对照表。

## 功能概述 | Overview { #overview }

- 用 `bcftools query -l` 自动提取 VCF 里的所有样本名
- 生成 `S1`、`S2`、`S3` 这样的序号新名（前缀可自定义）
- 用 `bcftools reheader` 完成重命名，不改动任何基因型数据
- 生成新旧样本名映射文件（可保留或删除）
- 自动为输出 VCF 创建索引（`.tbi`）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-renamer -i input.vcf.gz -o output.vcf.gz
```

把 `input.vcf.gz` 里的样本名按出现顺序改成 S1、S2、S3……，写到 `output.vcf.gz`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 样本名 | VCF 里每个被测个体在表头里的名字，像表格的列名 |
| reheader | 只改 VCF 的「表头」（样本名列），不动里面的数据，像只换列名不换内容 |
| 映射文件 | 一张「旧名 → 新名」的两列对照表，方便之后对号入座 |
| 索引 (.tbi) | 给结果文件做的「目录」，让程序能快速定位到某段 |
| 截断 | 某些软件对超长样本名只取前若干字符，导致两个样本被当成同一个 |

## 输入 | Input { #input }

输入一个 VCF 文件，支持 `.vcf` 或 `.vcf.gz`。样本名不限格式（无论多长、是否含特殊字符都能被读出来并重命名）。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 指定输入 VCF 和输出 VCF 的路径，每次必填。

### 命名规则 | Naming

**通俗理解|In plain words:** 决定新名字长什么样。默认前缀 `S`，结果是 S1、S2、S3……；改成 `--prefix IND` 就是 IND1、IND2……。**一般用默认 S 即可。**

### 映射文件 | Mapping file

**通俗理解|In plain words:** 控制「旧名对新名」的对照表怎么处理。默认保留一张映射文件方便追溯；如果不需要就加 `--no-mapping`（跑完自动删除）。想自己指定对照表保存位置用 `-m`。

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF
    │
    ▼
步骤1: 检查 bcftools 是否可用
    │
    ▼
步骤2: bcftools query -l 提取样本名
    │
    ▼
步骤3: 生成 S1/S2/S3 新旧名映射表
    │
    ▼
步骤4: bcftools reheader -s 映射表 完成重命名
    │
    ▼
步骤5: bcftools index -t 创建索引 -> output.vcf.gz.tbi
```

## 输出 | Output { #output }

以 `biopytools vcf-renamer -i input.vcf.gz -o /path/output.vcf.gz` 为例：

```text
/path/output.vcf.gz              # 重命名后的 VCF
/path/output.vcf.gz.tbi          # 索引(bcftools index -t 生成)
/path/vcf_renamer.log            # 运行日志(写在输出 VCF 同目录)
./output_name_mapping.txt        # 新旧样本名对照表(默认写在「当前工作目录」!)
```

- 映射文件默认名由输出文件名生成：`output.vcf.gz` → `output_name_mapping.txt`，两列制表符分隔，第一列旧名、第二列新名。
- 加了 `--no-mapping` 时，映射文件跑完会被删除。

## 结果解读 | Interpreting Results { #interpreting-results }

- 输出 VCF 的基因型数据与输入完全一致，只是样本名变成 S1、S2、S3……（按 VCF 里的出现顺序编号）。
- 打开映射文件核对：第一列是原始样本名，第二列是新的序号名。后续分析要用「哪个样本对应哪个新名」时查这里。
- 判断成功：日志末尾显示「重命名完成」及样本数量；`.tbi` 索引同时生成。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **`--prefix`**：默认 `S` 满足绝大多数需求；需可读性更强的名字（如样本名里有分组信息）再改前缀。
- **`-m/--mapping`**：建议显式指定映射文件路径（例如 `-m /path/output_name_mapping.txt`），避免它默认落到「当前工作目录」而找不到。
- **`--no-mapping`**：确定不需要追溯新旧名对应关系时再加。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | 必填 | Path | 输出VCF文件路径｜Output VCF file path |
| `--prefix, -p` | `S` |  | 新样本名前缀｜New sample name prefix |
| `--mapping, -m` | — | Path | 样本映射文件路径｜Mapping file path |
| `--no-mapping` | — |  | 不保留映射文件｜Do not keep mapping file |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | 必填 |  | 输出VCF文件路径｜Output VCF file path |
| `-p, --prefix` | `S` |  | 新样品名前缀｜New sample name prefix |
| `-m, --mapping` | — |  | 映射文件路径｜Mapping file path |
| `--no-mapping` | — | store_true | 不保留映射文件｜Do not keep mapping file |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- bcftools（外部命令行工具，通过 `which bcftools` 检查是否在 PATH 中）

## 常见问题 | FAQ { #faq }

**Q1：为什么重命名后样本名会「截断」？**
这通常发生在没做重命名、直接拿原始长样本名喂给下游软件时。本工具就是为了避免这种截断：统一改成短序号名 S1/S2/S3。

**Q2：映射文件在哪？**
默认写在**当前工作目录**（不是输出 VCF 的目录），文件名是 `输出名_name_mapping.txt`。想放别处用 `-m` 显式指定完整路径。

**Q3：有断点续传吗？**
没有。每次运行都会重新生成映射和输出，覆盖旧结果。

**Q4：重命名会改动基因型数据吗？**
不会。`bcftools reheader` 只改表头里的样本名，所有变异和基因型内容原样保留。
