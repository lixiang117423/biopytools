# Phylip 转 FASTA | phy2fa

把 Phylip 格式的多序列比对矩阵转成 FASTA 格式,一个命令的事——Phylip 每行一个样本"名称+序列",FASTA 是"&gt;名称"换行接序列,转换规则简单直接。

## 功能概述 | Overview { #overview }

- Phylip(.phy/.phylip,可 .gz)→ FASTA,自动处理首行头(`n_tax n_char`)
- 每行按空白拆分:第一个字段是样本名,其余拼接为序列;同名出现自动追加(兼容 interleaved 分块)
- 输出 FASTA 换行宽度可调(默认 60,0=不换行)
- 转换后按头行校验样本数与序列长度,不符报错(数据完整性保护)
- 断点续传:输出已存在则跳过重算
- 输出含 `00_pipeline_info/software_versions.yml` 与 `99_logs/`

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools phy2fa -i aln.phy -o fa_out/
```

得到 `fa_out/aln.fa`,可直接用于需要 FASTA 的下游工具。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| Phylip 格式 | 第一行写"样本数 序列长度",之后每行一个样本:名称 + 序列字符(可跨多列) |
| FASTA 格式 | 每样本两行:以 `>` 开头的名称行 + 序列行 |
| interleaved | Phylip 的一种排版:每行只放序列的一部分,多行拼成完整序列(同名重复出现) |

## 输入 | Input { #input }

- `.phy` / `.phylip`(可带 `.gz`)
- 首行:`样本数 序列长度`(如 `365 147433`)
- 数据行:名称 + 空白分隔的序列字符(名称在前,无严格列宽要求)

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** 只有输入输出和换行宽度三个参数。换行宽度只影响 FASTA 观感,不影响内容;下游工具对行宽有要求时才需要调。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | Phylip 文件(.phy/.phylip/.gz)｜Phylip file |
| `-o, --output-dir` | `./phy2fa_output` |  | 输出目录(默认./phy2fa_output)｜Output directory |
| `--line-width` | `60` | int | FASTA 换行宽度,0=不换行(默认60)｜FASTA line wrap width, 0=no wrap (default 60) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | Phylip 文件(.phy/.phylip/.gz)｜Phylip file |
| `-o, --output-dir` | `./phy2fa_output` |  | 输出目录(默认./phy2fa_output)｜Output directory |
| `--line-width` | `60` | int | FASTA 换行宽度,0=不换行(默认60)｜FASTA line wrap width, 0=no wrap (default 60) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 输出 | Output { #output }

```
output/
├── aln.fa                     # FASTA(>名称 换行 序列)
├── 00_pipeline_info/
│   └── software_versions.yml  # 输入/输出/样本数/序列长度/耗时
└── 99_logs/
    └── phy2fa.log             # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **样本数**:`grep -c '^>' aln.fa` 应等于 Phylip 头行的第一个数
- **序列长度**:FASTA 每条序列长度应一致,等于头行第二个数(模块已校验,不符会报错)
- **名称**:与 Phylip 中的名称逐字一致(含纯数字/字母等)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 建议<br>Recommendation |
|------|------|
| 常规转换 | 默认 `--line-width 60` 即可 |
| 下游要求单行序列 | `--line-width 0` |
| 超大文件 | 直接跑,解析是流式逐行的 |

**通俗理解|In plain words:** 不需要调参,默认就能用。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| Python 3.10+ | 无第三方依赖 |

## 常见问题 | FAQ { #faq }

- **报"样本数不符"?** 文件数据行数与头行第一个数不一致,检查文件是否被截断
- **报"序列长度不符"?** 某样本序列长度与头行第二个数不同,常见于手工编辑或样本缺失片段
- **`>365` 出现在输出里?** 说明头行没被跳过——文件首行可能不是 `n_tax n_char`,确认文件格式
- **支持 .gz 吗?** 支持,自动解压
