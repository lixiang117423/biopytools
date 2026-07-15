
## [1.13.0] - 2026-07-15

### Added
- `ps_gene_anno`：gap 验证报告新增 StringTie FPKM/TPM 定量——对 `gap_filled.gff3` 跑 `stringtie -e -G` 定量每个 gap 转录本，报告新增 `fpkm`/`tpm` 两列；transcript id 按 `{prefix}_gap_{N}.t1` 与 `build_gene_models` 建模型命名对齐；新增 `stringtie_bin` 配置（`get_tool_path`+`~`展开，默认 `~/.local/bin/stringtie`，env var `STRINGTIE_PATH`）；StringTie 命令走 `cmd_runner.run_command`（自动 conda 包装 + 记录完整命令）；`_parse_stringtie_fpkm` 兼容 GTF(`key "v"`) 与 GFF3(`key=v`) 两种属性格式

### Changed
- `ps_gene_anno`：gap 报告列从 `rnaseq_mean_depth/te_overlap_pct/te_family` 调整为 `rnaseq_mean_depth/fpkm/tpm/te_overlap_pct/te_family`；FPKM/TPM 与 raw depth 在无 RNA-seq BAM 或无 gap_filled 时置 0

## [1.12.0] - 2026-07-15

### Added
- **新模块**：`mga` — MGA 共识基因组组装（HiFi，consensusLJA 封装）；二进制不在 conda env 但运行时依赖 env 内 minimap2/samtools/python，故显式 `conda run -n <env> --no-capture-output` 包装；含断点续传（`5_polishing/assembly.fasta` 存在即跳过）、read name 空格预检、dry-run、`00_pipeline_info/software_versions.yml`
- `ps_gene_anno`：新增 gap 验证报告（`{prefix}.gap_report.tsv`）——每 gap 汇总蛋白证据 + RNA-seq mean depth（复用 braker 的 `compute_region_mean_depth`）+ TE 重叠% 与 family，report-only 不改 gap 结果
- `ps_gene_anno`/`braker4ps`：新增 `--exclude-te-gap`（质控排除 TE 区 gap；默认不排——疫霉效应子常落在 TE 区）
- `braker4ps`：阶段2 自动探测 braker 的 RepeatMasker `.out`（`01_repeat_masking/{genome}.out`）与 RNA-seq BAM（`03_short_reads/rnaseq.sorted.bam`）供 gap 报告使用（用户 `--repeat-out` 优先）

### Changed
- `ps_gene_anno`：`parse_repeat_out` 返回值增加 TE family 列（`(start,end)` → `(start,end,family)`，取 cols[10]）；header 跳过改为 `lstrip()` 后判 `SW`/`score`（修复原 `startswith('SW')` 因前导空格漏跳 header、靠 except 兜底的问题）；`qc_filter` 的 TE 排除默认关闭并适配 3 元组返回；`qc_filter` 类型注解同步为 `Tuple[int,int,str]`
- CLI 注册新命令 `mga`

## [1.11.1] - 2026-07-15

### Fixed
- `ps_gene_anno`：`detect_gaps` 命中覆盖判定 `>=` 改为 `>`，修复 `overlap_cutoff=0`（1.11.0 新默认值）下 `overlap_ratio 0.0 >= 0` 恒为真、所有 miniprot hit 被判"已覆盖"、漏检检测永远返回空的 bug（零重叠现在正确判为漏检）

### Added
- `fastp`：目录中无 `_1/_2` 配对文件时自动检测单末端文件（如 PacBio HiFi 的 `*.fq.gz`），命中即切换单末端模式（仅在用户未显式指定 `--single-end`/`--read1-suffix` 时生效）

## [1.11.0] - 2026-07-14

### Added
- **新模块**：`seq_extract` — seqkit 封装的序列提取，自动识别单 ID / ID 文件（一列）/ BED 文件（≥2 列），输出文件名自动推导（`{query}.{subject}.fa`）
- **新模块**：`raxml_ng` — RAxML-NG 最大似然系统发育树（`all`/`search`/`support` 三模式，原生断点续传不传 `--redo` 即续跑，支持 `--bs-trees`/`--bs-metric`/`--outgroup`/`--seed`）
- `transcript_assembly`：新增 step7 TransDecoder CDS 预测（`--predict-cds`，需 `-g`；输出基因组坐标 gene/mRNA/CDS GFF3 + pep + cds；复用既有 `cmd_runner` 记录完整命令）
- `bam2fastq`：`-o` 支持文件输出（自动识别 `.fq`/`.fastq`/`.fq.gz`/`.fastq.gz` 后缀为文件，否则为目录；多 read group 输出自动合并；目录输入遇单文件输出自动回退目录模式）

### Changed
- `ps_gene_anno`：全 prot 场景适配——同位置多 query 命中去重（`dedupe_hits`，CDS 重叠≥50% 合并保留 identity/coverage 最优）；错误合并基因改为按 query 分组检测（同一 query 的 ≥N 完整独立拷贝才算合并，避免不同 query 混合 pairwise 误判重叠）；`gap_min_cds_len` 100→300（过滤短蛋白片段）、`overlap_cutoff` 50→0（零重叠才算漏检）
- `braker4ps`：阶段2 输出目录 `gap_filling` → `05_gap_filling`（§12.2.2 序号前缀）
- `raxml_ng`：`get_conda_env` 对绝对路径（静态二进制）跳过全 env 搜索，修复 `os.path.join` 塌缩导致返回错误环境名的隐患（与 `indel_marker`/`braker`/`phobius` basename 防护同源）
- CLI 注册新命令 `seq-extract`、`raxml-ng`
- `biopytools.yml` 补 `argcomplete`/`colorama`/`pipx`/`platformdirs`/`tomli`/`userpath` 依赖

## [1.10.0] - 2026-07-13

### Added
- `transcript_assembly`：BAM 直入模式（`-b` 可多次，跳过索引/比对/排序步骤 1-3，直接组装）
- `transcript_assembly`：长读支持（`--read-type auto` 采样读长中位数、阈值 500bp 自动判短/长，StringTie 加 `-L`）
- `transcript_assembly`：guided 组装（`--guide-gff` 传参考注释）；主输出改 GFF3（gffread 转，无需基因组），cDNA 改为可选（`--transcripts`，需 `-g`）
- `transcript_assembly`：单样本自动跳过 merge；新增 `00_pipeline_info/software_versions.yml`

### Fixed
- `transcript_assembly`：StringTie `--merge -G` 原误传 `genome.fa`，改为传注释 GFF（仅提供时）

### Changed
- `transcript_assembly`：所有 `build_conda_command` 改传完整工具路径（原为裸命令名，违反 §13.6.1）；工具路径走 `get_tool_path`/`get_samtools_path`，`__post_init__` 统一 `expand_path`（§11.3.1）；BAM 读长检测管道用单 conda run + 系统 head/awk（§13.2.1）

## [1.9.0] - 2026-07-13

### Added
- **新模块**：`ps_gene_anno` — BRAKER 注释后用 miniprot 补回效应子多拷贝/漏注位点（config/evidence/gap_analysis/model_build/merge 分层）
- **新模块**：`braker4ps` — 端到端编排器，仅 import 调用 braker + ps_gene_anno（不改源码），统一日志（阶段1 注释 + 阶段2 查漏补缺）
- conda 环境定义：新增 `deeptmhmm_v.1.0` / `eggnog-mapper_v.2.1.15` / `phobius_v.1.0.1` / `starship_v.1.26.0` / `zsh`

### Changed
- CLI 注册新命令 `ps-gene-anno`、`braker4ps`
- `biopytools.yml` / `primer3_v.2.6.1.yml` 补 `primer3-py` 依赖

## [1.8.0] - 2026-07-13

### Added
- `oomycete_anno`：新增步骤9 效应子位点救援（Phase3，可选、失败不阻断主注释）——已知效应子全长 miniprot 比对当基因模型，替换/补回 Augustus 在效应子簇位点的错注（嵌合/截断）/漏注；全长判定靠结构完整性（Target 起=1 + stop_codon），不依赖高 identity
- 新增 `--effectors` / `--skip-rescue` / `--rescue-min-identity` / `--rescue-conflict-overlap` 参数与 `09_effector_rescue` 输出目录

### Fixed
- `braker`：repeat_refine 的 `prot_seq` 转绝对路径再喂 mmseqs（命令在 output_dir 作 cwd 执行，相对路径会失效）

## [1.7.0] - 2026-07-11

### Added
- **新模块**：`gene_table` — 基因信息+序列合并表（基因 DNA + CDS + 蛋白）
- **新模块**：`oomycete_anno` — 疫霉菌基因组注释（T2T Augustus 流程）
- `braker`：repeat_refine 新增 mmseqs 蛋白同源判据（consensus 翻译与近缘蛋白同源 → 假重复核心证据）；RxLR-EER 改紧邻结构扫描（下游 30aa 内含 EER），避免长 ORF 中随机 R-x-L-R 造成过度剔除；miniprot identity 增加 matches/alnlen 兜底（无 `id:f:` tag 时）

### Changed
- `braker`：`--skip_rescue` 改 `BooleanOptionalAction` 默认关闭（filter 库级已处理假重复）；safe_dir 残留时清理全新运行，避免 `--useexisting` 用旧状态
- CLI 注册新命令 `gene-table`、`oomycete-anno`

### Fixed
- `indel_marker`：恢复 Method2 兜底分支 basename 防护，与 braker/phobius 一致（修复绝对路径被 `os.path.join` 劫持、返回错误环境名的隐患）
- `bam-cov`：help 示例精简为单行

## [1.6.0] - 2026-07-09

### Added
- **新模块**：`indel-marker` — 抗病/感病INDEL共显性标记开发（VCF提取INDEL → 群体共显性判定 → 覆盖度/deletion骤降质控 → 侧翼提取 → primer3引物设计 → 候选主表/报告）
- **新模块**：`phobius` — Phobius跨膜拓扑+信号肽预测（与deeptmhmm互补的经典算法，输出 summary.tsv / 3line / gff3）
- `braker`：新增重复序列误伤修复——repeat库智能过滤（方案1：六框翻译 + hmmscan 扫 Pfam-A 的 TE域排除 + RxLR motif）+ 证据驱动还原（方案2：miniprot 蛋白 / samtools depth 转录组 unmask），解决多拷贝基因家族（疫霉RXLR/CRN效应子、抗病基因）被 RepeatModeler 误判为重复序列、导致 BRAKER 注释丢失的问题

### Changed
- `braker`：config 新增 repeat_refine 工具路径（hmmscan/miniprot/samtools）与过滤/还原阈值参数，samtools 走 `get_samtools_path()` 兜底，`__post_init__` 统一 `expand_path`
- CLI 注册新命令 `indel-marker`、`phobius`

## [1.5.0] - 2026-07-08

### Added
- **新模块**：`eggnog-mapper` — eggNOG 功能注释（GO/KEGG/COG/CAZy/Pfam），支持 mmseqs/diamond/hmmer 等搜索模式、断点续传、software_versions.yml，并附带注释结果重排版（中文 TSV + 双语 Excel）

### Fixed
- `ldblockshow`：GWAS P 值文件预过滤 NA/0/负值/Inf/非数值行，避免 ShowLDSVG 取 -log10(P) 时 Perl `log(0)` 崩溃；全文件无有效 P 值时自动禁用 GWAS 轨道；多 region 批处理按源文件去重，仅过滤一次

### Changed
- `ldblockshow`：BED 批处理改用标准 `LDBlockShowLogger` 替代 `print`，output_dir 创建失败时回退 print，超算 .out/.err 日志更规范

## [1.4.0] - 2026-07-02

### Added
- **新模块**：`deeptmhmm` — DeepTMHMM 跨膜螺旋/信号肽预测，输出整理后的 summary.tsv / topologies.3line / tmr.gff3
- **新模块**：`gene_density` — 基因密度计算
- **新模块**：`aliner` — a-liner 共线性可视化 pipeline（FASTA → minimap2 → 图）
- `assembly_qc`：LAI 流程新增 EDTA 序列名自动改名（兼容 EDTA/RepeatMasker 的 13 字符序列名限制，生成唯一短名副本并保留旧名→新名映射，已存在则复用）
- `longestmrna`：新增 GFF 与基因组序列名对齐，跳过基因组中缺失的序列
- `ldblockshow`：新增 `--no-snp-filter` 参数，默认禁用 SNP 过滤
- `ragtag`：合并 scaffolded 与 unscaffolded 序列为单一 FASTA 输出
- `fastp`：单文件输入自动检测单末端模式，无需手动指定 `--single-end`
- `smudgescope`：`get_conda_env` 支持 preferred 环境锁定版本（超算固定 genomescope_v.2.0.1）

### Changed
- **输出目录命名统一为下划线 `NN_xxx`**：`common/paths` 新增 `resolve_legacy_path` / `resolve_legacy_path_chain`，优先下划线规范名、找不到时回退点号老名以支持断点续传；转换 11 个点号派模块（assembly_qc / braker / hicanu / hifi_hic_workflow / dual_rnaseq / fastq2vcf_gtx / rnaseq / hifi_hic / haphic / fastq2vcf_parabricks / wgsim），并同步更新对应 docs
- `smudgescope`：用 conda 环境自动检测替代手动定位 smudgeplot 目录

### Fixed
- `tmhmm`：路径导入增加 fallback，兼容 common 模块不可用场景
- `rmvp`：输入 VCF/表型路径解析为绝对路径，修复 `-o` 子目录时 PLINK 找不到文件
- `hic_heatmap`：补 samtools faidx 命令日志（§7.5）

### Housekeeping
- 移除代码中的 emoji（§5.2），保留数学/科学符号
- 修复多模块硬编码绝对路径、conda 缺 `--no-capture-output` 等规范违规


## [1.3.0] - 2026-06-25

### Added
- `admixture`：新增 LD 剪枝步骤（`--indep-pairwise`）
- `rmvp`：多表型时自动按显著性合并 GWAS 结果；kinship/PCA 改用 LD 去连锁 SNP 计算（GWAS 仍用全部 SNP）；GWAS 输出改流式转发并修复断点续传检测

### Changed
- `bwa_gatk`、`fastq2vcf_gtx`：输出目录改用 §12 规范命名；`fastq2vcf_gtx` 生成软件版本元数据

### Fixed
- `rmvp`：适配 rMVP 1.4.6 真实输出格式，改进断点续传可靠性
- `admixture`：`build_conda_command` 补 `--no-capture-output`（§13.2.0）、修复 PlinkProcessor 命令构造
- `bwa_gatk`：修复 CommandRunner 成功返回值

### Housekeeping
- 同步代码规范 CLAUDE.md 至 v2.14，整理 .gitignore
- 更新 conda 环境配置文件，清理 prefix 硬编码路径


## [1.2.0] - 2026-06-24

### Added
- **新模块**：`qiime2` — QIIME2 微生物组多样性分析
- **新模块**：`subgenome_assign` — 基于亲本比对的亚基因组归属
- **新模块**：`faprotaxtax` — FAPROTAX 微生物群落功能注释

### Fixed
- `kaks`：改用输出内容检测替代退出码检测 KaKs_Calculator 安装
- `cphasing`：移除 conda run 包装改用 `activate_cphasing` 激活；加 stderr 异常扫描防御上游吞异常 bug；修复透传参数 bug 并默认开启亚基因组聚类

### Housekeeping
- README 添加 AI 辅助开发声明


## [1.1.0] - 2026-06-17

### Added
- **新模块**：`minibwa` 短读长比对（标准 / Hi-C / BS-seq / 长读四种模式），封装 Minibwa 工具
- 补提 4 个被注册命令依赖的底层模块目录（之前未入库）：
  - `tgsgapcloser/`：`gap-fill` 命令的底层
  - `hicpro_qc/` 和 `pairtools_qc/`：`hic-qc` 命令的底层
  - `fasta_id_renamer/`：`rename-genome-id` 命令的底层

### Changed
- `picrust2`：将原本只注释通路表的功能扩展到通路/EC/KO 三类功能丰度表，改用 PICRUSt2 官方 `add_descriptions.py` 添加描述列

### Housekeeping
- 移除仓库中遗留的 backup 文件（`*.bak`、`*.backup`、`_kmer_count_backup/` 等）


## [1.0.0] - 2026-06-15

### Added
- 大规模模块扩充：新增 150+ 个生物信息学模块，覆盖完整基因组学研究流程
- 模块按功能分为 19 个类别：数据下载与质控、基因组组装、组装评估与QC、
  Hi-C与挂载、基因组后处理、比对与BAM处理、变异检测与VCF、泛基因组、
  注释与功能预测、转座子与重复序列、RNA-seq与转录组、共线性与比较基因组、
  系统发育、群体遗传、GWAS与BSA、微生物组与k-mer、甲基化、效应子与抗病、
  其他工具
- 补全 59 个模块的中文文档（docs/），所有 175 个 active CLI 命令均有对应文档
- README 重构：按类别组织所有模块的文档索引
- biopytools/common/ 和 biopytools/core/ 通用工具模块
- 完善的 .gitignore，排除常见的临时/输出/备份文件

### Changed
- biopytools/cli/main.py: 注册所有新增模块的 CLI 命令
- biopytools/__init__.py: 更新包入口
- Python 最低版本要求提升至 3.10

### Notes
- 1.0.0 标志着 BioPyTools 进入稳定发行阶段，模块体系趋于完整
- 后续版本将以模块优化、bug 修复、新功能增量为主


## [0.27.2] - 2026-06-05

### Changed
- fix(fastp): sync CLI args with module interface, fix simulated data detection and wgsim quality compatibility
- Updated files: biopytools/fastp


## [0.27.1] - 2026-06-04

### Changed
- fix(rnaseq): 大基因组拆分HISAT2+samtools管道并取消默认超时限制

   大基因组（>1GB）的SAM输出巨大，管道缓冲区溢出导致所有样本
   比对超时。改为运行时根据基因组文件大小自动选择策略：大于阈值
   走两步拆分（hisat2 -S输出SAM，再samtools sort转BAM）；
   小基因组保持原管道方式。同时将sample_timeout默认值从21600秒
   改为None，不再强制超时限制。
- Updated files: biopytools/rnaseq


## [0.27.0] - 2026-06-04

### Changed
- add braker module
- Updated files: biopytools/braker/,docs/braker.md,biopytools/cli/commands/braker.py,biopytools/cli/main.py


## [0.26.3] - 2026-06-04

### Changed
- feat(blast): auto-detect blast-type from input file sequence types

  Scan all input files individually to determine nucleotide/protein type
  by majority vote, then infer blastn/blastp/blastx/tblastn accordingly.
  Falls back to blastn on failure. Warns on mixed sequence types.
  Also update CLI wrapper to pass None defaults for auto-detection.
- Updated files: biopytools/blast


## [0.26.2] - 2026-06-03

### Changed
- feat(blast): auto-detect blast-type from input file sequence types

  When --blast-type is not specified, read first few sequences from -i
  and -r files to determine nucleotide/protein, then infer blastn,
  blastp, blastx, or tblastn accordingly. Falls back to blastn on failure.
- Updated files: biopytools/blast


## [0.26.1] - 2026-06-03

### Changed
- feat(annovar): add reference and mutant protein sequences to exonic variant results
- Updated files: biopytools/annovar


## [0.26.0] - 2026-05-26

### Changed
- add deeploc module
- Updated files: biopytools/deeploc,docs/deeploc.md,biopytools/cli/commands/deeploc.py,biopytools/cli/main.py,README.md


## [0.25.3] - 2026-05-14

### Changed
- update blast module
- Updated files: biopytools/blast


## [0.25.2] - 2026-05-14

### Changed
- fix(annovar): 修复gff3ToGenePred因错误上限导致大部分染色体基因模型丢失|Fix gff3ToGenePred silent data loss from default error limit
- Updated files: biopytools/annovar,biopytools/cli/commands/annovar.py


## [0.25.1] - 2026-05-12

### Changed
- fix(annovar): 修复gff3ToGenePred因错误上限导致大部分染色体基因模型丢失|Fix gff3ToGenePred silent data loss from default error limit
- Updated files: biopytools/annovar,biopytools/cli/commands/annovar.py


## [0.25.0] - 2026-05-06

### Changed
- add gff-renamer module
- Updated files: biopytools/gff_renamer,docs/gff_renamer.md,biopytools/cli/commands/gff_renamer.py,biopytools/cli/main.py,README.md


## [0.24.6] - 2026-04-29

### Changed
- update README
- Updated files: README.md


## [0.24.5] - 2026-04-29

### Changed
- update cim module
- Updated files: biopytools/cim


## [0.24.4] - 2026-04-28

### Changed
- update annovar module
- Updated files: biopytools/annovar


## [0.24.3] - 2026-04-27

### Changed
- update cim module
- Updated files: biopytools/cim,biopytools/cli/commands/cim.py


## [0.24.2] - 2026-04-27

### Changed
- add fixGenoError to CIM pipeline for genotype correction
- Updated files: biopytools/cim


## [0.24.1] - 2026-04-27

### Changed
- update conda env yaml
- Updated files: conda_env


## [0.24.0] - 2026-04-27

### Changed
- add cim module
- Updated files: biopytools/cim,docs/cim.md,biopytools/cli/commands/cim.py,biopytools/cli/main.py,conda


## [0.23.3] - 2026-04-13

### Changed
- update fastp module: replace repair.sh with seqkit pair
- Updated files: biopytools/fastp,biopytools/cli/commands/fastp.py


## [0.23.2] - 2026-04-13

### Changed
- update fastp module
- Updated files: biopytools/fastp,biopytools/cli/commands/fastp.py


## [0.23.1] - 2026-04-11

### Changed
- Updated files: README.md,update README


## [0.23.0] - 2026-04-11

### Changed
- add merge_deepbsa module
- Updated files: biopytools/merge_deepbsa,docs/merge_deepbsa.md,biopytools/cli/commands/merge_deepbsa.py,biopytools/cli/main.py


## [0.22.2] - 2026-04-02

### Changed
- fix ss file bug for dual-rnaseq
- Updated files: biopytools/dual_rnaseq


## [0.22.1] - 2026-04-02

### Changed
- fix ss file bug for dual-rnaseq
- Updated files: biopytools/dual_rnaseq


## [0.22.0] - 2026-04-02

### Changed
- add smudgescope module
- Updated files: biopytools/smudgescope,docs/smudgescope.md,biopytools/cli/commands/smudgescope.py,biopytools/cli/main.py


## [0.21.1] - 2026-03-31

### Changed
- update README
- Updated files: README.md


## [0.21.0] - 2026-03-31

### Changed
- add deepbsa module
- Updated files: biopytools/deepbsa,docs/deepbsa.md,biopytools/cli/commands/deepbsa.py,biopytools/cli/main.py,README.md


## [0.20.0] - 2026-03-26

### Changed
- add hifi_hic module
- Updated files: biopytools/hifi_hic,docs/hifi_hic.md,biopytools/cli/commands/hifi_hic.py,biopytools/cli/main.py,README.md


## [0.19.0] - 2026-03-25

### Changed
- add agp2table module
- Updated files: biopytools/agp2table,docs/agp2table.md,biopytools/cli/commands/agp2table.py,biopytools/cli/main.py,README.md


## [0.18.0] - 2026-03-23

### Changed
- add vcf2pca module
- Updated files: biopytools/vcf_pca,docs/vcf2pca.md,biopytools/cli/commands/vcf2pca.py,biopytools/cli/main.py,README.md


## [0.17.1] - 2026-03-23

### Changed
- update README
- Updated files: README.md


## [0.17.0] - 2026-03-23

### Changed
- add dual rnaseq module
- Updated files: biopytools/dual_rnaseq,docs/dual_rnaseq.md,biopytools/cli/commands/dual_rnaseq.py,biopytools/cli/main.py


## [0.16.2] - 2026-03-15

### Changed
- update fastp module for --read1-suffix and --read2-suffix
- Updated files: biopytools/fastp,biopytools/cli/commands/fastp.py


## [0.16.1] - 2026-03-05

### Changed
- update genomescope module using conda env
- Updated files: biopytools/genome_analysis,biopytools/cli/commands/genome_analysis.py


## [0.16.0] - 2026-03-05

### Changed
- add genomescope module
- Updated files: biopytools/genome_analysis,docs/genomescope.md,biopytools/cli/commands/genome_analysis.py,biopytools/cli/main.py


## [0.15.2] - 2026-03-05

### Changed
- update BUSCO module
- Updated files: biopytools/busco,biopytools/cli/commands/busco.py


## [0.15.1] - 2026-02-27

### Changed
- update README
- Updated files: README.md


## [0.15.0] - 2026-02-27

### Changed
- add vcf2phylip module
- Updated files: biopytools/vcf2phylip,docs/vcf2phylip.md,biopytools/cli/commands/vcf2phylip.py,biopytools/cli/main.py


## [0.14.0] - 2026-02-24

### Changed
- add bwa module
- Updated files: biopytools/bwa,docs/bwa.md,biopytools/cli/commands/bwa.py,biopytools/cli/main.py


## [0.14.0] - 2026-02-24

### Changed
- add bwa module
- Updated files: biopytools/bwa,docs/bwa.md,biopytools/cli/commands/bwa.py,biopytools/cli/main.py


## [0.13.0] - 2026-02-24

### Changed
- add bam stat module
- Updated files: biopytools/bam_stats,docs/bam_coverage_stats.md,biopytools/cli/commands/bam_stats.py,biopytools/cli/main.py


## [0.12.1] - 2026-02-24

### Changed
- update README
- Updated files: README.md


## [0.12.0] - 2026-02-24

### Changed
- add conda env yml files
- Updated files: conda_env


## [0.11.0] - 2026-02-24

### Changed
- add conda env yml files
- Updated files: conda_env


## [0.10.0] - 2026-02-24

### Changed
- add iseq module and update README
- Updated files: biopytools/iseq,docs/iseq.md,biopytools/cli/commands/iseq.py,biopytools/cli/main.py,README.md


## [0.9.0] - 2026-02-24

### Changed
- add busco module
- Updated files: biopytools/busco,docs/busco.md,biopytools/cli/commands/busco.py,biopytools/cli/main.py,README.md,docs/busco.md


## [0.8.1] - 2026-02-09

### Changed
- update README
- Updated files: README.md


## [0.8.0] - 2026-02-09

### Changed
- add bam-cov module
- Updated files: biopytools/bam_cov,docs/bam_coverage_stats.md,biopytools/cli/commands/bam_cov.py,biopytools/cli/main.py


## [0.7.1] - 2026-02-05

### Changed
- update README
- Updated files: README.md


## [0.7.0] - 2026-02-05

### Changed
- add sra2fastq module
- Updated files: biopytools/sra2fastq,docs/sra2fastq.md,biopytools/cli/commands/sra2fastq.py,biopytools/cli/main.py


## [0.6.0] - 2026-02-05

### Changed
- add vcf2phylip module
- Updated files: biopytools/vcf2phylip,docs/vcf2phylip.md,biopytools/cli/commands/vcf2phylip.py,biopytools/cli/main.py


## [0.5.2] - 2026-02-05

### Changed
- fix version management


## [0.5.1] - 2026-02-05

### Changed
- fix version management


## [0.5.0] - 2026-01-24

### Changed
- add admixture module
- Updated files: biopytools/admixture,biopytools/cli/commands/admixture.py


## [0.4.6] - 2026-01-24

### Changed
- update fastp module
- Updated files: biopytools/fastp,biopytools/cli/commands/fastp.py


## [0.4.5] - 2026-01-24

### Changed
- upload rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py


## [0.4.4] - 2026-01-24

### Changed
- upload annovae module
- Updated files: biopytools/annovar,biopytools/cli/commands/annovar.py


## [0.4.3] - 2026-01-24

### Changed
- upload blast module
- Updated files: biopytools/blast,biopytools/cli/commands/blast.py


## [0.4.2] - 2026-01-13

### Changed
- update blast module
- Updated files: biopytools/blast,biopytools/cli/commands/rnaseq.py,biopytools/cli/commands/blast.py


## [0.4.1] - 2026-01-11

### Changed
- update blast annovar fastp rnaseq module
- Updated files: biopytools/annovar,biopytools/blast,biopytools/rnaseq,biopytools/fastp,biopytools/cli/commands/annovar.py,biopytools/cli/commands/blast.py,biopytools/cli/commands/rnaseq.py,biopytools/cli/commands/fastp


## [0.4.0] - 2026-01-08

### Changed
- add fastp module
- Updated files: biopytools/fastp,biopytools/cli/commands/fastp.py,biopytools/cli/main.py


## [0.3.5] - 2026-01-06

### Changed
- update rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py,./docs/rnaseq.md


## [0.3.5] - 2026-01-06

### Changed
- update rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py,./docs/rnaseq.md


## [0.3.5] - 2026-01-06

### Changed
- update rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py,./docs/rnaseq.md


## [0.3.4] - 2026-01-06

### Changed
- update blast module
- Updated files: biopytools/blast,biopytools/cli/commands/blast.py,README.md,./docs/blast_v2.md


## [0.3.3] - 2026-01-06

### Changed
- update annovar module
- Updated files: biopytools/annovar,biopytools/cli/main.py,biopytools/cli/commands/annovar.py,README.md,./docs/annovar.md


## [0.3.2] - 2026-01-03

### Changed
- update rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py,biopytools/cli/main.py,docs/rnaseq.md,README.md


## [0.3.1] - 2026-01-03

### Changed
- update README
- Updated files: README.md


## [0.3.0] - 2026-01-03

### Changed
- add rnaseq module
- Updated files: biopytools/rnaseq,biopytools/cli/commands/rnaseq.py,biopytools/cli/main.py,docs/rnaseq.md,README.md


## [0.2.5] - 2026-01-03

### Changed
- update blast cli and README
- Updated files: biopytools/cli/commands/blast.py,README.md


## [0.2.4] - 2026-01-03

### Changed
- update blast cli and README
- Updated files: biopytools/cli/commands/blast.py,README.md


## [0.2.3] - 2026-01-03

### Changed
- update blast cli and README
- Updated files: biopytools/cli/commands/blast.py,README.md


## [0.2.2] - 2026-01-03

### Changed
- update blast cli and README
- Updated files: biopytools/cli/commands/blast.py,README.md


## [0.2.2] - 2026-01-03

### Changed
- update blast cli and README
- Updated files: biopytools/cli/commands/blast.py,README.md


## [0.2.1] - 2026-01-03

### Changed
- update annovar module
- Updated files: biopytools/annovar,biopytools/cli/commands/annovar.py


## [0.2.0] - 2026-01-02

### Changed
- add annovar module
- Updated files: biopytools/annovar,biopytools/cli/commands/annovar.py,biopytools/cli/main.py,docs/annovar.md,README.md


## [0.1.0] - 2026-01-02

### Changed
- add blast module
- Updated files: biopytools/blast,biopytools/cli/commands/blast.py,biopytools/cli/main.py,docs/blast_v2.md,README.md


