
## [1.32.1] - 2026-08-14

### Added
- `insert2locus`（新模块，CLI `insert2locus`）：转基因插入位点提取——fastq 配对发现 → 参考比对步移（walking）→ 完整 locus 重建与验证 → 报告输出；配置与工具路径走 `expand_path`/`get_tool_path` 规范，含断点续传与 conda 包装调用

### Changed
- `mixrace`：核心方法论重构为单倍体方案（导师 v2）——变异检测由 bcftools mpileup+call 改为 `freebayes -p 1 --min-alternate-fraction 0.02 --min-coverage 30`（保低频等位，AF 用 AO/RO 字段）；判读引擎由四指标投票（VAF/Hw/多等位/Fws）改为全基因组杂合率 het_rate 主判据（导师 <0.01%/0.1%/1% 阈值）+ AFS 形态 + 优势株占比，Fws 弃用（疟原虫二倍体指标不适用）；流程 8 步改 7 步，新增 `--clean-fastq-dir`（跳过 QC）与 `--min-alt-fraction`，输出目录改 by-step（`02_alignment/`、`04_filtered/`、`05_vaf/` 等）
- `cim`：physical block 移除 est.map（v1.31.1 补的 est.map 在 r2=0.1 低 LD 标记上图谱膨胀至百万 cM、置换检验 13h+ 永不完成——有意回退），R 脚本显式标注 physical block 为 Mb pseudo-cM、与 mstmap block（真实 cM）尺度不同、LOD/阈值仅同 block 内可比；estimate 模式 est.map 保留
- `bwa`：help 示例精简为单条最简用法（§3.3）

### Fixed
- `mixrace` 15 项修复：repeat 过滤管道失败不再静默回退 raw（中止 + 清残留 + 不建断点）；无变异数据不再误判疑似纯（no_data → uncertain）；`software_versions.yml` 废弃参数改现行 13 个字段；HTML 报告 rationale/样名 `html.escape` 全路径转义；config 死字段清理 + step 注释 1..7；AFS 中间频率区间 0.05–0.95 真正生效（VAF≈0 噪声不再拉向 smeared）；`run_align` 顺序执行器失败短路；`_genome_size` 支持 .gz；`_parse_sn` 精确字段匹配；`--min-alt-fraction` argparse/Click 双层透传；`run_qc`/`run_kmer` 断点自愈（`_done` must_exist）；md 报告指标百分比格式化；测试同步更新 + 6 个回归测试

## [1.32.0] - 2026-08-12

### Added
- `mixrace`（新模块，CLI `mixrace`）：WGS 混合小种检测——fastq → bwa-mem2 比对（含正确 name-sort→fixmate -m→coord-sort→markdup 链，修复仓库既有 `bwa/alignment.py` 缺 fixmate 的 landmine）→ bcftools mpileup+call（保留 multiallelic + FORMAT/AD/DP）→ 过滤（QUAL/DP/SNP + 去 repeat + ALT reads 级联）→ VAF/Hw/多等位统计 + 队列 Fws + 可选 k-mer 谱（复用 `smudgescope`）→ 四指标透明投票判读（疑似纯/混合/灰色 + 置信 + 人话依据）。8 步 `--step` 调度 + 断点续传 + `--pure-samples` 已知纯样品阈值校准 + 跨样品汇总表（TSV/HTML）；conda 管道统一用单 `conda run -n env --no-capture-output bash -c 'pipe'`（同 env 工具一次激活，规避 §13.2.2）；`CommandRunner` 以独立进程组启动，SIGTERM/超时/Ctrl-C 整组杀（无孤儿 conda/samtools）
- `selective_sweep`：新增 XP-CLR 跨群体选择信号分量——每两群体×每染色体计算 XP-CLR 并入 composite_score；新增 `--xpclr-maxsnps`/`--minsnps`/`--ld`/`--min-samples`/`--include-xpclr-low-n`；filtered.vcf.gz 自动补建 `.tbi`（allel.read_vcf tabix 需要）；单染色体失败优雅降级（不阻塞其他染色体）；`software_versions.yml` 用 `python -c "import xpclr.__version__"` 探版本（xpclr 无 `--version` CLI）

### Changed
- `interproscan`：`build_command`/`CommandRunner.run` 由 str+`shell=True` 改 list+`shell=False`（§13 安全）；`config.__post_init__` 输入/输出路径改 `expand_path`（原 `os.path.abspath` 不展开 `~`，`~/input` 校验失败）；TSV 解析加行级容错（单行 int/score 畸形只跳过该行，不丢整份结果）；CLI `--format` 支持逗号分隔多格式（默认 `TSV,XML`）
- `resistify`：`CommandRunner.run` 改 list+`shell=False`；`--threads` 接线传给 resistify/nlrexpress（原未传）；日志持久化到 `{output_dir}/99_logs/resistify.log`（原 `ResistifyLogger()` 无参不落盘）；批量模式预创建样本子目录（避免 resistify 写文件失败）；删除未使用的 `parse_annotations` 死代码
- `nlr_annotator`：`build_command` 改用 `build_conda_command` 包装 java（原裸 `['java', ...]` 在 conda env 调不到正确解释器，属死代码激活），新增 `--java-path`；批量模式 `NLRLogger` 改每样本独立命名 logger（原共享 `"nlr_annotator"` logger 被后续样本劫持 handler）；`merge_only`/`collect_input_files` 异常转 `sys.exit(1)`；CLI 透传 jar/mot/store/java 路径与三项距离参数

### Fixed
- `interproscan`：`--goterms`/`--pathways` 死开关 bug——原反向开关 `--no-goterms`（`action='store_true'`, `default=True`）恒为 True，`goterms = args.no_goterms is False` 永远 False，GO/Pathway 从不获取；改读正向 `args.goterms`/`args.pathways`；GO 类别匹配从缩写 `|BP|`/`|CC|`/`|MF|` 改全称 `|Biological Process|` 等（匹配 TSV 实际填充格式）
- `resistify`：序列提取文件（nlr.fasta 等）预检条件收紧为仅 `--skip-resistify` 模式；原 `if not self.is_directory` 在 non-skip 模式会检查尚未由 resistify 生成的文件而误报缺失
- `gtx_joint`：移除参考基因组 `.fai` 索引强制预检（索引由工具自建，原校验误拦合法输入）

## [1.31.1] - 2026-08-12

### Fixed
- `cim`：修复 physical block 下 CIM 的 LOD 整体压低问题。`build_tidy_files` 写入的 Mb 占位位置原注释"est.map will override"，但 mstmap 模式从未调用 `est.map`，导致 Mb 被当成 cM 喂给 `cim()`——`window=10` 被读成 10Mb，协因子剔除窗口占整条染色体约 24%，引发协因子饥饿、LOD 被整体压低。现对 physical block 补 `est.map()` + `replace.map()`，用真实重组率重估 cM，使 window/step 与 mstmap block 处于同一真实 cM 尺度
- `vcf_sample_hete`：`config.__post_init__` 改用 `expand_path` 展开 `~`/`$VAR`（原 `os.path.abspath` 不展开 `~`，`~/input.vcf` 在后续 `os.path.exists` 校验时失败）；`__init__.py` 文档示例导入路径误写 `biopytools.vcf_stats` → 更正为 `biopytools.vcf_sample_hete`

### Changed
- `vcf_sample_hete`：删除从未生效的死参数 `--dry-run`/`force`（CLI/config/main 全链路移除，原参数从未被处理逻辑读取）；日志输出落 `99_logs/` 子目录（§12）；`VCFStatsLogger` 支持显式 `--log-file`/`--log-level`，级别优先级 log_level > quiet > verbose > INFO，并加 `propagate=False`；大数字用 `format_number`（>1M 显示为 M 单位，§5.3）；生成 `00_pipeline_info/software_versions.yml`（纯文本手写、无 yaml 依赖）；DEBUG 级补输出过滤诊断（qual/depth 丢弃数）与逐样本统计

## [1.31.0] - 2026-08-12

### Added
- `centier`：新增 Hi-C FASTQ 自动模式——提供 `-1`/`-2` 即自动跑 HiC-Pro 产出 100kb/20kb 矩阵再进 CentIER 全流程（复用 `hic_heatmap.HiCProPipeline.run_hicpro_only()`，`require_plothic=False`）；新增 ChrN 命名预检（CentIER 实测要求 id 以 `ChrN` 结尾，`--strict-chrname` 不符即中止，默认仅 WARNING）；Hi-C 模式 CentIER 结果落 `02_centier/`（HiC-Pro 输出在 `01_hic_mapping/`）；生成 `00_pipeline_info/software_versions.yml`；运行后清理 CentIER 在 CWD 生成的临时文件到 `tmp/`；CLI 新增 `-1`/`-2`/`-g/--genome-id`/`--restriction-enzyme`/`--bowtie2-idx`/`--bin-sizes`/`--max-memory`/`--force-hicpro`/`--hic-matrix-type`/`--strict-chrname`

### Changed
- `centier`：`main()` 改用 `CentIERRunner` 统一管理 config/logger/analyzer（kwargs 转发）；CentIER 执行命令补 `命令|Command:` INFO 日志（§2.2.1）；`traceback.print_exc()` 改 `logger.error(traceback.format_exc())`；`--step` 类型 `Choice(['1'..'6'])`→`IntRange(1,6)`；`run_command` 命令日志 DEBUG→INFO
- `centier`：依赖检查增强——CentIER 自带二进制补可执行权限检查、`gt`/`LTR_retriever` 改 PATH 搜索（原硬编码 `~/miniforge3/envs/centier/bin/`）+ conda env 兜底安装建议、Hi-C 自动模式额外校验 HiC-Pro/bowtie2-build、Python 包（pyfastx/numpy/pandas/scipy）存在性检查

## [1.30.0] - 2026-08-12

### Added
- `admixture`：新增 `--dry-run` 模拟运行模式（打印各步骤计划命令但不执行）；分层输出目录（`00_pipeline_info`/`01_preprocessing`/`02_plink`/`03_admixture`/`04_results`/`99_logs`）；per-K 与各预处理/质控步骤支持断点续传；生成 `00_pipeline_info/software_versions.yml`；CLI 补全 LD 剪枝参数（`--ld-prune`/`--ld-window`/`--ld-step`/`--ld-r2`）
- `hic_heatmap`：拆分 `run_hicpro_only()`（仅跑 HiC-Pro 不画热图）供其他模块复用；新增 `require_plothic` 配置（仅在画热图时校验 plothic，默认 True 向后兼容）

### Changed
- `admixture` 全面合规化重构：所有 plink/bcftools/admixture/adamixture 调用改走 `build_conda_command` 传**完整工具路径**（修复原 `build_conda_command('admixture', ...)` 裸命令名、以及 `os.path.basename(adamixture_path)` 提取命令名的 §13.6 违规）；管道命令（`zcat|grep|cut|sort|uniq`）改 Python `Counter` 实现，规避 conda-run-in-pipe（§13.2.1）；复制 bed/fam 改 `shutil`；`AdmixtureLogger` 改命名 logger + `propagate=False`（隔离 root 避免污染全局）；依赖检查改非阻断（仅 warning，登录节点 PATH 不可靠）；R 可视化缺包则跳过不自动装（超算无网）、读 Q/fam/cv 改绝对路径；工具路径走 `get_tool_path`；threads 默认 64→12
- `vcf2tree`：`run_pipeline` 异常改 re-raise（不再 `sys.exit`），由 CLI 层 catch 转退出码，保持库可用性；`get_conda_env` 加绝对路径特判（避免 `os.path.join` 塌缩误判静态二进制所属环境）；FastTree 后端加 WARNING（不支持 ASC 校正，SNP 分支长度可能低估，需 ASC 请用 IQ-TREE）；入口提示"不做 LD 剪枝/位点 QC，假设输入已质控"；`vcf_converter` 统计并记录跳过的 InDel/非单碱基位点数
- `cim`：CIM R 块剔除孤立单标记 LG（<3 markers），避免抬高置换检验阈值

### Fixed
- `admixture`：`config.validate` 补 `ld_step > 0` 校验

## [1.29.0] - 2026-08-11

### Added
- `selective_sweep`（新模块，CLI `selective-sweep`）：选择性扫荡检测——VCF 过滤（双等位 SNP + MAF + 缺失率）→ 群体拆分 → 每群体 π / Tajima's D / RAiSD μ / SweeD CLR + 两两 Weir-Cockerham Fst → composite_score（各分量经验百分位均值）→ top 分位数候选窗口合并输出。支持断点续传、低样本群体默认排除 MU/CLR 分量（`--include-mu-low-n`/`--include-sweed-low-n` 可强制加入）、`.gz` 输入、`00_pipeline_info/software_versions.yml`、`tmp/` 中间文件清理

### Changed
- `cim`：`CommandRunner` 由 `subprocess.run` 改为 `Popen` + `start_new_session=True`，并注册 SIGTERM 处理；用户中断（Ctrl-C）/ 超时 / 作业调度取消（SIGTERM）时整组杀掉子进程（含 sh/conda/Rscript 全部后代），消除超算孤儿进程

### Fixed
- `deepbsa`：`vcf2csv` 支持 `.gz` 输入（`gzip.open(..., 'rt')`）；修复 `VCF` 迭代器首条记录丢失——`__init__` 已预读首条 `record`，原 `__next__` 首次调用直接读下一行致首条遗漏，加 `_first_yielded` 标志先返回已读记录
- `vcf2tree`：`config.__post_init__` 路径标准化（`abspath`）提前到 `mkdir` 之前——原顺序下相对路径 `-o` 产生的 `snps_fa` 仍为相对路径，FastTree 以 `step2_dir` 为 cwd 时找不到文件

## [1.28.2] - 2026-08-10

### Removed
- `cim`：移除 `fixGenoError` 基因型纠错步骤——删除 `fix_geno_error` 函数、`config.fix_geno_error_size` 字段、`--fix-geno-error-size` CLI/argparse 选项、`Step 3.5` 调用块。该步骤原默认开启（fix_size=10），移除后 CIM 流程不再做 RLE 短片段基因型纠错，标记基因型与下游 CIM 结果可能略有差异

## [1.28.1] - 2026-08-09

### Fixed
- `vcf2tree`：FastTree 后端改走 `build_conda_command`（默认 `~/.local/bin/FastTree`，原裸 `'fasttree'` 靠 PATH，conda env 未激活则 `FileNotFoundError`）；移除 stdout 空时误读 stderr 当 Newick 的 fallback（FastTree 从不把树写到 stderr）
- `vcf2tree`：`get_iupac_code` 入口统一 `ref/alt.upper()`（原仅大写键，小写 ref/alt 被误判为缺失 N）；序列存储 `list`→`array('u')`（50 样本×1M SNP 内存从 ~2.5GB 降到 ~200MB）
- `vcf2tree`：IQ-TREE 默认开启 SNP 数据 ASC 校正（`-m MFP+ASC`，Lewis 2001，避免分支长度低估），加 `--no-asc` 开关；补 `00_pipeline_info/software_versions.yml`；`expand_path` 改用 `common.paths`（删本地实现）
- `cim`：MSTmap 自动调优全部 p.value 无标记时 `best_pvalue=None` → `f"{None:.0e}"` TypeError 崩溃（v1.28.0 引入的回归），改为清晰 `sys.exit(1)` 报错
- `filter_snp_indel`：`FilterLogger` 补 `propagate=False`（避免向 root 重复输出，既有问题）

## [1.28.0] - 2026-08-09

### Added
- `vcf2tree`（新模块，CLI `vcf2tree`）：VCF→系统发育树——SNP 转 IUPAC FASTA 对齐矩阵，再用 FastTree 或 IQ-TREE2 建树（默认 IQ-TREE2，ModelFinder + UFBoot 1000）；两步流程 + 断点续传。**注**：FastTree 后端尚未走 conda 包装（默认 IQ-TREE2 可用）
- `cim`：`build_mstmap_linkage_map` 新增"核心 LG"判定（≥3 markers，排除孤立单标记）+ "最少核心 LG 优先"兜底

### Changed
- `cim`：`_calc_pairwise_rf_matrix` 由 O(n²) 循环重构为矩阵向量化；RF 均值分母从固定 `n-1` 改为"实际有效配对数"（修 bug，仅当存在无共同样本的标记对时结果不同）；plink/bcftools 改 `build_conda_command` 包装
- `filter_snp_indel`：所有 bcftools 调用 `shell=True`→`shell=False`+conda 包装；删死码（注释旧方法 + 零调用方方法，过滤逻辑逐字保留）；`bcftools_path` 走 `get_tool_path`；默认线程 64→12

### Fixed
- `cim`：R/qtl `est.map` bug——`est.map()` 返回 list，旧代码 `cross <- est.map(...)` 覆盖 cross 对象致后续 `calc.genoprob/cim` 崩溃，改 `replace.map(cross, newmap)`（`map_mode=estimate` 此前基本跑不通）
- `cim`：`keep_sample_mask` 应用时机从 `save_filtered_vcf` 之后移到之前（旧顺序保存的过滤 VCF 带已丢弃的全缺失样本）
- `filter_snp_indel`：`--no-snp-biallelic` 帮助文案中英文自相矛盾（"不禁用" vs "Do not"）→ 改为"禁用|Disable"；删从未接通到 config 的死参数 `--log-file`

## [1.27.1] - 2026-08-06

### Fixed
- `cli`（全仓 35 个命令）：`_validate_path_exists` 加 `os.path.expanduser`，用户传 `~/...` 输入路径不再在 Click 校验阶段被提前拒（`bam2fastq`/`fof` 用 `lexists` 支持软链接）
- `fof`：`extract_sample_name` 双扩展剥层修复（`.fq.gz`→`sample1_R1`，原只剥 `.gz` 留 `.fq`）；`expand_path` 导入改 `from ..common.paths`（原顶层导入靠 fallback 兜底）；`FofLogger` 改标准 3-way（named logger + `propagate=False` + stdout/stderr/file）
- `genome2sv`：`parse_fof` 展开 fof 内路径的 `~`/环境变量（原 validate 存在性检查误报）；`check_dependencies` 对 SURVIVOR（子命令式 CLI 无 `--version`）按"能执行"判定而非 `rc==0`，避免误报缺失终止流程
- `vcf2pca`：`_drop_zero_genotype_samples` 的 `--remove` 改绝对路径（`.resolve()`，原相对路径在 PLINK `cwd` 下找不到→降级跳过剔除→`--pca` 仍崩）；删除一段不可达死代码
- `assembly_qc`：NGS 单文件输入不再回退父目录 glob 全目录（避免误收无关样本），改为只取推断的 read2，找不到则警告返回；`docs/assembly_qc.md` 同步 `--skip-lai`→`--enable-lai`
- `psvcp`：`9gff_split.py` 加空 gff 守卫（`EmptyDataError`/`gff.empty`/`not out_dfs` 三重，对齐 R 原版空输出）

## [1.27.0] - 2026-08-06

### Added
- `fof`（新模块，CLI `fof`）：FOF 文件生成——扫描文件/目录生成 `sample<TAB>绝对路径` 映射清单，支持 `-s` 后缀过滤、`-r` 递归、单文件输入
- `genome2sv`（新模块，CLI `genome2sv`）：assembly-to-assembly SV calling——minimap2 → samtools 排序索引 → svim-asm → SURVIVOR 合并多样本 VCF → bcftools stats + SVTYPE 统计；支持断点续传、by-step 输出目录、`software_versions.yml`
- `vcf2pav`（新模块，CLI `vcf2pav`）：将 SURVIVOR 合并后的多样本 VCF 转 PAV（Presence/Absence）矩阵 TSV（行=SV，列=样本，值=0/1/NaN）+ 每样本统计摘要；纯 Python 解析
- `psvcp`：新增 `pav_table.py`——从 `pan.pav.gff` 第 9 列解析来源样本，生成 `pan.pav.info.tsv` + `pan.pav.matrix.tsv`（`_finalize` 惰性导入，失败降级非致命）；R 步骤 8/9/10 的 Python 向量化端口（`8gff_update.py`/`9gff_split.py`/`10gene_in_pv.py`，`np.searchsorted`），pipeline 改 `_run_py` 调用，R 版保留作参考

### Changed
- `vcf2pca`：PLINK 扶正为默认后端（`v2p`→`plink`，支持 SNP/INDEL）；`--skip-qc`（默认过滤）→ `--apply-qc`（默认不过滤）；默认不过滤时自动剔除 `F_MISS==1.0` 零基因型样本（避免 PLINK `--pca` 崩溃）；新增 `config.base_name` 输出前缀、`-g/--group-column` 分组着色（贯通 visualization）、合并 Excel `pca_results.xlsx`
- `assembly_qc`：`--lai-repeatmasker-species`（默认 `Viridiplantae`），EDTA 回退路径给 RepeatMasker 显式 `-species`；`ngs_reads`/`long_reads` 放宽单文件输入

### Fixed
- `psvcp`：R 脚本 9/10/12 `read.table` 加 `sep='\t', quote=''`，修 gff3 属性含空格被默认空白切分导致列错位；R9 补空 chr 守卫（query-only chr 触发 `seq(0)` 索引 0 行 → NA → 中断）

### ⚠️ BREAKING
- `assembly_qc`：`skip_lai` 默认 `False`→`True`（CLI 改 `--enable-lai` opt-in，因 EDTA 耗时长）。**旧脚本不带任何 LAI flag 时，原会跑 LAI，现会跳过**——依赖默认跑 LAI 的自动化脚本需显式加 `--enable-lai`
- `cli`：`vcf-pca` + `vcf2pca` 合并为单条 `vcf2pca`（独立 `vcf-pca` 入口移除）

## [1.26.0] - 2026-08-05

### Added
- `seq_len`（新模块，CLI `seq-len`）：FASTA 序列长度统计——输入单个文件或文件夹，流式计算每条序列长度（不载入整条序列，支持 `.gz`），输出主表 TSV；文件夹模式合并所有文件并加 `source_file` 列，同时输出含 N50/总长/最短/最长/平均 的汇总表（每文件一行 + ALL 全局行）。智能 `-o`（目录→`{prefix}.seq_len.tsv` / 文件→直接当主表）、`--min-length` 过滤、`--sort` 降序
- `nlr_annotator`：新增 `--merge-only`——批处理中途被杀（超时/OOM）未生成汇总时，直接合并已有 `*.nlr_annotator.tsv` 为 `nlr_annotator_summary.tsv`，不调用 Java、不重跑 NLR-Annotator；`collect_result_files` 自动识别 by-sample 子目录与平铺两种形态（按 abspath 去重、按样本名排序）
- `psvcp`：基因组注释支持 `.gff3`——新增 `PSVCPConfig.find_genome_gff`（`.gff` > `.gff3` 顺序查找），`validate()` 与 pipeline 统一改用（原 `_gff_name` 只认 `.gff`）

## [1.25.1] - 2026-08-04

### Fixed
- `kmertools` extract：unikmer 路径统一交由 `ExtractConfig.__post_init__` 解析（环境变量 > 配置文件 > 默认，再展开 `~`），移除 `cmd_extract` 里 `shutil.which` 的重复/不可靠查找（未在 PATH 时会返回 None 导致回退行为不一致）
- `vcf2pca` V2P 后端：`-InVCF`/`-OutPut` 改用绝对路径（原相对路径在 `subprocess.run(cwd=output_dir)` 下会被错误解析）；未请求聚类时显式传 `-ClusterMethod None`（工具默认 EM 聚类 ClusterMethod=3，数据含 NaN 时会在 EM 阶段崩溃）

## [1.25.0] - 2026-08-03

### Added
- **新模块 `needle_identity`**：EMBOSS needle 序列两两 identity 计算（`biopytools needle-identity`）——输入 FASTA 内所有序列对并行比对，输出 `{prefix}.needle_identity.tsv`（identity/similarity/gaps/score），`--gapopen`/`--gapextend`/`--threads` 可调；conda 自动包装、`99_logs/`+`software_versions.yml`、临时文件落 `output_dir/tmp` 并清理
- `gene_table`：新增 `--upstream`(默认3000)/`--downstream`(默认1000) 侧翼长度参数，输出表新增 `Region` 列 + `{prefix}.region.fa`（链定向 上游+基因+下游 区间，越界夹紧并警告）
- `annorefine`：小蛋白通道新增**强同源直通**——identity ≥ `--small-strong-homology-identity`(默认95%) 的命中视为铁证（近乎自比对），绕过 TE 排除与表达过滤（效应子常在 TE 区且低表达，不该被辅助证据拦）；新增 `--small-min-expression-depth`/`--small-min-coverage-breadth`/`--no-small-exclude-te` CLI 参数
- `phyto_effector`：Protease HMM 从 paper 版扩为 203 个 Pfam peptidase 家族（`protease_pfam.hmm`），命中面变广后加 score≥20 + 必须有信号肽双重控误报；CLI 改为单命令 + 目录输入自动逐样本鉴定并合并结果（原独立 `merge` 子命令移除）
- `msaviz`：MAFFT 默认参数加 `--preservecase`；颜色方案默认自动选择（DNA→Nucleotide，蛋白→Zappo），配色查询与 DNA/蛋白判定大小写不敏感
- `nlr_annotator`：支持 DNA/CDS 输入，目录模式默认匹配后缀 `*.cds.fa` → `*.fa`
- 新增 `conda_env/psvcp_v.1.0.1.yml`（PSVCP v1.0.1 完整环境锁定文件）

### Changed
- `phyto_effector`：SignalP 3/6/both 逻辑统一收口到 `run_signalp_pipeline`（原 main/rxlr/crn 三处重复实现，存在版本分支覆盖、`both` 模式尾部残留死代码等问题）

### Fixed
- `phyto_effector`：`parse_domtblout_details` 域 E-value/score 列序修正——HMMER3 domtblout 第 13/14 列才是 i-Evalue/domain score（原误取 12/13 列的 c-Evalue/i-Evalue）
- `phyto_effector`：TMHMM 改 `build_conda_command` 包装 + stdin 传文件句柄（原裸路径 + 整文件读内存）；`makeblastdb` 直接由 blastp 同 bin 目录推导；blastp/tmhmm 结果解析不再按 `/` 截断序列 ID
- `blast` CLI：修复 `--no-auto-detect-samples` 不可达问题（现在仅在关闭时透传该开关）

## [1.24.3] - 2026-07-30

### Added
- `tgsgapcloser`：新增 gap 填充前后对比报告——对比处理前(`scaff_file`)与处理后(`gapclosed`)的 gap 数量/总长/长度分布/top10/填充率，产出 `{prefix}.gap_report.txt`（人类可读报告）与 `{prefix}.gap_table.tsv`（per-gap 明细：处理前长度 + 填充状态 Filled/Remaining/Unknown + 处理后残留长度）；`dry_run` 或无最终输出时自动跳过。`quartet_filler` 新增 `analyze_gaps`/`track_gaps` 支撑方法（基于 flanking 锚点定位，重复区非唯一锚点极少数 gap 可能误判）

## [1.24.2] - 2026-07-30

### Added
- `annorefine`：新增**小蛋白回收通道**（`--recover-small-proteins`，默认关）——`gap_min_cds_len` 硬悬崖会把小蛋白整类丢掉，本通道放宽长度（CDS < `gap_min_cds_len` 且 ≤ 450bp/150aa），用 完整ORF① + 同源 + 表达②⑤ + 强制 TE 排除 四重门控找回；无表达数据时退化为 ORF + 严同源(70/80)。输出独立 ID 前缀 `{prefix}_small_gap_{N}`
- `tgsgapcloser` / `gap-fill` CLI：新增 `-f/--force`（忽略断点续传）、`--dry-run`（只打印命令）、`-mgl/--min-gap-length`（第2轮最小 gap 长度，原硬编码 100）；工具路径走路径管理系统（§11）

### Changed
- `annorefine`：`detect_gaps` / `detect_merged_genes` 性能优化——按染色体建基因 span / 命中 start 的 bisect 索引，把 O(hits×genes) 暴力扫描降为对数级；结果与原实现完全一致（CDS⊆span 假设下，span 相交预过滤 + 精确 CDS 判定），命中输出顺序还原保证确定性
- `annorefine`：写 `00_pipeline_info/software_versions.yml`（miniprot/samtools/stringtie 版本 + 关键参数，§12.5）

### Fixed
- `tgsgapcloser`：logger 改 named logger + `propagate=False`（原 `logging.basicConfig` 配置 root logger，biopytools 作为库被 import 时会与其他模块串扰/重复输出，§2.3.3）
- `tgsgapcloser`：第1轮 TGS-GapCloser2 改流式执行（`run_with_progress`，Popen+读线程），原 `capture_output` 把大基因组输出全缓冲进内存致 OOM（§13.2.0）
- `tgsgapcloser`：`minimap2` 改 `build_conda_command` + `shell=False`（原 `subprocess.run(shell=True)`），记完整命令（§2.2.1）；`unitig_dict` 索引加 None 防护（原直接索引会 KeyError）
- `tgsgapcloser`：`__post_init__` 展开所有工具路径（原漏了 unitig/ngs/racon/pilon/samtools/java 的 `~` 展开，§11.3.1）；`read_fasta` 改逐行解析省内存 + 兼容 `\r\n`
- `gap-fill` CLI：`tgsgapcloser_path` 改为始终透传（修 env/config.yml 探测的默认值到不了 main 的 bug）

## [1.24.1] - 2026-07-30

### Fixed
- `rnaseq2vcf`：`__post_init__` 先 `expand_path(output_dir)` 再 `mkdir`（§11.3.1，原顺序对 `-o ~/x` 会建出字面 `~` 目录）
- `rnaseq2vcf`：hisat2|samtools 管道加 `set -o pipefail`——原写法 hisat2 崩溃时 samtools 仍返回 0，静默产出残缺 BAM；所有路径 `shlex.quote` 化（空格/特殊字符安全）；支持 `env=None` 回退直调（非 conda 安装的 hisat2）
- `rnaseq2vcf`：calling 断点续传须同时存在 gVCF + `.tbi`（HaplotypeCaller 正常退出必产两者；仅 `.g.vcf.gz` 缺 `.tbi` 是中断残缺产物，不再误判已完成）
- `rnaseq2vcf`：`--force` 现绕过共享 `genome_index` 的 checkpoint（原对共享索引无效）
- `rnaseq2vcf`：变异计数（`bcftools view -H`）检查返回码，失败返回 None 让报告显示「不可用」，而非不可信的部分计数
- `rnaseq2vcf`：GATK 临时文件落 `output_dir/tmp`（`--java-options -Djava.io.tmpdir`，§12.4.1，避免超算 `/tmp` 爆满）

### Changed
- `rnaseq2vcf` CLI：暴露 GATK VariantFiltration 阈值 `--fs-threshold`(30)/`--qd-threshold`(2)/`--cluster-window`(35)/`--cluster-size`(3)（原硬编码）；`--gff3` 改可选（HISAT2 可不带已知剪接位点 de novo 运行）；`--step` 简化为仅 `0`(仅建索引)/省略(全流程)

## [1.24.0] - 2026-07-29

### Added
- **新模块 `predgpi`**：PredGPI GPI 锚定蛋白预测。Python-import 方式调用 PredGPI（`sys.path` + `PREDGPI_HOME`，无需 conda 包装），延迟加载 HMM/SVM 模型；逐序列预测含非标准氨基酸替换 + 短序列(≤40aa)防护 + **逐序列异常兜底**（predgpi 内部 numpy/HMM 报错时标记为非 GPI，不中断整批）+ numpy 标量强转（`float()`/`int()` 避免 truth-value 歧义）；FPR 四级分类（Highly/Probable/Weakly/Improbable）；断点续传；3-handler 日志分离；`software_versions.yml`（PredGPI 无 `--version`，以发表文献标注版本）。CLI：`biopytools predgpi -i proteins.fa -o out/ [--conservative] [--predgpi-home ~/software/predgpi]`

## [1.23.3] - 2026-07-29

### Added
- `gff_renamer`：新增 `-v/--version <tag>` 选项——非空时在所有生成的基因 ID 最前面加 `v{tag}-` 前缀（如 `-v 1` → `v1-CDRT_Ov12g000010`），下游 mRNA/exon/CDS 派生 ID 自动继承。`config.validate()` 校验版本号不含空白字符（避免破坏 GFF 属性 ID）

### Changed
- `gff_renamer`：`--version` 语义变更——原为 argparse 内置 `action='version'`（打印程序版本号 `1.0.0`），现改为接收版本标签参数用于 ID 前缀；依赖 `--version` 查程序版本的用法不再生效

## [1.23.2] - 2026-07-28

### Docs
- CLAUDE.md v2.17 → **v2.18**：新增 §10.4「临时测试/调试输出位置」——探索性/ad-hoc 测试（跑函数看输出、临时脚本验证、肉眼检查）产物**严禁写入仓库 cwd**，统一放 `~/tmp/<描述性子目录>/`；与 §11.A 正式 pytest 单测（`tests/`）、§12.4 模块运行时临时文件（`output_dir/tmp`）明确区分，避免污染仓库/误提交/`.gitignore` 漏网

## [1.23.1] - 2026-07-28

### Fixed
- `psvcp`：补提交 1.23.0 漏掉的 `psvcp/scripts/construct_pan/` 下 **13 个 vendored R/Python helper**（`pipeline.py` 逐轮调用它们，缺失则模块在 fresh checkout 上不可用）。根因：`.gitignore` 的 `scripts/` 规则过宽（匹配任意深度），已改为锚定根目录 `/scripts/`（根 `scripts/` 仍忽略，模块内嵌 `scripts/` 可跟踪）

## [1.23.0] - 2026-07-28

### Added
- **新模块 `psvcp`**：PSVCP 线性泛基因组构建（MUMmer + Assemblytics 检测 PAV 并逐轮并入参考）。Python 重做原始 PSVCP 脚本的逐轮循环（步骤序列与命名 1:1 对齐原始），外部工具走 `build_conda_command` + 命令日志，vendored R/Python helper 原样调用；整体 + 逐轮断点续传；assemblytics numpy 预检（env python 被 GraalPy 顶掉时明确报错而非黑箱失败）；写 `software_versions.yml`。CLI：`biopytools psvcp -i genome_gff_dir/ -l genome_list.txt -o out/`

### Changed
- `annorefine` 升级为**端到端**：原仅做 BRAKER 后查漏补缺，现一条命令跑完 BRAKER 注释 + 同源查漏补缺 → 整合 GFF3（吸收原 braker4ps 的角色）；CLI 参数相应扩展（`-s species` / `--rnaseq-dirs` / `--fungus` / `--skip-*` 等）
- **弃用 `braker4ps`**：功能已完全并入 annorefine（两者均为 `BrakerPipeline + AnnorefineRunner` 端到端，100% 重复），从 CLI 注册表移除并删除模块 + CLI + 文档。**破坏性|Breaking**：原 `biopytools braker4ps` 命令不再存在，改用 `biopytools annorefine`

### Fixed
- `annorefine`：表达深度计算失败时返回 `None` → `qc_filter` 跳过表达过滤（原返回 `{}` 会让所有基因被判无表达而**静默全丢**）；hit 索引从 `id(hit)` 改为内容键 `hit_key`（含 `cds_exons`，防 miniprot 同 span 多 CDS 结构碰撞，survives 对象复制）；每 hit 在自身 CDS 区间独立求深度（消除重叠 CDS 的归属不确定性）；解析 `braker.gff3` 跳过畸形坐标行不崩

## [1.22.0] - 2026-07-28

### Changed
- **模块重命名 `ps_gene_anno` → `annorefine`**：原定位 BRAKER 后「效应子查漏补缺」(疫霉专用)，现重新定位为通用「注释精修」(同源补漏 + 合并拆分 + ORF/表达质控，适用任意物种)
  - 类名 `PsGeneAnnoConfig/Runner/Logger` → `AnnorefineConfig/Runner/Logger`
  - CLI 命令 `biopytools ps-gene-anno` → `biopytools annorefine`
  - `braker4ps` 改为包装 annorefine（braker4ps 命令名不变）
  - `docs/ps_gene_anno.md` → `docs/annorefine.md`
  - 逻辑无变化（v1.21.8 的真实 ORF / 表达证据 / 坐标零重叠 / 合并拆分门控等完整保留）
  - **破坏性|Breaking**：原 `ps-gene-anno` 命令不再存在，改用 `annorefine`
  - 注：超算上若仍残留旧 `biopytools/ps_gene_anno/` 目录（`copybiopytools` 无 `--delete` 不会自动清），需手动 `rm -rf biopytools/ps_gene_anno biopytools/cli/commands/ps_gene_anno.py`，否则下次同步又会推回本地

## [1.21.8] - 2026-07-27

### Added
- `ps_gene_anno`：查漏补缺新增多层生物学质控——① **真实完整 ORF**（翻译验证：CDS 3 倍数 + ATG + 无内部终止 + 末终止，比 miniprot 覆盖率更严）；②⑤ **表达证据**（新 `expression.py`：用唯一比对 reads `NH==1`（旧版 samtools 回退 MAPQ）过滤 BAM，算每 hit 平均深度 + 覆盖广度；多比对 reads 多落 TE/重复区会虚假抬升表达故排除）；③ gap **坐标零重叠**模式（与任一 BRAKER 基因 span 有坐标交集即不算新基因）；路径分治 `enable_gap_fill`/`enable_split`；★ **合并拆分门控**（某 merged 基因的 split copies 全未过 QC 则回退保留原 BRAKER 基因，不误删）
- `braker4ps` CLI：新增 `--no-real-orf`/`--no-coord-zero-overlap`/`--no-unique-reads`/`--min-unique-mapq`/`--min-expression-depth`/`--min-coverage-breadth`/`--no-gap-fill`（透传 ps_gene_anno 新参数）
- `swave`：call 步骤断点续传（主 VCF 存在则跳过）；`SWAVE_CONDA_ENV` 环境变量覆盖 conda 环境名（默认 swave_v.1.2）；上游 swave 跳过 snarl 的 dotplot 时扫描 swave.log 以 WARNING 提示（数据完整性透明）；passthrough 子命令加命令日志
- `cactus`：默认输出新增 `odgi`；`get_genome_dirs` 解析 seqfile 挂载所有基因组文件所在目录（基因组分散多目录时容器内 cactus 才能读取）

### Fixed
- `swave`：找不到 conda 环境时改为 **明确 RuntimeError**（原静默回退系统 python 会因缺 pysam 等崩溃、报错极不直观）
- `cactus`：`bind_paths` **缩进 bug**——seqfile/输出目录的绑定原误嵌在 `if work_dir not in bind_paths` 块内（仅新增 work_dir 时才执行），改为始终执行
- `cactus`：参考基因组必须是 seqfile 首条（minigraph-cactus 硬性要求 `--reference` 与首条匹配），不符则明确报错
- `braker4ps`：argparse help 中字面 `%` 转义为 `%%`（原被当作格式符）

## [1.21.7] - 2026-07-27

### Added
- `annovar`：外显子变异结果表与全变异结果表新增 **`VCF坐标` 列**——回填每条 ANNOVAR 记录对应的原始 VCF POS。新增 `VcfCoordinateMapper`：读原始 VCF（支持 `.gz`），按「归一化等位基因」（trim REF/ALT 共同前后缀，消除 `convert2annovar.pl` 对 indel 锚碱基归一化导致的 1bp 偏差）建索引，反查 ANNOVAR 记录的原始 VCF POS（同键多位点取最近）；懒加载、两 processor 共享、读取失败仅告警且列填 `NA`
- `annovar`：`ANNOVARResultsProcessor` 新增 `vcf_file` 参数，由 `main.py` 透传 `config.vcf_file`

### Changed
- `ps_gene_anno`：GFF3 source 列（第 2 列）`ps_gene_anno` → `psfill`（gene/mRNA/exon/CDS 行一致）

## [1.21.6] - 2026-07-27

### Added
- `annovar`：蛋白序列与 DNA 变异效果处理支持 **dup（重复）** 变异型——`ProteinSeqModifier` / `ExonicVariantProcessor` 新增 `c.N_Mdup`（区间重复）/ `c.Ndup`（单碱基重复）解析（ANNOVAR cDNA 记法），在重复区间后插入被复制的碱基（无显式碱基时复制源区间）；此前 dup 变异落入兜底分支未正确应用

## [1.21.5] - 2026-07-27

### Added
- `rnaseq2vcf`：分析报告 `ANALYSIS_REPORT.txt` 增加变异过滤统计——过滤前/后变异数、被滤掉占比（FS>30 / QD<2 / cluster 3@35bp）、每样本 PASS 非参考基因型计数（解析 `bcftools stats` 的 PSC 行：nNonRefHom+nHet）；并标注下游送 annovar 的最终 VCF 及其索引/被滤变异复查路径。`bcftools view -H` 流式计数（内存安全），失败仅告警不阻断
- `conda_env/gene-anno.yml`：新增 gene-anno 注释环境定义

## [1.21.4] - 2026-07-25

### Changed
- `rnaseq2vcf`：变异检测从**单样本 VCF** 改为 **多样本联合 calling**——每样本 `HaplotypeCaller -ERC GVCF` 生成 gVCF，再 `CombineGVCFs → GenotypeGVCFs → VariantFiltration → bcftools PASS` 合并为一个多样本 VCF（GATK 最佳实践）。`VariantFilter` 类更名为 `JointCaller`，输出目录 `04_filter` → `04_joint`，最终产物 `all_samples.pass.vcf.gz`；断点续传保留
- 预测模块（`deeploc`/`deeptmhmm`/`tmhmm`/`signalp`）日志统一改**三 handler 分离**（stdout INFO→.out / stderr WARNING+→.err / file DEBUG+，§2.3），日志文件集中到 `99_logs/`；CLI 可选参数改为**始终显式透传**（避免改 config 默认值时 CLI 层静默失效）

### Added
- `deeploc`/`deeptmhmm`/`tmhmm`/`signalp`：运行成功后生成 `00_pipeline_info/software_versions.yml`（§12.5，best-effort 检测工具版本，失败仅告警不阻断）
- `deeploc`/`signalp`：主输出已存在时**断点续传跳过**昂贵预测步骤（§10.2）
- `signalp`：`signalp_path` 支持 `SIGNALP6_PATH` 环境变量 / config.yml 覆盖（§11.B）

### Fixed
- `deeploc`/`signalp`：conda 包装从 `os.path.basename(工具路径)` 改为传**完整路径**（§13.6.1），否则 `get_conda_env` 丢 `/envs/` 无法注入 `-n <env>`、在错误环境运行；signalp 同时从 `shell=True` 字符串改为 `shell=False` 列表（防注入）
- `deeploc`/`deeptmhmm`/`signalp`/`tmhmm`：`__post_init__` **先 `expand_path` 再 `mkdir`**（§11.3.1），原顺序会对字面 `~` 输出目录建出名为 `~` 的目录、`output_path` 与 `output_dir` 分家
- `deeptmhmm`：搬运结果后误查源目录致**每次误报全部输出缺失**；改为与实际搬运目标文件名比较
- `tmhmm`：输出改**原子写**（`.tmp` → `os.replace`，失败清理），避免断点续传解析截断文件；以 `output_dir` 为 cwd 防 plot 文件散落
- `signalp`：结果统计 `sp_count/total_count` 加 `total_count>0` 防护，空输入不再 ZeroDivision
- `blast`：比对可视化**反链坐标修复**——`sstart>send`（反链命中）时 subject 坐标应递减，原代码恒递增致反链比对坐标错误（html/text 两处对齐生成器同步修复）
- `blast`：区分**零命中**（rc=0 输出为空=合法结果，优雅保留）与真正失败（无输出文件），原代码把零命中当失败
- `blast`：`samples_count` 始终为 0——改为取 `len(sample_stats)`；无法自动检测类型时默认 `blastn`+告警（原留 None 致下游崩溃）
- `blast` CLI：`--auto-detect-samples` 从 `is_flag+default=True` 改为 `--x/--no-x` 开关（原写法使 `--no-auto-detect-samples` 不可达）
- `core`：`BaseAnalyzer` 透传 `config.log_file` 给 `setup_logger`，否则 `--log-file` 是死参数

### Docs
- CLAUDE.md v2.16 → **v2.17**：瘦身拆分（76.6KB→24.5KB），核心规则保留，完整模板/paths.py 实现/命名示例/conda 故障排查下沉到 `docs/dev-standards/` 五个按需参考文档（该次文档提交未单独 bump 包版本，在此补记）

## [1.21.3] - 2026-07-24

### Docs
- CLAUDE.md v2.15 → **v2.16**：§12.2.1 目录结构默认改为 **by-step**（多样本共享带序号步骤目录 + 文件名前缀 `{sample}_xxx` 区分，扁平便于统一查看某步骤所有样本结果），by-sample 降为可选（单样本产出大量独立文件、需按样本隔离时）；§12.6.1 加注释标注 genomescope 既有 by-sample 结构为新模块规范前的示例

## [1.21.2] - 2026-07-24

### Changed
- `rnaseq2vcf`：输出目录结构从 **by-sample**（`output_dir/{sample}/01_qc/...`）改为 **by-step**（`output_dir/01_qc/{sample}_1.clean.fq.gz`，4 个步骤目录 `01_qc`/`02_align`/`03_calling`/`04_filter` 共享、文件名带 `{sample}` 前缀区分多样本）。扁平结构便于统一查看某步骤的所有样本结果；断点续传不受影响（checkpoint 仅用于共享的 `genome_index`）。注：与 CLAUDE.md §12.2.1「多样本建议 by-sample」相反，按实际使用偏好选择 by-step

## [1.21.1] - 2026-07-24

### Changed
- `rnaseq2vcf`：样本发现改为**自动识别**（对齐 `fastq2vcf_gtx`）——`discover_samples` 的 R1 按优先级取首命中后缀（`_1.clean.fq.gz` > `_1.fq.gz` > `_1.fastq.gz`）、R2 多后缀尝试、`realpath` 解析软链；`read1_pattern`/`read2_pattern` 默认改 `None`（自动），CLI 新增可选 `--read1-pattern`/`--read2-pattern`。免去用户手动指定后缀，并杜绝 S1 误匹配 S10

### Chore
- `.gitignore`：新增 `/blast_results*.xlsx`（根目录 blast 测试样例；正式输出在 `output_dir/blast_dir/blast_results.xlsx`。该文件反复被 `copybiopytools` 从超算同步带回、干扰提交，ignore 后不再出现在 git status）

## [1.21.0] - 2026-07-24

### Added
- **新模块 `rnaseq2vcf`**：RNA-seq 变异检测（到 VCF）。HISAT2 比对（支持已知剪接位点或 de novo 发现 junction）→ GATK 变异检测 → VCF 过滤（FS/QD/聚类过滤）。HISAT2+samtools 管道用单个 `conda run -n ENV --no-capture-output bash -c '{pipeline}'` 包裹（§13.2.1，避免 conda run|conda run）；断点续传；by-sample 输出（`{sample}/01_qc/.../04_filter/{sample}.pass.vcf.gz`）；tmp 用 `output_dir/tmp`（§12.4.1）；命令执行日志记完整命令到 INFO（§2.2.1）；支持 `--dry-run`/`--step`/`--force`/`--skip-qc`/`--clean-fastq-dir`
- CLI：`main.py` 注册 `rnaseq2vcf` 命令

### Fixed
- `braker`：从头运行时清空 `braker_safe_dir`（flock 已保证独占），避免上次 GeneMark-ETP 失败残留的 `arx/proteins.fa` 被复用、因重复 ID 再次崩溃（实发）
- `braker`：`find_protein_files_in_directory` 排除自身中间文件（`cleaned_*`/`*_merged_proteins.fa`/`*_fixed_proteins.fa`）+ 按真实路径去重，避免合并结果滚雪球累积、ID 碰撞致 GeneMark-ETP "duplicated protein ID" 崩溃（实发）

## [1.20.3] - 2026-07-24

### Fixed
- `blast`：Excel 导出的 `raw_results` sheet 修复——BLAST outfmt 6 **无表头**，原 `pd.read_csv` 会把第一行数据当作列名；改为 `header=None, names=<15 标准字段>`（与 outfmt 配置 `main.py` line 353 精确匹配）+ 从文件名提取 sample 新增 `Sample` 列以区分合并的多样本

### Changed
- `blast`：Excel 导出美化——sheet 标签色（raw_results/summary/sorted/high_quality 四色区分）、冻结首行、列宽自适应（采样前 100 行估算）、表头深蓝底白字加粗居中、全表细边框

## [1.20.2] - 2026-07-24

### Fixed
- `braker`：全新训练（非 `--useexisting`）前清理上次遗留的 augustus `species/` 目录——否则 `braker.pl` 在 line 3490 因 `"species/<name> already exists"` **每次重跑必崩**。此前仅 `--useexisting` 分支（pipeline.py line 691）清理了 species 目录，而两个全新训练分支（`braker.gtf` 残留清理后重跑 / 完全从头运行）未清理。新增逻辑对称覆盖这两个分支；`use_singularity` 条件用于保护 `augustus_config_dir` 的定义作用域（该变量仅在 singularity 分支内定义）

## [1.20.1] - 2026-07-24

### Fixed
- `braker`：`braker_safe_dir` 的 md5 哈希改用**绝对路径**——原实现用相对路径参与哈希，而中文路径场景下 `smart_normalize_path` 会主动保持相对路径（`config.py` line 55-58），导致不同项目用相同的相对 `-o`/`-g` 参数时**哈希碰撞**、共用同一个 `~/tmp/braker_run_xxx` 工作目录，跨项目结果互相污染（v4/v5 即因此被污染）。修复后 `os.path.abspath()` 保证不同项目（不同绝对路径）得到唯一哈希。⚠️ 副作用：含中文路径的已有项目其 `braker_safe_dir` 目录名会改变，断点续传（`--useexisting`）会视为新项目重新运行

## [1.20.0] - 2026-07-24

### Added
- `blast`：结果导出多 sheet Excel（`blast_results.xlsx`，含 raw_results/summary/sorted/high_quality）；pandas/openpyxl 为**软依赖**（缺失则 warning + 跳过，TSV 照常输出，优雅降级）；单 sheet 超 Excel 行限（1048576）自动跳过该 sheet
- `braker`：项目级 `flock` 运行锁——防止同一项目并发运行在共享的 `braker_safe_dir` 内互踩（braker.pl 的 `ln -s` 等非原子操作会触发 "File Exists" 致命错误，白跑数小时后才崩）；非阻塞（`LOCK_NB`，锁被占用立即失败）；fd 绑定、进程退出（含 kill）自动释放，无死锁/残留锁；锁文件置于 safe_dir 同级（safe_dir 可能被清理重建，锁文件不随之删除）

## [1.19.0] - 2026-07-23

### Added
- `func_anno`：多样本支持——`-i` 传目录自动识别为多样本（by-sample，每文件一子目录），传单文件为单样本（by-step 不嵌套）；新增 `--by-sample` 强制单文件也嵌套（往同一 `-o` 多次跑不覆盖）；多样本时单样本失败不阻断其他
- `func_anno`：KEGG 通路过滤——新增 `--kegg-exclude-keywords`（name 黑名单，`None`=内置植物无关词 cancer/estrogen 等）与 `--kegg-exclude-categories`（按 BRITE 大类排除，如 Human Diseases）

### Changed
- **临时目录统一改造（CLAUDE.md v2.15）**：35 个模块临时文件/目录从系统 `/tmp` 改为 `output_dir/tmp` 子目录，消除超算 `/tmp` 爆满风险。两种向后兼容 Pattern：直接 `dir=output/tmp`（fastq_stats/centier/insert_detection 等）；函数新增可选 `tmp_dir`/`base_dir` 参数、不传则回退系统 `/tmp`（phyto_effector `run_signalp3`/`split_fasta_chunks`、haphic `create_temp_dir` 等）。phyto_effector `run_signalp3` 顺带加 `try/finally` 异常清理
- `func_anno`：建表阶段去除断点续传（纯解析秒级，每次重建确保过滤参数生效）
- CLAUDE.md：v2.14 → v2.15，更新 §12.4.1/§12.6/§12.7（临时目录改 `output_dir/tmp`）

### Fixed
- `func_anno`：eggnog annotations 路径修正为 `eggnog_dir/01_emapper/{sample}.emapper.annotations`（匹配 eggnog_mapper 实际输出结构，原路径找不到文件导致断点续传误判）

## [1.18.0] - 2026-07-22

### Added
- **新模块 `trimal`**：多序列比对自动修剪（trimAl）封装。支持 6 种修剪方法（automated1/gappyout/strict/strictplus/gt/cons）、8 种输出格式、`-colnumbering` 新旧列号映射、`-complementary` 互补比对、`-backtrans` AA→NT 反向翻译。trimal 作为独立编译 C++ 二进制（RPATH 解析动态库），直接以绝对路径调用、不经 `conda run` 包装（实测 0.01s）。断点续传；by-step 目录（`01_trimal`/`00_pipeline_info`/`99_logs`）+ `software_versions.yml`
- **新模块 `phylo_trim`**：整合 mafft-fasttree + trimal，自动产出 trimal 前后两棵 FastTree 系统发育树（序列→MAFFT→FastTree→trimal→FastTree）。仅 import 复用 `mafft_fasttree`（`PhyloTreeBuilder`/`FastTreeBuilder`/`SequenceProcessor`/`CommandRunner`）与 `trimal`（`TrimalRunner`），不修改两模块源码；序列类型自动检测；`--skip-trimal` 仅出前树；三段式断点续传（前流程/trimal/后树）；by-step 目录（`01_mafft_fasttree`/`02_trimal`/`03_fasttree_trimmed`）
- CLI：`main.py` 注册 `trimal` 与 `phylo-trim` 两个命令

## [1.17.1] - 2026-07-22

### Fixed
- `get_link_from_CNCB`（CLI 包装层）：移除遗留的 `--threads` Click 选项 / `threads` 形参 / `max_threads=threads` 传参——1.16.0 已从 `CNCBConfig` 删除 `max_threads` 死字段，但 Click 包装层漏删，导致 `CNCLinkExtractor(**kwargs)` 经 `**kwargs` 转发给 `CNCBConfig` 时抛 `TypeError: unexpected keyword argument 'max_threads'`，`biopytools get-link-from-CNCB` 启动即崩溃

### Changed
- `get_link_from_CNCB`（docs）：清理 `docs/get_link_from_CNCB.md` 中 6 处线程相关残留——参数表 `--threads` 行 1 处、示例命令 `--threads` 5 处（其中 4 处 `python -m ...main --threads` 本就与 argparse 实现不符）、Python 配置 dict 中 `'max_threads': 6` 1 处（照抄会经 `CNCLinkExtractor(**config)` 直接触发上述 TypeError）

## [1.17.0] - 2026-07-21

### Added
- **新模块**：`func_anno` — 蛋白功能注释流水线（interproscan 结构域 + eggnog-mapper GO/KEGG → 标准 GO/KEGG 表衔接下游 R）；braker4ps 模式（不改 interproscan/eggnog_mapper 源码，仅 import 调用）；支持 `--ips-result`/`--eggnog-result` 复用已跑结果、`--skip-ips`/`--skip-eggnog`、断点续传、外部 KEGG 映射补 category；by-sample 目录（01_interproscan/02_eggnog/03_tables/99_logs）

### Changed
- `blast`：`BLASTAnalyzer` 改继承 `core.BaseAnalyzer`（复用通用 run_command/依赖检查），删除 utils.py 中 8 个冗余函数（check_tool_availability/get_input_files/run_command/check_blast_dependencies 等，均 0 调用方）；main.py 重构精简
- `hite`：新增 `hite_runner.py`/`main.py` 作为单基因组入口（`HiteRunner` 用 `singularity_path`+`build_singularity_command` 直接挂载，支持断点续传），`HiteConfig` 字段重构（singularity_path）；panHiTE 流程保留（PanHiteConfig + SingularityContainerManager）

## [1.16.0] - 2026-07-21

### Added
- `sra2fastq`：断点续传（已转换的 SRA 跳过）+ `.fastq.gz`→`.fq.gz` 统一重命名 + skipped 统计
- `vcf2genotype`：`--dry-run` 试运行模式（只扫描统计不写数据）+ 变异/染色体分布统计
- `poplddecay`：`-o` 指向目录时自动从输入文件名推导前缀（避免产出 `.stat.gz`/`.png` 等隐藏文件）
- conda 环境：新增 `foldseek`、`mga`

### Changed
- 多模块路径标准化加 `expand_path`（dual_rnaseq/fastp/get_link_from_CNCB/smudgescope/vcf2genotype/rnaseq）
- `rnaseq`/`dual_rnaseq`：CommandRunner 用 `conda run -n ENV bash -c` 包装整条含管道的命令（§13.2.1，避免 conda run|conda run 双重包装）；`fastp` 单工具自动 conda 包装
- `rnaseq`：删除 alignment/config/quantification/utils 文件顶部大段注释死代码；CommandRunner 新增 dry_run
- `rmvp`：`--r-env` 默认改环境名 `rMVP`（与 config 对齐）
- `smudgescope`：CLI `_lazy_import_genome_analysis_main`→`_lazy_import_smudgescope_main`（模块名修正）；software_versions 用 build_conda_command 取版本

### Fixed
- `smudgescope`：`fmt_range` 定义顺序导致 `UnboundLocalError`（het_point 为 None 时引用未定义函数）；openpyxl 新版 `.copy()` 弃用→`Font` 相加
- `fastq2vcf_gtx`：GTX 索引清单错误（原硬编码 .gtx/.gtx.bwt/.gtx.sa，实际产物是 .amb/.ann/.pac + .bwt.* 后缀随基因组大小变化）→改 glob 匹配

### Removed
- `vcf2genotype`：excel 输出（流式难以正确生成 xlsx，且基因型表常超 Excel 行数上限），保留 txt/csv
- `get_link_from_CNCB`：`max_threads`/`--threads`（CNCB 为生成下载链接，无并发需求）

## [1.15.2] - 2026-07-20

### Fixed
- `annovar`：`GFF3Validator.clean_and_fix_gff3` 不再原地覆盖用户输入的 GFF3（原实现备份+替换输入文件，危险副作用）——改为读输入、写 output_dir 内的工作副本，输入保持不动；`gff3_to_genepred` 配套使用工作副本（`{build_ver}.cleaned.gff3`）并重定向后续步骤
- `annovar`（CLI）：`--step`（int）传 sys.argv 未转 str → `str(step)`

### Changed
- `annovar`：删除 `database_path` 参数（改用 output_dir 作为 ANNOVAR 数据库目录，refGene 在此生成）；工具路径 gffread/seqkit 改 `get_tool_path`（§11.3）+ `build_conda_command` 传完整路径（§13.6.1）；`data_processing.py` 缩进 Tab→4 空格统一

## [1.15.1] - 2026-07-17

### Fixed
- `longestmrna`：过滤最长转录本无 CDS 的非编码基因——保留会导致 gene_info 与蛋白输出条目数不一致；现跳过并统计（WARNING + summary 的 `noncoding_skipped`）

### Changed
- `longestmrna`：`GFF3Parser` 合并进 `CDSCalculator.calculate_from_gff`（DRY，删除冗余类）；工具路径 gffread/seqkit 改 `get_tool_path`（§11.3，支持 `GFFREAD_PATH`/`SEQKIT_PATH`）+ `build_conda_command` 传完整路径（§13.6.1）；`__init__.py` 文档导入路径 `longest_mrna`→`longestmrna`；`TempFileManager.cleanup` 末尾清空列表避免复用重复删除

## [1.15.0] - 2026-07-17

### Added
- `gff_renamer`：新增 `--prefer-mrna`（默认开）——对含 mRNA 的基因丢弃冗余 transcript（misc_RNA 变体）及子特征，针对 NCBI/Gnomon EGAPx 注释；`--no-prefer-mrna` 关闭
- `fastq2vcf_gtx`：变异过滤新增 `--snp-qual`/`--indel-qual` 质量过滤；`CLUSTER_MODE` 哨兵（集群手动投递模式不计为失败）；faketime 依赖检查（GTX 命令依赖它绕过 license 时间校验）

### Changed
- `fastq2vcf_gtx`：`GenomeIndexer.build_genome_index` 合并原 main.py 的 `_force_build_gtx_index` 为单一入口（断点续传 + 索引已存在跳过 + 统一 GTX 索引清单）；断点续传默认启用；`--input` 改可选（与 `--clean-fastq-dir` 二选一）；`_checkpoint_manager` 工厂函数统一 CheckpointManager 构造（DRY）
- `gff_renamer`：UTR 始终收集并重编号（保证输出 ID 唯一），`include_utr` 现仅影响 UTR 重排序/Name 清理

### Fixed
- `fastq2vcf_gtx`：`run_with_progress` 管道死锁修复——改独立 reader 线程持续排空 stdout（原 readline 阻塞 + communicate 反模式可能死锁）；超时优雅降级（terminate→wait→kill）
- `gff_renamer`：`id_mapping` 的 exon/intron/CDS/UTR/codon key 改 `(line_index, old_id)` 元组，修复 AGAT/合并共享 old_id 导致输出重复 ID（查询方双重 try 兼容 str/tuple key）
- `gtx/utils`：裸 `except:` 改 `except Exception:`

## [1.14.2] - 2026-07-17

### Changed
- `hifi_hic`：命令执行统一走 conda 包装（§13）——单命令 `run_command`(list+shell=False)、shell 命令 `run_shell_command`(shell+LD_LIBRARY_PATH，含 tee/awk/重定向)、管道 `run_pipeline_command`(剥 conda run+LD_LIBRARY_PATH)；`assembler`/`ngs_polisher` 的 hifiasm|tee、awk|seqkit、samtools merge|sort 跨 env 管道改完整路径+LD_LIBRARY_PATH（§13.2.1）；工具路径(hifiasm/seqkit/samtools)改 `get_tool_path`+`__post_init__` expand（§11.3）；日志改命名 logger + stdout/stderr/file 三 handler 分离（§2.3）；目录 `05_statistics`→`00_pipeline_info`、`06_logs`→`99_logs`（§12）；新增 `software_versions.yml`（§12.5）；CLI threads 默认 12→88，移除冗余 `--resume`/`--verbose`
- `hifi_hic`：`_get_fasta_summary`/read names 提取/`wc -l` 改 Python 原生（移除 `os.popen`/grep|cut/wc 管道）；`purge_dups_wrapper` 抽出 `_find_purged_output`(DRY)+异常改 `exc_info=True`

### Fixed
- `hifi_hic`：`assembler`/`purge_dups_wrapper` 新增断点续传（fasta/去冗余输出已存在则跳过，§10.2）；错误消息前导空格清理

## [1.14.1] - 2026-07-16

### Changed
- `haphic`：工具路径（haphic/bwa/samtools/matlock/samblaster/filter_bam/3d-dna 系列）全面改用 `get_tool_path`（env>config>默认）+ `__post_init__` `_expand_tool_path` 展开（§11.3）；日志与报告落 `99_logs/`；Juicebox 生成流程重构拆分为 `_resolve_filtered_bam`/`_ensure_mnd_file`/`_sort_mnd_file`/`_generate_hic_file` 小方法；BWA/filter_bam 管道补 `LD_LIBRARY_PATH`（§13.2.1 方案A）；新增 `software_versions.yml`（§12.5）+ `00_pipeline_info`/`99_logs` 目录（§12.2.3）；调试日志 info→debug
- `haphic`（CLI）：samblaster 默认 env 改 `haphic`、matlock 默认改完整路径 `~/miniforge3/envs/juicer_v.1.6/bin/matlock`

### Fixed
- `haphic`：`get_step_status` 用错 key `scaffolds_fa`→`scaffolds_fasta`（原 KeyError）；删除 `main.py` 冗余 `main()`（CLI 经 `HapHiCProcessor` 类入口，不依赖）
- `chr_rename`：CLI 示例命令名 `chr_rename`→`chr-rename`（与注册命令一致）
- `deeptmhmm`：`predict.py` 加 `cwd=deeptmhmm_dir`（相对路径依赖）

## [1.14.0] - 2026-07-16

### Added
- `kmertools`：所有子命令（build/query/split-fasta/gen-fof/import-db/extract/count）成功后生成 `00_pipeline_info/software_versions.yml`（§12.5，含工具版本探测 + 关键参数 + 运行时间）；新增 `ensure_pipeline_dirs`、`generate_software_versions_yml` 辅助函数

### Changed
- `kmertools`：工具路径（kmtricks/kmindex/bgzip/unikmer/jellyfish）全面改用 `get_tool_path` 按优先级解析（环境变量 > 配置文件 > 默认）+ `__post_init__` `expand_path` 展开 `~`（§11.3）；所有用户路径（input/output/fof/header/tmp/rocksdb/bed/fasta 等）统一 `expand_path`
- `kmertools`：`run_command`（utils 模块级 + `CommandRunner.run`）命令日志由 DEBUG 改 INFO（§2.2.1），并改用 `shell=False`（字符串按空白拆分）；调用处删除重复的 `运行命令` 日志（命令记录已统一到 `run_command` 内部）
- `kmertools`：`file_processor` 的 `zcat|head|wc` 管道改 Popen 链（无 shell，§13.2.1）；`extract/main.py` unikmer 调用 `shell=True` → 列表 `shell=False`
- `kmertools`：`KmerToolsLogger` 改命名 logger + stdout(INFO)/stderr(WARNING+)/file(DEBUG+) 三 handler 分离（§2.3）；各子命令日志集中到 `99_logs/`
- `kmertools`：`query` 矩阵转置由 awk+sed 改纯 Python（`zip(*rows)`），修复原 `awk脚本.split()` 喂 `shell=False` 导致转置失效的 bug

### Fixed
- `kmertools/kmer_count/main.py`：修复 `f"...总样本数...: \1` 未闭合字符串导致的 SyntaxError（kmer_count 模块完全无法导入/运行）
- `kmertools/kmer_count/jellyfish_processor`：删除误引入的 `from modelscope.hub.errors import FileIntegrityError`（无关依赖，未安装时阻断导入）
- `kmertools/kmer2vcf/utils`：`zstandard` 改延迟导入（仅处理 .zst 时需要），缺失时给清晰安装提示，避免硬依赖阻断模块导入

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
