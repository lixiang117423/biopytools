# 软件调用与安装决策清单|Tool Invocation & Installation Policy

> 用途|Purpose: **AI 写代码时**的软件调用决策清单——从哪调、缺了装哪、何时新建环境。
> 配套|Companions: 代码侧映射 biopytools/common/env_map.py(TOOL_DOMAIN_MAP)、
> 调用公共层 biopytools/common/conda_runner.py、超算速查表 docs/conda_env_software_map.md。
> 更新|Updated: 2026-08-18 (1.33.0)

---

## 一、决策树|Decision Tree

AI 在新模块/新代码中需要调用外部软件时, 按以下顺序决定:

1. **优先查域环境映射表**(下方第二部分): 工具在表中 → 代码只写
   `get_domain_tool_path('<tool>', '<旧默认>', '<TOOL>_PATH')` 取路径 +
   `build_conda_command(path, args)` 包装执行。**代码里不写死环境名**。
2. **工具不在映射表**: 查超算速查表第二部分「保留独立环境」——
   在某个保留环境里(如 braker/cphasing/EDTA/qiime2/picrust2/EGAPx/telocomp/fanc 等) →
   用 `get_tool_path` 指向该环境路径, `build_conda_command` 自动包装(它从路径提取环境名)。
3. **非 conda 软件**: `~/software/` 第三方二进制、java -jar、singularity/apptainer 容器、
   wget/ascp、python/perl/R 解释器跑脚本 → **直接调用, 不包装**(公共层会自动跳过非 envs 绝对路径)。
4. **系统工具**(grep/awk/zcat/wc/df/which/cut/sort/head/tail/sed/tr/uniq/cat) → 裸名直调, 不属于 conda 范畴。
5. **管道命令**: 严禁 `conda run | conda run`(§13.2.2)——用公共层
   `build_pipeline_command` / `run_pipeline`(方案B: 管道段提取实际命令后直调域环境二进制)。

### 缺失怎么办|When a tool is missing

- 域环境缺工具(存在性检查失败): **提示用户安装到所属域环境**, 不要新建环境:
  ```bash
  conda install -n <域环境> -c bioconda -c conda-forge <工具>
  # 同步更新 envs/<域>.yml 配方, 保证环境可重建
  ```
- 域环境装不了(依赖冲突): 才允许使用保留独立环境; 若两者都不可行 → 见下条。
- **新建 conda 环境的唯一条件**: 工具与所有现有域环境依赖冲突、且不在保留独立环境清单中。
  新建前必须先: (a) 在 docs/conda_env_software_map.md 登记; (b) 在 envs/ 留 yml 配方;
  (c) 说明为何不能并入域环境。**禁止默默新建环境**。
- 模块运行时检测到工具缺失: 用 `check_tools` 报错并给出安装命令(如现有模块的
  `conda install -c bioconda bcftools` 提示), 不要静默失败。

---

## 二、工具 → 域环境清单|Tool → Domain Env List

代码调用方式统一为: `get_domain_tool_path(工具名, 旧默认, 环境变量) + build_conda_command`
缺失安装目标为所属域环境(`conda install -n <域> <工具>`)。

| 域环境<br>Domain | 工具<br>Tools |
|---|---|
| **align** | gatk, bcftools, bgzip, tabix, samtools, bwa, freebayes, bedtools, wgsim, minimap2 |
| **pop** | plink, vcftools, admixture, treemix, pixy, PopLDdecay, poplddecay, RAiSD, xpclr |
| **asm** | canu, hifiasm, kmc, kmc_tools, jellyfish, merqury, meryl, purge_dups, genomescope2, tidk, get_organelle_from_reads.py, spades.py, spades |
| **hic** | haphic, pairtools, yahs, juicer, matlock, samblaster, filter_bam |
| **annot** | augustus, bam2hints, etraining, gff2gbSmallDNA.pl, gtf2gff.pl, agat_convert_sp_gxf2gxf.pl, gffcompare, miniprot, TransDecoder.LongOrfs, TransDecoder.Predict, emapper.py, orthofinder, blastn, blastp, makeblastdb, diamond, mmseqs, gt, eviann |
| **repeat** | RepeatMasker, BuildDatabase, RepeatModeler, LTR_retriever, ltr_harvest_parallel, ltr_finder_parallel |
| **rna** | hisat2, hisat2-build, extract_exons.py, extract_splice_sites.py, stringtie, fastp, gffread, infer_experiment.py |
| **protein** | signalp6, signalp, resistify, hmmsearch, needle, needleall, embossversion, meme, phobius.pl, tmhmm, tmhmm2 |
| **phylo** | iqtree, iqtree2, mafft, trimal, nw_reroot, nw_display, wgdi, KaKs_Calculator, raxml-ng |
| **pan** | pggb, vg, kmtricks, kmindex, panman, nucmer, delta-filter, show-coords, minigraph |
| **viz** | samplot, pycirclize |
| **misc** | iseq, primer3_core, bbmap.sh, repair.sh, seqkit, sra-tools, fasterq-dump, prefetch, axel, pigz |
| **r** | mstmap |
| **busco** | busco(独立域环境) |

## 三、例外清单|Exception List

### 3.1 保留独立环境(域环境容纳不了, 继续用旧环境路径)
braker / fanc / deeptmhmm / EDTA / cphasing / jcvi / qiime2 / picrust2 / EGAPx /
telocomp / rnaseq_val / singularity / GATK / Genome_dedup / HiC-Pro / LTR_retriever /
merqury / pggb / vg / treemix / signalp 等(完整表见速查表第二部分)。
→ 代码用 `get_tool_path('<tool>', '~/miniforge3/envs/<旧环境>/bin/<tool>', '<TOOL>_PATH')`,
`build_conda_command` 自动从路径提取环境名并包装。

### 3.2 ~/software 第三方(直调, 不包装)
interproscan.sh / EGAPx / JanusX / GTX(外裹 faketime) / PanDepth / Dsuite /
FAPROTAX(collapse_table.py) / VCF2PCACluster / TeloComp / GenomeSyn / KGGSee(java -jar) /
gemma(~/.local/bin)。

### 3.3 容器与下载/解释器(直调)
cactus(singularity+sif) / parabricks(apptainer) / wget / ascp / java -jar /
python/perl/R 解释器跑脚本。

### 3.4 禁止使用
scripts/delete_list.txt 中的 154 个待退役环境(见速查表第四部分)。

---

## 四、验收检查|Review Checklist

- [ ] 新代码工具路径走 `get_domain_tool_path`/`get_tool_path`, 无硬编码绝对路径
- [ ] 外部命令执行前记录完整命令到 INFO(`命令|Command: ...`)
- [ ] 管道无 `conda run | conda run`, 走方案B
- [ ] 缺失工具给安装命令提示(指向域环境), 不新建环境(除非登记过)
- [ ] 若新建环境: 速查表已登记 + envs/*.yml 配方存在 + 理由写明
