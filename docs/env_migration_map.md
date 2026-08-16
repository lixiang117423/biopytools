# 新旧 Conda 环境对应关系|Old → New Conda Env Mapping

> 用途|Purpose: debug 时快速查「旧环境去哪儿了」；删除旧环境前核对。
> 数据源|Source: docs/conda_env_audit.csv + scripts/migrate_env_paths.py + 超算实测核验
> 更新|Updated: 2026-08-16

## 一、14 个功能域环境|The 14 Domain Envs

| 域|Domain | 吸收的旧环境|Absorbed Old Envs |
|---|---|---|
| align | GATK_v.4.6.2.0, bcftools_v.1.22, sv_calling, freebayes, Genome_dedup |
| pop | Population_genetics, selective_sweep, adamixture_v.1.0.2, treemix_v.1.13, pixy_v.2.0.0, poplddecay_v.3.43 |
| asm | canu_v.2.3, hifiasm_v.0.25.0, kmc_v.3.2.4, K-mer, merqury_v.1.3, purge_dups_v.1.2.6, genomescope_v.2.0.1, genomescope2_v.2.1.0, tidk_v.0.2.65, getorganelle_v.1.7.71, spades_v.4.3.0 |
| hic | haphic, pairtools_v.1.1.3, yahs_v.1.2.2, juicer_v.1.6 |
| annot | Augustus_v.3.5.0, agat_v.1.7.0, gffcompare_v.0.12.10, miniprot_v.0.18, transdecoder_v.5.5.0, eggnog-mapper_v.2.1.15, orthofinder_v.3.1.5, Blast_v.2.16.0, genometools_v.1.6.5, eviann_v.2.0.5 |
| repeat | repeatmodeler_v.2.0.7, repeat_identiy, ltr_retriever_v.3.0.1, ltr_harvest_parallel_v.1.2, ltr_finder_parallel_v.1.3 |
| rna | RNA_Seq, RSeQC_v.5.0.4 |
| protein | signalp6, resistify_v.1.3.0, needle_v.1.0.3, meme_v.5.5.9, phobius_v.1.0.1, tmmhmm_v.2.0c |
| phylo | iqtree_v.3.0.1, mafft_v.7.525, trimal_v.1.5.0, newick_utils_v.1.6, wgdi_v.0.75, kakscalculator2_v.2.0.1 |
| pan | pggb_v.0.7.4, kmtricks_v.1.5.1, kmindex_v.0.6.0, pan-blocks, panman_v.0.1.4 |
| viz | samplot_v.1.3.0, pycirclize_v.1.10.1 |
| misc | iseq_v.1.9.8, primer3_v.2.6.1, bbmap_v.39.81, BioinfTools |
| r | rMVP, WGCNA_v.1.73 |
| busco | BUSCO_v.6.0.0 |

## 二、保留不动的旧环境|Kept Old Envs（勿删|DO NOT DELETE）

### Tier3 例外|Tier-3 Exceptions（巨无霸/锁生态/基础设施）
| 环境|Env | 原因|Reason |
|---|---|
| qiime_v.2024.10.1 | qiime2 生态锁死 |
| picrust_v.2.6.3 | qiime2 依赖 |
| EDTA_v.2.2.2 | 428 包巨无霸 |
| EGAPx_v.0.4.0-alpha | 完整应用环境 |
| cphasing | 394 包巨无霸 |
| jcvi_v.1.5.7 (JCVI_v.1.5.6) | 354 包巨无霸 |
| telocomp | 含 flye, 观察期 |
| rnaseq_val | 291 包大环境 |
| singularity_v.3.8.7 | 基础设施 |
| base | 基础设施 |
| psvcp_v.1.0.1 | 179 包混合栈(R+python) |

### Tier2 legacy|Tier-2 Legacy（py3.7/3.8 老栈，升级验证后并入）
| 环境|Env | 升级去向|Upgrade Into |
|---|---|
| braker_v.3.0.8 | annot（需 py3.7→3.13 重装验证） |
| fanc_v.0.9.23b | hic |
| deeptmhmm_v.1.0 | protein |
| DeepBSA | r |
| HiC-Pro_v3.1.0 | 与 SubPhaser 合成 hic-legacy |
| SubPhaser | 同上 |
| pasa_v.2.5.3 | 保持独立（python2 依赖） |

### 源码/手工安装（不在 conda）|Source / Manual Installs
| 环境|Env | 说明|Note |
|---|---|
| gctb | 自编译二进制 |
| plothic_v.1.0.0 | 源码安装 |
| genomesyn2 | 源码安装（perl 脚本） |
| alignoth | 源码安装 |
| a-liner | 固定 conda 环境（aliner_env） |
| kmeriaenv | kmeria 模块专用 |
| vcf2gwas_v.0.8.9 | 独立工具 |

### 黑名单|Blacklist（同名二进制含义不同，禁止合并）
| 环境|Env | 原因|Reason |
|---|---|
| signalp_v.3.0b | bin/signalp 是 SignalP 3（Perl），与 protein 域的 SignalP 6 别名完全不同 |

### 待修复|Pending Fix（修复后并入域环境）
| 环境|Env | 问题|Issue |
|---|---|
| adamixture_v.1.0.2 | 代码 4 处引用 bin/adamixture，疑似笔误（应为 admixture），需超算核实 |
| vg_v.1.7.0 | vg 模块 env 名字段，pan 域已有 vg（版本 1.63.1 vs 旧 1.70.0 待功能验证） |
| kmc_v.3.2.4 / merqury_v.1.3 | 仍被裸环境目录引用（conda_env 字段），验证后可切 |
| xxx | nlr_annotator 的 Java 注释占位符 |

## 三、可删除旧环境清单|Deletable Old Envs（155 个）

> 标准|Rule: **模块（代码）未引用的环境一律可删**（含手工工具孤儿环境）。
> ⚠️ **删除前提**：先把新代码同步到超算（GitHub 拉取），确认模块跑通后再删。
> 删除脚本：scripts/delete_obsolete_envs.sh（默认 dry-run，清单文件 scripts/delete_list.txt）

### 分组统计|Groups

| 组|Group | 数量|Count | 说明|Note |
|---|---|---|---|
| 已迁移到域环境 | 42 | 代码引用已全部切换（51 个中 9 个仍被裸目录/字段引用，转入保留） |
| 版本分裂对（旧版） | 14 | DeepTMHMM / EDTA / edta_v.2.3.0 / Orthofinder_v.3.0.1b1 / R_v.4.5.1 / clustalo / gfatools / jellyfish_v.1.1.3 / muscle_v.5.3 / quast_v.5.3.0 / rnabloom_v.2.0.1 / tesorter_v.1.4.7 / trinity_v.2.15.2 / vg_v.1.67.0 |
| 仅文档提及 | 4 | biopytools / centier / pbbam_v.2.4.0 / smudgeplot |
| 纯孤儿（手工工具） | 94 | 下载工具/质检/转录组/宏基因组/GWAS 等，从未被模块引用 |
| 怪异名 | 1 | Name（疑似误导出产生的空壳环境） |

### 删除前人工确认项|Double-Check Before Deleting

- **同事共享风险**：这些环境在共享路径 /share/org/YZWL/yzwl_lixg，若同事手工使用过某些环境（尤其 aria2/axel/aspera 下载类、checkm/quast 质检类），请先沟通
- **PATH 查找类模块**：rnabloom 模块默认从 PATH 找 rnabloom 二进制（不硬编码环境），删除 rnabloom_v.2.0.1 后需用 --rnabloom-path 指定或重新装进某域
- **centier**：centier 模块用 PATH 搜索 + conda 兜底，删除后依赖检查会提示重建

### 全量 155 个（scripts/delete_list.txt 为准）

3d-dna_v.201008, Baidu, BioinfTools, CentroMiner, DeepTMHMM, EDTA, Effector_annotation,
GAPIT, GTDB-Tk_v.2.4.1, InterProScan.v.5.75-106.0, JAVA_11, K-mer,
LinearPangenomeBuilder, Mumemto_v.1.3.0, NLR_Annotation_Pipeline, Name,
Orthofinder_v.3.0.1b1, Python_v.3.13.5, RNA_Seq, RSeQC_v.5.0.4, R_v.4.5.1,
RagTag_v2.10., Syri_v.1.7.1, Tiberius_v.1.1.1, agat_v.1.7.0, allhic_v.0.9.14,
annevo_v.2.1, aria2, aspera_v.3.9.6, assembly_stats_v.1.0.1, axel_v.2.17.13,
bamtools, bbtools_v.37.62, bcftools_v.1.22, biomformat, biopytools, bismark_v.0.24.2,
braker4_v.0.5.0, cafe_v.5.1.0, canu_v.2.3, caster_v.1.23, centier, checkm_v.1.1.0,
clustalo_v.1.2.4, cooler_v.0.10.2, csvkit, deeploc20, deeploc_v.2.1, diamond_v.2.1.13,
edta_v.2.3.0, effectR, egepx_v.0.5.1, eggnog-mapper_v.2.1.15, emboss,
entrez-direct_v.24.0, evidencemodeler, faketime, fastani_v.1.34, fastme_v.2.1.6.3,
fastqc_v.0.12.1, fcs-gx, fithic-v.2.0.8, foldseek, freebayes, gapcloser_v.1.2.1,
gemma_v.0.98.5, gemoma_v.1.9, gene-anno, gene_family, genehapr_v.1.2.4,
genomethreader_v.1.7.1, genometools_v.1.6.5, getorganelle_v.1.7.71, gfatools_v.0.5.5,
gffcompare_v.0.12.10, glnexus_v.1.4.1, haphic, hicexplorer, hicexplorer_v.3.7.6,
hifiasm_v.0.25.0, hmmer_v.3.4, iqtree_v.3.0.1, jellyfish_v.1.1.3, julia_v.1.12.2,
julia_v.1.7.2, kakscalculator2_v.2.0.1, kat_v.2.4.2, kmeria_v.2.0.1, kmindex_v.0.6.0,
kmtricks_v.1.5.1, kofamscan_v.1.3.0, kraken_v.2.17, mafft_v, mafft_v.7.525,
mcscanx_v.1.0.0, meme_v.5.5.9, metaWRAP_v.1.2, metagraph, methylkit, miniprot_v.0.18,
mmseqs2_v.16.747c6, mrna_prediction, mummer_v.3.23, mummer_v.4.0.1, muscle_v.5.3,
mutmap, ncbi-datasets-cli, needle_v.1.0.3, newick_utils_v.1.6, orthofinder_v.3.1.5,
pairtools_v.1.1.3, paml_v.4.10.9, pan-blocks, pbbam_v.2.4.0, pbsv_v.2.11.0,
phobius_v.1.0.1, pilon_v.1.24, plotsr_v.1.1.1, poplddecay_v.3.43, primer3_v.2.6.1,
puzzle-hi-c, qtlseq, quast_v.5.3.0, racon_v.1.5.0, raxml_v.8.2.13, repeat_identiy,
repeatmodeler_v.2.0.7, resistify_v.1.3.0, rnabloom_v.2.0.1, ruby_v.3.4.7,
salmon_v.1.10.3, samplot_v.1.3.0, selective_sweep, seqtk, signalp6, smcpp_v.1.15.2,
smudgeplot, sniffles_v.2.6.3, spades_v.4.3.0, sratoolkit_v.2.5.7, starship_v.1.26.0,
sv_calling, swave_v.1.2, tassel, taxonkit_v.0.20.0, tesorter_v.1.4.7, tidk_v.0.2.65,
tmmhmm_v.2.0c, trimal_v.1.5.0, trinity_v.2.15.2, vg_v.1.67.0, wgdi_v.0.74,
wgdi_v.0.75, yahs_v.1.2.2, zsh

### 保留 59 个（scripts/protect_list.txt 为准）

14 个域环境 + base + 44 个仍被模块引用的旧环境（Tier2/Tier3/黑名单/待修复）。


## 四、死引用修复记录|Dead-Ref Fixes（已修）

| 代码引用|Cited | 修复|Fixed |
|---|---|
| kakscalculator | → kakscalculator2_v.2.0.1（后并入 phylo） |
| flye (envs/flye/) | → telocomp 环境（config_template） |
| vg_v.1.7.0「写错」疑云 | 已核实环境存在（内含 vg 1.70.0），非 bug |
| raxml-ng_v.2.0.1 | 仅文档/测试引用，生产代码无引用 |

## 五、Debug 指引|How to Debug

1. **模块报「命令找不到」** → 查模块 config 的默认路径属于哪个域（本文件第一节）→ 超算确认 `ls ~/miniforge3/envs/<域>/bin/<工具>`
2. **工具在旧环境能用、新环境不能** → 用 `conda run -n <域> <工具> --version` 复现；版本差异见 `conda env list` 与 `mamba list -n <域>`
3. **怀疑误删** → 本文件第二节的「保留不动」清单是白名单；第三节之外的环境一律不删
