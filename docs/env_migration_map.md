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

## 三、可删除旧环境清单|Deletable Old Envs（51 个）

> ⚠️ **删除前提**：先把新代码同步到超算（GitHub 拉取），确认模块跑通后再删。
> 删除脚本：scripts/delete_obsolete_envs.sh（默认 dry-run）

BUSCO_v.6.0.0, BioinfTools, GATK_v.4.6.2.0, Genome_dedup, K-mer, RNA_Seq,
RSeQC_v.5.0.4, agat_v.1.7.0, bcftools_v.1.22, canu_v.2.3, eggnog-mapper_v.2.1.15,
eviann_v.2.0.5, freebayes, genomescope2_v.2.1.0, genometools_v.1.6.5,
getorganelle_v.1.7.71, gffcompare_v.0.12.10, haphic, hifiasm_v.0.25.0,
iqtree_v.3.0.1, kakscalculator2_v.2.0.1, kmindex_v.0.6.0, kmtricks_v.1.5.1,
ltr_retriever_v.3.0.1, mafft_v.7.525, meme_v.5.5.9, miniprot_v.0.18,
needle_v.1.0.3, newick_utils_v.1.6, orthofinder_v.3.1.5, pairtools_v.1.1.3,
pan-blocks, phobius_v.1.0.1, pixy_v.2.0.0, poplddecay_v.3.43, primer3_v.2.6.1,
purge_dups_v.1.2.6, pycirclize_v.1.10.1, repeat_identiy, repeatmodeler_v.2.0.7,
resistify_v.1.3.0, samplot_v.1.3.0, selective_sweep, signalp6, spades_v.4.3.0,
sv_calling, tidk_v.0.2.65, tmmhmm_v.2.0c, trimal_v.1.5.0, wgdi_v.0.75,
yahs_v.1.2.2

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
