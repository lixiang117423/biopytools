# Conda 环境软件速查表|Conda Environment Software Quick Reference

> 用途|Purpose: 超算上的 AI/开发者找现成软件时, 按本表选择调用环境。
> 更新|Updated: 2026-08-16

## 使用规则|Usage Rules

1. **优先**使用 14 个功能域环境(第一部分)|Prefer the 14 domain envs (Part 1)
2. 域环境没有的软件, 查保留独立环境(第二部分)|Fall back to kept standalone envs (Part 2)
3. **禁止**在 scripts/delete_list.txt(154 个待退役环境)中找软件, 也不要让新模块依赖它们
4. 新模块引入新软件时, 优先并入现有域环境(yaml 配方在 envs/*.yml), 不要新建环境
5. 调用方式: conda run -n <env> <tool> --no-capture-output

---

## 第一部分: 14 个功能域环境|Part 1: The 14 Domain Envs

### align — 比对与变异核心|Alignment & variant calling

```bash
conda run -n align <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `bcftools` | `bedtools` | `bgzip` |
| `bwa` | `freebayes` | `gatk` |
| `minimap2` | `samtools` | `tabix` |
| `wgsim` |  |  |

### pop — 群体遗传|Population genetics

```bash
conda run -n pop <工具> --no-capture-output
```

> easyhap: 纯 Python pip 安装(上游 Guangqi-He/EasyHap v1.0)

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `PopLDdecay` | `RAiSD` | `admixture` |
| `pixy` | `plink` | `poplddecay` |
| `treemix` | `vcftools` |  |
| `easyhap` |  |  |

### asm — 基因组组装|Genome assembly

```bash
conda run -n asm <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `canu` | `genomescope2` | `get_organelle_from_reads.py` |
| `hifiasm` | `jellyfish` | `kmc` |
| `kmc_tools` | `merqury` | `meryl` |
| `purge_dups` | `spades` | `spades.py` |
| `tidk` |  |  |

### hic — Hi-C 三维|Hi-C & 3D genome

```bash
conda run -n hic <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `filter_bam` | `haphic` | `juicer` |
| `matlock` | `pairtools` | `samblaster` |
| `yahs` |  |  |

### annot — 注释与功能预测|Annotation

```bash
conda run -n annot <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `TransDecoder.LongOrfs` | `TransDecoder.Predict` | `agat_convert_sp_gxf2gxf.pl` |
| `augustus` | `bam2hints` | `blastn` |
| `blastp` | `diamond` | `emapper.py` |
| `etraining` | `eviann` | `gff2gbSmallDNA.pl` |
| `gffcompare` | `gt` | `gtf2gff.pl` |
| `makeblastdb` | `miniprot` | `mmseqs` |
| `orthofinder` |  |  |

### repeat — 重复序列|Repeats

```bash
conda run -n repeat <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `BuildDatabase` | `LTR_retriever` | `RepeatMasker` |
| `RepeatModeler` | `TEsorter` | `ltr_finder_parallel` |
| `ltr_harvest_parallel` |  |  |

### rna — 转录组|RNA-seq

```bash
conda run -n rna <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `extract_exons.py` | `extract_splice_sites.py` | `fastp` |
| `gffread` | `hisat2` | `hisat2-build` |
| `infer_experiment.py` | `stringtie` |  |

### protein — 蛋白分析|Protein

```bash
conda run -n protein <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `embossversion` | `hmmsearch` | `meme` |
| `needle` | `needleall` | `phobius.pl` |
| `resistify` | `signalp` | `signalp6` |
| `tmhmm` | `tmhmm2` |  |

### phylo — 系统发育|Phylogenetics

```bash
conda run -n phylo <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `KaKs_Calculator` | `iqtree` | `iqtree2` |
| `mafft` | `nw_display` | `nw_reroot` |
| `raxml-ng` | `trimal` | `wgdi` |

### pan — 泛基因组|Pan-genome

```bash
conda run -n pan <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `delta-filter` | `kmindex` | `kmtricks` |
| `minigraph` | `nucmer` | `panman` |
| `pggb` | `show-coords` | `vg` |

### viz — 可视化|Visualization

```bash
conda run -n viz <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `pycirclize` | `samplot` |  |

### misc — 杂项工具|Misc utilities

```bash
conda run -n misc <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `axel` | `bbmap.sh` | `fastq-dump` |
| `fasterq-dump` | `iseq` | `pigz` |
| `prefetch` | `primer3_core` | `repair.sh` |
| `seqkit` |  |  |

### r — R 生态|R ecosystem

```bash
conda run -n r <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `mstmap` | `ggtree`(R包,进化树可视化)<br>`ape`(R包,树处理)<br>`ggplot2`(R包,画图) |

### busco — BUSCO 评估|BUSCO assessment

```bash
conda run -n busco <工具> --no-capture-output
```

| 工具<br>Tool | 工具<br>Tool | 工具<br>Tool |
|---|---|---|
| `busco` |  |  |

---

## 第二部分: 保留的独立环境(强依赖)|Part 2: Kept Standalone Envs

| 环境<br>Env | 关键软件<br>Key software |
|---|---|
| Augustus_v.3.5.0 | (见 backup yaml) |
| BUSCO_v.6.0.0 | (见 backup yaml) |
| Blast_v.2.16.0 | (见 backup yaml) |
| DeepBSA | (见 backup yaml) |
| EDTA_v.2.2.2 | (见 backup yaml) |
| edta_v.2.3.0 | (见 backup yaml) |
| EGAPx_v.0.4.0-alpha | (见 backup yaml) |
| GATK_v.4.6.2.0 | (见 backup yaml) |
| Genome_dedup | (见 backup yaml) |
| HiC-Pro_v3.1.0 | (见 backup yaml) |
| JCVI_v.1.5.6 | (见 backup yaml) |
| LTR_retriever_v.3.0.4 | (见 backup yaml) |
| Population_genetics | (见 backup yaml) |
| Rqtl | (见 backup yaml) |
| SubPhaser | (见 backup yaml) |
| WGCNA_v.1.73 | (见 backup yaml) |
| a-liner | (见 backup yaml) |
| adamixture_v.1.0.2 | (见 backup yaml) |
| alignoth | (见 backup yaml) |
| bbmap_v.39.81 | (见 backup yaml) |
| braker_v.3.0.8 | (见 backup yaml) |
| cphasing | (见 backup yaml) |
| deeptmhmm_v.1.0 | (见 backup yaml) |
| eviann_v.2.0.5 | (见 backup yaml) |
| fanc_v.0.9.23b | (见 backup yaml) |
| gctb | (见 backup yaml) |
| genomescope2_v.2.1.0 | (见 backup yaml) |
| genomescope_v.2.0.1 | (见 backup yaml) |
| genomesyn2 | (见 backup yaml) |
| hicpro | (见 backup yaml) |
| iseq_v.1.9.8 | (见 backup yaml) |
| jcvi_v.1.5.7 | (见 backup yaml) |
| juicer_v.1.6 | (见 backup yaml) |
| kmc_v.3.2.4 | (见 backup yaml) |
| kmeriaenv | (见 backup yaml) |
| ltr_finder_parallel_v.1.3 | (见 backup yaml) |
| ltr_harvest_parallel_v.1.2 | (见 backup yaml) |
| ltr_retriever_v.3.0.1 | (见 backup yaml) |
| merqury_v.1.3 | (见 backup yaml) |
| mga | (见 backup yaml) |
| panman_v.0.1.4 | (见 backup yaml) |
| pasa_v.2.5.3 | (见 backup yaml) |
| pggb_v.0.7.4 | (见 backup yaml) |
| picrust_v.2.6.3 | (见 backup yaml) |
| pixy_v.2.0.0 | (见 backup yaml) |
| plothic_v.1.0.0 | (见 backup yaml) |
| psvcp_v.1.0.1 | (见 backup yaml) |
| purge_dups_v.1.2.6 | (见 backup yaml) |
| pycirclize_v.1.10.1 | (见 backup yaml) |
| qiime_v.2024.10.1 | (见 backup yaml) |
| rMVP | (见 backup yaml) |
| rnaseq_val | (见 backup yaml) |
| signalp_v.3.0b | (见 backup yaml) |
| singularity_v.3.8.7 | (见 backup yaml) |
| telocomp | (见 backup yaml) |
| transdecoder_v.5.5.0 | (见 backup yaml) |
| treemix_v.1.13 | (见 backup yaml) |
| vcf2gwas_v.0.8.9 | (见 backup yaml) |
| vg_v.1.7.0 | (见 backup yaml) |
| selective_sweep | xpclr 1.1.2(已打 scipy/numpy 兼容补丁,补丁源码 ~/software/xpclr;pop 域环境内同名包为坏版本,勿用) |

---

## 第三部分: 项目运行环境|Part 3: Project Runtime Env

| 环境<br>Env | 内容<br>Content |
|---|---|
| biopytools | 项目 CLI(pip 装 biopytools) + kmtricks/rocksdb/modelscope — 用户日常命令行入口 |
| base | conda 基础环境(conda/mamba 本体) |

---

## 第四部分: 禁止使用|Part 4: Forbidden (待退役)

scripts/delete_list.txt 中的 153 个环境已被域环境取代或从未被模块引用, **不要**再调用它们, 也不要让新模块依赖它们; 如需其中的软件, 先并入对应域环境。

> 完整新旧对应关系见 docs/env_migration_map.md; 域环境重建配方见 envs/*.yml, 保留环境重建配方见 envs/legacy/*.yaml。
