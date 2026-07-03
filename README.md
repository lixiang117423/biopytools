# BioPyTools

A Python toolkit for bioinformatics analysis and computational biology.

一个用于生物信息学分析和计算生物学的Python工具包。

## 简介 | Overview

BioPyTools 是一个专为生物信息学研究设计的Python工具包，提供了一系列常用的生物数据分析功能。

BioPyTools is a Python toolkit designed for bioinformatics research, providing a series of commonly used biological data analysis functions.

## 系统要求 | Requirements

- Python >= 3.10
- NumPy >= 1.19.0
- Pandas >= 1.2.0
- Matplotlib >= 3.3.0
- pyfastx >= 0.8.4
- scikit-learn >= 1.0.0
- seaborn >= 0.11.0
- click >= 8.0.0

## 环境配置 | Environment Setup

相关的Conda 环境配置文件位于 [`conda_env/`](conda_env/) 目录下。

Conda environment files can be found in the [`conda_env/`](conda_env/) directory.

## 安装方法 | Installation

### 从 PyPI 安装 | Install from PyPI

```bash
pip install biopytools
```

### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biopytools.git
cd biopytools
pip install -e .

# or
pip install .
```

## 使用方法 | Usage

### 查看帮助 | Getting Help

```bash
biopytools -h
```

查看所有可用命令 | List all commands:

```bash
biopytools
```

查看特定命令的帮助 | Help for specific command:

```bash
biopytools <命令>-h    # 例如: biopytools fastp -h
```

## 模块文档 | Module Documentation

BioPyTools 提供175+个模块，按功能分为以下类别。

BioPyTools provides 175+ modules organized by functionality into the following categories.

### 数据下载与质控 | Data Download & QC

- [ena_downloader](./docs/ena_downloader.md) - ENA数据下载工具（`biopytools ena-downloader`）
- [fastp](./docs/fastp.md) - FASTQ数据质量控制（`biopytools fastp`）
- [fastq_gc_filter](./docs/fastq_gc_filter.md) - FASTQ文件GC含量和序列长度过滤（`biopytools fastq-gc-filter`）
- [fastq_stats](./docs/fastq_stats.md) - FASTQ文件统计工具（`biopytools fq-stats`）
- [get_link_from_CNCB](./docs/get_link_from_CNCB.md) - 从CNCB批量获取测序数据下载链接（`biopytools get-link-from-CNCB`）
- [iseq](./docs/iseq.md) - 公共测序数据下载工具（`biopytools iseq`）
- [pair_fastq](./docs/pair_fastq.md) - FASTQ配对修复工具（`biopytools pair-fastq`）
- [sra2fastq](./docs/sra2fastq.md) - SRA转FASTQ转换工具（`biopytools sra2fastq`）
- [wgsim](./docs/wgsim.md) - Wgsim基因组测序数据模拟（`biopytools wgsim`）

### 基因组组装 | Genome Assembly

- [genomeasm](./docs/genomeasm.md) - 三代基因组组装流程（`biopytools genomeasm`）
- [get_plastome](./docs/get_plastome.md) - 叶绿体基因组组装工具（`biopytools get-plastome`）
- [hicanu](./docs/hicanu.md) - HiCanu基因组组装工具（`biopytools hicanu`）
- [hifi_hic](./docs/hifi_hic.md) - 使用HiFi和Hi-C数据进行基因组组装（`biopytools hifi-hic`）
- [hifi_hic_workflow](./docs/hifi_hic_workflow.md) - HiFi+Hi-C基因组组装与挂载完整流程（`biopytools hifi-hic-workflow`）
- [hifiasm](./docs/hifiasm.md) - 运行hifiasm基因组组装（`biopytools hifiasm`）

### 组装评估与QC | Assembly QC

- [assembly_qc](./docs/assembly_qc.md) - 基因组组装质量综合评估（`biopytools assembly-qc`）
- [assembly_qv](./docs/assembly_qv.md) - 装配质量QV值计算（`biopytools assembly-qv`）
- [assembly_stats](./docs/assembly_stats.md) - 基因组装配统计（`biopytools assembly-stats`）
- [busco](./docs/busco.md) - BUSCO质量评估分析工具（`biopytools busco`）
- [centier](./docs/centier.md) - CentIER着丝粒鉴定工具（`biopytools centier`）
- [gap_stat](./docs/gapstat.md) - 基因组Gap统计工具（`biopytools gap-stat`）
- [hic_qc](./docs/hic_qc.md) - Hi-C数据质量评估工具（`biopytools hic-qc`）
- [lai](./docs/lai.md) - LAI组装质量指数计算工具（`biopytools lai`）

### Hi-C与挂载 | Hi-C & Scaffolding

- [allhic](./docs/allhic.md) - 使用ALLHiC进行染色体挂载（`biopytools allhic`）
- [cphasing](./docs/cphasing.md) - CPhasing基因组分相和挂载工具（`biopytools cphasing`）
- [find_telomere](./docs/find_telomere.md) - 端粒识别分析工具（`biopytools find-telomere`）
- [genome_mount_rate](./docs/genome_mount_rate.md) - 基因组挂载率统计（`biopytools genome-mount-rate`）
- [haphic](./docs/haphic.md) - HapHiC基因组scaffolding工具（`biopytools haphic`）
- [hic_heatmap](./docs/hic_heatmap.md) - Hi-C全基因组热图分析（`biopytools hic-heatmap`）
- [subphaser](./docs/subphaser.md) - SubPhaser异源多倍体亚基因组分离（`biopytools subphaser`）
- [telocomp](./docs/telocomp.md) - TeloComp端粒鉴定工具（`biopytools telocomp`）
- [yahs](./docs/yahs.md) - YaHS Hi-C scaffolding流程（`biopytools yahs`）

### 基因组后处理与工具 | Polishing & Utilities

- [agp2table](./docs/agp2table.md) - AGP转表格工具（`biopytools agp2table`）
- [assembly2agp](./docs/assembly2agp.md) - Assembly文件转AGP格式工具（`biopytools assembly2agp`）
- [chr_rename](./docs/chr_rename.md) - 基于minimap2的染色体重命名工具（`biopytools chr-rename`）
- [gap_fill](./docs/gap_fill.md) - TGS-GapCloser Gap填充工具（`biopytools gap-fill`）
- [genebank2fasta](./docs/genebank2fasta.md) - GenBank序列提取工具（`biopytools genebank2fasta`）
- [gffconverter](./docs/gffconverter.md) - GFF文件整理工具（`biopytools renamegff`）
- [gtf2gff](./docs/gtf2gff.md) - GTF到GFF文件转换工具（`biopytools gtf2gff`）
- [promoter_extractor](./docs/promoter_extractor.md) - 启动子提取工具（`biopytools promoter-extractor`）
- [purge_dups](./docs/purge_dups.md) - Purge_Dups基因组去冗余工具（`biopytools purge-dups`）
- [ragtag](./docs/ragtag.md) - RagTag基因组scaffolding工具（`biopytools ragtag`）
- [rename_chromosomes](./docs/rename_chromosomes.md) - 染色体重命名工具（`biopytools rename-chromosomes`）
- [rename_genome_id](./docs/rename_genome_id.md) - 基因组ID重命名工具（`biopytools rename-genome-id`）
- [split_fasta_id](./docs/split_fasta_id.md) - 分割FASTA文件ID（`biopytools split-fasta-id`）
- [subseq](./docs/subseq.md) - 序列子集提取工具（`biopytools subseq`）

### 比对与BAM处理 | Alignment & BAM

- [bam2fastq](./docs/bam2fastq.md) - BAM to FASTQ批量转换工具（`biopytools bam2fastq`）
- [bam_cov](./docs/bam_coverage_stats.md) - BAM覆盖度统计（`biopytools bam-cov`）
- [bam_stats](./docs/bam_stats.md) - BAM文件批量统计分析（`biopytools bam-stats`）
- [bam_view](./docs/bam_view.md) - BAM比对可视化工具（`biopytools bam-view`）
- [blast](./docs/blast.md) - BLAST序列比对分析（`biopytools blast`）
- [bwa](./docs/bwa.md) - 全基因组比对工具（`biopytools bwa`）
- [bwa_gatk](./docs/bwa.md) - 全基因组比对和变异检测（`biopytools bwa-gatk`）
- [coverage](./docs/coverage_filter.md) - BAM覆盖度分析工具（`biopytools coverage`）
- [coverage_filter](./docs/coverage_filter.md) - 基于覆盖度的序列质量过滤（`biopytools coverage-filter`）
- [extract_reads](./docs/extract_reads.md) - 基于contig-reads对应关系提取fastq reads（`biopytools extract-reads`）
- [longestmrna](./docs/longest_mrna.md) - 提取最长转录本（`biopytools longest-mrna`）
- [mafft_fasttree](./docs/mafft_fasttree.md) - 系统发育树构建工具（`biopytools mafft-fasttree`）
- [minibwa](./docs/minibwa.md) - Minibwa短读长比对（标准/Hi-C/BS-seq/长读）（`biopytools minibwa`）
- [minimap2](./docs/minimap2.md) - Minimap2比对与区域提取（`biopytools minimap2`）
- [msa](./docs/msaviz.md) - 多序列比对分析工具（`biopytools msa`）
- [msaviz](./docs/msaviz.md) - MSA可视化工具（自动比对+可视化）（`biopytools msaviz`）
- [parse_gene_dna](./docs/parse_gene_dna.md) - 基因DNA序列提取工具（`biopytools parse-gene-dna`）
- [parse_seq](./docs/parse_sequence_vcf.md) - 核酸或蛋白序列提取工具（`biopytools parse-seq`）
- [pep2genome](./docs/pep2genome.md) - 蛋白质到基因组比对工具（`biopytools pep2genome`）
- [samplot](./docs/samplot.md) - Samplot结构变异可视化工具（`biopytools samplot`）
- [seq2genome](./docs/seq2genome.md) - 序列到基因组比对工具（支持DNA/蛋白质自动检测）（`biopytools seq2genome`）

### 变异检测与VCF | Variant & VCF

- [annovar](./docs/annovar.md) - ANNOVAR变异注释（`biopytools annovar`）
- [fastq2vcf_gtx](./docs/fastq2vcf_gtx.md) - Fastq到VCF (GTX) 全流程分析（`biopytools fastq2vcf-gtx`）
- [fastq2vcf_parabricks](./docs/fastq2vcf_parabricks.md) - Fastq到VCF (Parabricks) 全流程分析（`biopytools fastq2vcf-parabricks`）
- [filter_snp_indel](./docs/filter_snp_indel.md) - SNP和INDEL过滤工具（`biopytools filter-snp-indel`）
- [gatk_joint](./docs/gatk_joint.md) - GATK Joint Genotyping工具（`biopytools gatk-joint`）
- [geneinfo](./docs/parse_gene_info.md) - 从GFF文件提取基因信息（`biopytools geneinfo`）
- [gffcompare](./docs/gffcompare.md) - GFF/GTF文件两两比较分析（`biopytools gffcompare`）
- [gtx](./docs/gtx-joint.md) - 运行GTX WGS流程（`biopytools gtx`）
- [gtx_joint](./docs/gtx-joint.md) - GTX Joint Calling命令生成工具（`biopytools gtx-joint`）
- [hap_type](./docs/hap_type.md) - 单倍型可视化工具（`biopytools hap-type`）
- [indelpav](./docs/indelpav.md) - INDEL PAV分析工具（`biopytools indelpav`）
- [insert_detection](./docs/insert_detection.md) - 插入序列位点检测（`biopytools insert-detection`）
- [parabricks](./docs/parabricks.md) - 基于GPU的全基因组流程（`biopytools parabricks`）
- [swave](./docs/swave.md) - Swave结构变异检测工具（`biopytools swave`）
- [vcf2gene](./docs/vcf2gene.md) - VCF变异基因注释工具（`biopytools vcf2gene`）
- [vcf2genotype](./docs/vcf2genotype.md) - VCF基因型提取（`biopytools vcf2genotype`）
- [vcf_filter](./docs/vcf_filter.md) - VCF文件筛选（`biopytools vcf-filter`）
- [vcf_merger](./docs/vcf_merger.md) - VCF按染色体合并工具（`biopytools vcf-merger`）
- [vcf_renamer](./docs/vcf_renamer.md) - VCF样品名称重命名工具（`biopytools vcf-renamer`）
- [vcf_sample_hete](./docs/vcf_sample_hete.md) - VCF样本基因型统计（`biopytools vcf-sample-hete`）
- [vcf_sampler](./docs/vcf_sampler.md) - VCF文件SNP抽样工具（`biopytools vcf-sampler`）
- [vcf_sequence](./docs/vcf_sequence.md) - 从基因组和VCF提取序列（`biopytools vcf-sequence`）
- [vg](./docs/vg.md) - VG变异图分析工具（construct/index/giraffe/deconstruct）（`biopytools vg`）

### 泛基因组 | Pan-genome

- [cactus](./docs/cactus.md) - Cactus泛基因组构建和分析工具（`biopytools cactus`）
- [minigraph](./docs/minigraph.md) - Minigraph泛基因组图构建和分析工具（`biopytools minigraph`）
- [pan_blocks](./docs/pan_blocks.md) - 泛基因组Block构建工具（`biopytools pan-blocks`）
- [pandepth](./docs/pandepth.md) - PanDepth覆盖度计算工具（`biopytools pandepth`）
- [panman](./docs/panman.md) - Panman泛基因组构建和分析工具（`biopytools panman`）
- [panvar](./docs/panvar.md) - 泛基因组变异分析（`biopytools panvar`）
- [pggb](./docs/pggb.md) - PGGB泛基因组图构建工具（`biopytools pggb`）

### 注释与功能预测 | Annotation

- [braker](./docs/braker.md) - BRAKER3基因组注释工具（`biopytools braker`）
- [deeploc](./docs/deeploc.md) - DeepLoc 2.1蛋白质亚细胞定位预测工具（`biopytools deeploc`）
- [egapx_batch](./docs/egapx_batch.md) - EGAPx批量运行配置生成工具（`biopytools egapx-batch`）
- [eviann](./docs/eviann.md) - EviAnn基因组注释工具（`biopytools eviann`）
- [gff_renamer](./docs/gff_renamer.md) - GFF文件ID规范化工具（`biopytools gff-renamer`）
- [hmmsearch](./docs/hmmsearch.md) - HMMsearch结果处理工具（`biopytools hmmsearch`）
- [interproscan](./docs/interproscan.md) - InterProScan蛋白质功能注释（`biopytools interproscan`）
- [meme_parser](./docs/meme_parser.md) - MEME Motif发现和解析工具（`biopytools meme-parser`）
- [nlr_annotator](./docs/nlr_annotator.md) - NLR基因预测工具（`biopytools nlr-annotator`）
- [signalp](./docs/signalp.md) - SignalP 6.0信号肽预测工具（`biopytools signalp`）
- [tmhmm](./docs/tmhmm.md) - TMHMM跨膜螺旋预测（`biopytools tmhmm`）
- [deeptmhmm](./docs/deeptmhmm.md) - DeepTMHMM跨膜螺旋/信号肽预测（`biopytools deeptmhmm`）
- [transcript_assembly](./docs/transcript_assembly.md) - 转录本从头组装(HISAT2+StringTie)（`biopytools transcript-assembly`）

### 转座子与重复序列 | TE & Repeats

- [edta](./docs/edta.md) - EDTA转座子注释（`biopytools edta`）
- [hite](./docs/hite.md) - HiTE转座子检测与注释（`biopytools hite`）
- [panedta](./docs/edta.md) - PanEDTA泛基因组转座子注释（`biopytools panedta`）
- [panhite](./docs/hite.md) - panHiTE群体基因组TE分析（`biopytools panhite`）
- [repeat_analyzer](./docs/repeat_analyzer.md) - 重复序列分析模块（`biopytools repeat-analyzer`）
- [repeatmask](./docs/repeatmask.md) - 重复序列屏蔽工具（`biopytools repeatmask`）

### RNA-seq与转录组 | RNA-seq

- [dual_rnaseq](./docs/dual_rnaseq.md) - 互作转录组分析（`biopytools dual-rnaseq`）
- [gene_rnaseq_check](./docs/rnaseq.md) - 候选基因RNA-seq转录验证（`biopytools gene-rnaseq-check`）
- [longrnaseq](./docs/longrnaseq.md) - 三代转录组比对工具（`biopytools longrnaseq`）
- [rnabloom](./docs/rnabloom.md) - RNA-Bloom转录组从头组装工具（`biopytools rnabloom`）
- [rnaseq](./docs/rnaseq.md) - RNA-seq表达定量流程（`biopytools rnaseq`）
- [rnaseq_val](./docs/rnaseq.md) - 转录组验证注释（`biopytools rnaseq-val`）

### 共线性与比较基因组 | Synteny

- [genomesyn](./docs/genomesyn2.md) - 基因组共线性分析（`biopytools genomesyn`）
- [genomesyn2](./docs/genomesyn2.md) - GenomeSyn2比较基因组学可视化工具（`biopytools genomesyn2`）
- [jcvi](./docs/jcvi.md) - JCVI共线性分析工具集（`biopytools jcvi`）
- [kaks](./docs/kaks.md) - Ka/Ks计算（`biopytools kaks`）
- [microsynteny](./docs/microsynteny.md) - 微观共线性分析工具（`biopytools microsynteny`）
- [ngenomesyn](./docs/ngenomesyn.md) - NGenomeSyn可视化工具（`biopytools ngenomesyn`）
- [orthofinder](./docs/orthofinder.md) - OrthoFinder泛基因组分析工具包（`biopytools orthofinder`）
- [plotsr](./docs/plotsr.md) - 多基因组共线性可视化工具（`biopytools plotsr`）
- [wgdi](./docs/wgdi.md) - WGDI比较基因组学分析工具（`biopytools wgdi`）
- [aliner](./docs/aliner.md) - a-liner共线性可视化pipeline（FASTA→minimap2→图）（`biopytools aliner`）

### 系统发育 | Phylogenetics

- [iqtree](./docs/iqtree.md) - IQ-TREE系统发育树分析工具（`biopytools iqtree`）
- [phylo_selector](./docs/phylo_selector.md) - 系统发育树样品选择工具（`biopytools phylo-selector`）
- [raxml](./docs/raxml.md) - RAxML系统发育树（`biopytools raxml`）
- [vcf2nj](./docs/vcf2nj.md) - VCF构建NJ进化树（`biopytools vcf2nj`）
- [vcf2phylip](./docs/vcf2phylip.md) - VCF转phylip格式（`biopytools vcf2phylip`）

### 群体遗传 | Population Genetics

- [admixture](./docs/admixture.md) - ADMIXTURE群体结构分析（`biopytools admixture`）
- [dsuite](./docs/dsuite.md) - Dsuite D统计量分析工具（`biopytools dsuite`）
- [fst](./docs/fst.md) - Fst遗传分化计算工具（`biopytools fst`）
- [pi](./docs/pixy.md) - 核苷酸多样性计算工具（`biopytools pi`）
- [pi4gene](./docs/pi4gene.md) - 基因分组核苷酸多样性计算（`biopytools pi4gene`）
- [pixy](./docs/pixy.md) - Pixy群体遗传学统计工具（`biopytools pixy`）
- [poplddecay](./docs/poplddecay.md) - 连锁不平衡衰减分析工具（`biopytools poplddecay`）
- [smudgescope](./docs/smudgescope.md) - GenomeScope2+Smudgeplot基因组评估工具（`biopytools smudgescope`）
- [treemix](./docs/treemix.md) - TreeMix群体历史与基因流分析（`biopytools treemix`）
- [vcf2pca](./docs/vcf2pca.md) - VCF主成分分析 (PCA)（`biopytools vcf2pca`）
- [vcf_pca](./docs/vcf_pca.md) - VCF主成分分析 (PCA)（`biopytools vcf-pca`）

### GWAS与BSA | GWAS & BSA

- [atomm](./docs/atomm.md) - 双物种混合效应模型关联分析（`biopytools atomm`）
- [cim](./docs/cim.md) - R/qtl复合区间作图(CIM)分析（`biopytools cim`）
- [deepbsa](./docs/deepbsa.md) - DeepBSA批量分析工具（`biopytools deepbsa`）
- [gctb](./docs/gctb.md) - GCTB全基因组复杂性状贝叶斯分析（`biopytools gctb`）
- [gemma_gwas](./docs/gemma-gwas.md) - GEMMA GWAS批量分析工具（`biopytools gemma-gwas`）
- [gwas2gene](./docs/gwas2gene.md) - GWAS候选基因筛选工具（`biopytools gwas2gene`）
- [gwas_gec](./docs/gwas_gec.md) - GWAS基因组范围多重检验校正（`biopytools gwas-gec`）
- [gwas_lambda](./docs/gwas_lambda.md) - GWAS Lambda GC计算工具（`biopytools gwas-lambda`）
- [janusx](./docs/janusx.md) - JanusX GWAS和基因组选择分析（`biopytools janusx`）
- [ldblockshow](./docs/ldblockshow.md) - 连锁不平衡热图分析（`biopytools ldblockshow`）
- [ocbsa](./docs/ocbsa.md) - BSA分析工具套件（`biopytools ocbsa`）
- [plinkgwas](./docs/plinkgwas.md) - PLINK GWAS分析（`biopytools plink-gwas`）
- [rmvp](./docs/rmvp.md) - rMVP GWAS批量分析工具（GLM/MLM/FarmCPU）（`biopytools rmvp`）
- [snp_index](./docs/snp_index.md) - SNP index计算和分析工具（`biopytools snp-index`）
- [snp_region_gene](./docs/snp_region_gene.md) - SNP区域基因提取工具（`biopytools snp-region-gene`）
- [tassel_gwas](./docs/tassel_gwas.md) - TASSEL GWAS分析工具（`biopytools tassel-gwas`）
- [vcf2gwas](./docs/vcf2gwas.md) - vcf2gwas GWAS分析工具（`biopytools vcf2gwas`）

### 微生物组与k-mer | Microbiome & k-mer

- [kmc](./docs/kmc.md) - KMC k-mer统计和分析工具（`biopytools kmc`）
- [kmeria](./docs/kmeria.md) - K-mer GWAS全流程分析工具（`biopytools kmeria`）
- [kmertools](./docs/kmertools.md) - K-mer工具集 - 建库、查询和分析（`biopytools kmertools`）
- [kraken2](./docs/kraken2.md) - Kraken2宏基因组分类工具（`biopytools kraken2`）
- [mcyc](./docs/mcyc.md) - 甲烷循环基因丰度分析工具（`biopytools mcyc`）
- [ncbi_taxo](./docs/ncbi_taxo.md) - NCBI分类学注释工具（`biopytools ncbi-taxo`）
- [picrust2](./docs/picrust2.md) - PICRUSt2微生物群落功能丰度预测（`biopytools picrust2`）

### 甲基化 | Methylation

- [bismark](./docs/bismark.md) - 全基因组甲基化分析（`biopytools bismark`）

### 效应子与抗病 | Effectors

- [phyto_effector](./docs/phyto_effector.md) - Phytophthora效应子鉴定(rxlr/crn)（`biopytools phyto-effector`）
- [resistify](./docs/resistify_parser.md) - Resistify NLR分析工具（`biopytools resistify`）
- [rxlr_scanner](./docs/rxlr_scanner.md) - RxLR效应蛋白扫描工具（`biopytools rxlr-scanner`）

### 其他工具 | Other Utilities

- [primer3](./docs/primer3.md) - Primer3引物设计工具（`biopytools primer3`）
- [protein_stats](./docs/protein_stats.md) - Protein Stats理化性质分析工具（`biopytools protein-stats`）
- [gene_density](./docs/gene_density.md) - 基因密度计算（每窗口基因数+基因/Mb）（`biopytools gene-density`）

## 许可证 | License

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 作者信息 | Author

**李详 (Xiang Li)**
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 致谢 | Acknowledgments

- 感谢所有为本项目做出贡献的开发者 | Thanks to all developers who contributed to this project
- 感谢开源社区的支持 | Thanks to the open source community for support
- AI 辅助开发 | AI-assisted development: Claude Code (Anthropic) + GLM 5.2

## 问题反馈 | Issues

如果遇到问题或有建议，请在GitHub上提交issue：

If you encounter problems or have suggestions, please submit an issue on GitHub:

https://github.com/lixiang117423/biopytools/issues
