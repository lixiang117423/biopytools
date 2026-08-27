# 模块分组草案 | Draft Module Grouping

> Phase 3 分组草案,**待确认**。边界与归属可调整;确认后用于首页分类表与 mkdocs.yml 导航。

> 范围:仅 CLI 注册的 210 个模块;非 CLI 模块不写文档。

## 群体遗传|Population(30)|30 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| admixture | admixture | ADMIXTURE群体结构分析|ADMIXTURE Population Structure Analysis |
| atomm | atomm | 双物种混合效应模型关联分析|Two-organism mixed model association (ATOMM) |
| cim | cim | R/qtl复合区间作图(CIM)分析|R/qtl Composite Interval Mapping (CIM) Analysis |
| deepbsa | deepbsa | DeepBSA批量分析工具|DeepBSA batch analysis tool |
| dsuite | dsuite | Dsuite D统计量分析工具|Dsuite D-statistics analysis tool |
| fastani | fastani | 全基因组ANI计算(FastANI)|Whole-genome ANI (FastANI) |
| fst | fst | Fst遗传分化计算工具|Fst genetic differentiation calculation tool |
| gctb | gctb | GCTB全基因组复杂性状贝叶斯分析|GCTB Genome-wide Complex Trait Bayesian Analysis |
| gemma_gwas | gemma-gwas | GEMMA GWAS批量分析工具|GEMMA GWAS batch analysis tool |
| gwas_gec | gwas-gec | GWAS基因组范围多重检验校正|GWAS genome-wide error correction |
| gwas_lambda | gwas-lambda | GWAS Lambda GC计算工具|GWAS Lambda GC calculation tool |
| gwas2gene | gwas2gene | GWAS候选基因筛选工具|GWAS candidate gene finder |
| janusx | janusx | JanusX GWAS和基因组选择分析|JanusX GWAS and Genomic Selection Analysis |
| kmeria | kmeria | K-mer GWAS全流程分析工具|K-mer GWAS Complete Pipeline Tool |
| ldblockshow | ldblockshow | 连锁不平衡热图分析|LD Heatmap Analysis |
| mixrace | mixrace | WGS混合小种检测|WGS mixed-race detection |
| ocbsa | ocbsa | BSA分析工具套件|BSA Analysis Tool Suite |
| pi | pi | 核苷酸多样性计算工具|Nucleotide diversity calculation tool |
| pi4gene | pi4gene | 基因分组核苷酸多样性计算|Calculate nucleotide diversity per gene group |
| pixy | pixy | Pixy群体遗传学统计工具|Pixy population genetics statistics tool |
| plinkgwas | plink-gwas | PLINK GWAS分析|PLINK GWAS analysis |
| poplddecay | poplddecay | 连锁不平衡衰减分析工具|Linkage disequilibrium decay analysis tool |
| rmvp | rmvp | rMVP GWAS批量分析工具（GLM/MLM/FarmCPU）|rMVP batch GWAS analysis (GLM/MLM/FarmCPU) |
| selective_sweep | selective-sweep | 选择性扫荡检测|Selective sweep detection |
| snp_index | snp-index | SNP index计算和分析工具|SNP index calculation and analysis tool |
| tassel_gwas | tassel-gwas | TASSEL GWAS分析工具|TASSEL GWAS analysis tool |
| treemix | treemix | TreeMix群体历史与基因流分析|TreeMix Population History & Gene Flow Analysis |
| vcf2gwas | vcf2gwas | vcf2gwas GWAS分析工具|vcf2gwas GWAS Analysis Tool |
| vcf2pca | vcf2pca | VCF主成分分析 (PCA, 默认PLINK后端支持SNP/INDEL)|VCF PCA (default PLINK backend, SNP/INDEL) |
| xpclr | xpclr | XP-CLR跨群体选择信号扫描|XP-CLR Cross-Population Selection Scan |

## 基因组组装|Assembly(24)|24 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| agp2table | agp2table | AGP转表格工具|AGP to table converter |
| assembly2agp | assembly2agp | Assembly文件转AGP格式工具|Assembly to AGP format converter |
| assembly_qc | assembly-qc | 基因组组装质量综合评估|Genome Assembly Quality Control |
| assembly_qv | assembly-qv | 装配质量QV值计算|Assembly Quality QV Calculation |
| assembly_stats | assembly-stats | 基因组装配统计|Genome Assembly Statistics |
| centier | centier | CentIER着丝粒鉴定工具|CentIER centromere identification tool |
| chr_rename | chr-rename | 基于minimap2的染色体重命名工具|Chromosome rename tool based on minimap2 |
| find_telomere | find-telomere | 端粒识别分析工具|Telomere identification analysis tool |
| gap_fill | gap-fill | TGS-GapCloser Gap填充工具|TGS-GapCloser gap filling tool |
| gap_stat | gap-stat | 基因组Gap统计工具|Genome gap statistics tool |
| genome_mount_rate | genome-mount-rate | 基因组挂载率统计|Genome mount rate statistics |
| genomeasm | genomeasm | 三代基因组组装流程|Third generation genome assembly pipeline |
| get_plastome | get-plastome | 叶绿体基因组组装工具|Chloroplast genome assembly tool |
| hicanu | hicanu | HiCanu基因组组装工具|HiCanu genome assembly tool |
| hifi_hic | hifi-hic | 使用HiFi和Hi-C数据进行基因组组装|Genome assembly using HiFi and Hi-C data |
| hifi_hic_workflow | hifi-hic-workflow | HiFi+Hi-C基因组组装与挂载完整流程|Complete HiFi+Hi-C Genome Assembly and Scaffolding Workflow |
| hifiasm | hifiasm | 运行hifiasm基因组组装|Run hifiasm genome assembly |
| mga | mga | MGA共识基因组组装(HiFi)|MGA consensus genome assembly (HiFi) |
| purge_dups | purge-dups | Purge_Dups基因组去冗余工具|Purge_Dups genome deduplication tool |
| ragtag | ragtag | RagTag基因组scaffolding工具|RagTag genome scaffolding tool |
| rename_chromosomes | rename-chromosomes | 染色体重命名工具|Chromosome renaming tool |
| rename_genome_id | rename-genome-id | 基因组ID重命名工具|Genome ID renaming tool |
| telocomp | telocomp | TeloComp端粒鉴定工具|TeloComp telomere identification tool |
| yahs | yahs | YaHS Hi-C scaffolding流程|YaHS Hi-C scaffolding pipeline |

## 基因组注释|Annotation(30)|30 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| annorefine | annorefine | BRAKER+查漏补漏端到端→整合GFF3(基因组+转录组+同源蛋白)|End-to-end: BRAKER + homology gap-filling → integrated GFF3 |
| braker | braker | BRAKER3基因组注释工具|BRAKER3 genome annotation tool |
| braker4phyto | braker4phyto | BRAKER3疫霉菌基因组注释(默认不屏蔽重复)|BRAKER3 Phytophthora annotation (no repeat masking by default) |
| busco | busco | BUSCO质量评估分析工具|BUSCO quality assessment tool |
| egapx_batch | egapx-batch | EGAPx批量运行配置生成工具|EGAPx batch run configuration generator |
| eggnog_mapper | eggnog-mapper | eggNOG功能注释(GO/KEGG/COG/CAZy/Pfam)|eggNOG functional annotation (GO/KEGG/COG/CAZy/Pfam) |
| eviann | eviann | EviAnn基因组注释工具|EviAnn genome annotation tool |
| func_anno | func-anno | 蛋白功能注释(IPS+eggnog→GO/KEGG标准表,衔接下游R)|Protein annotation (IPS+eggnog→GO/KEGG tables) |
| genebank2fasta | genebank2fasta | GenBank序列提取工具|GenBank sequence extraction tool |
| gene_density | gene-density | 基因密度计算(每窗口基因数与基因/Mb)|Gene density (genes/window and genes/Mb) |
| gene_rnaseq_check | gene-rnaseq-check | 候选基因RNA-seq转录验证|Candidate gene RNA-seq transcriptional validation |
| gene_table | gene-table | 基因信息+序列合并表(基因DNA+CDS+蛋白)|Gene info + sequence merged table (gene DNA + CDS + Protein) |
| geneinfo | geneinfo | 从GFF文件提取基因信息|Extract gene information from GFF file |
| gff_renamer | gff-renamer | GFF文件ID规范化工具|GFF file ID standardization tool |
| gffcompare | gffcompare | GFF/GTF文件两两比较分析|GFF/GTF pairwise comparison analysis |
| gffconverter | renamegff | GFF文件整理工具|GFF file organization tool |
| gtf2gff | gtf2gff | GTF到GFF文件转换工具|GTF to GFF file converter |
| insert_detection | insert-detection | 插入序列位点检测|Insert sequence insertion site detection |
| insert2locus | insert2locus | 转基因插入位点提取(步移+完整locus+验证)|Transgenic insertion locus extraction |
| interproscan | interproscan | InterProScan蛋白质功能注释|InterProScan protein function annotation |
| longestmrna | longest-mrna | 提取最长转录本|Extract longest transcript |
| nlr_annotator | nlr-annotator | NLR基因预测工具|NLR gene prediction tool |
| oomycete_anno | oomycete-anno | 疫霉菌基因组注释(T2T Augustus流程)|Oomycete genome annotation (T2T Augustus pipeline) |
| parse_gene_dna | parse-gene-dna | 基因DNA序列提取工具|Gene DNA sequence extraction tool |
| pep2genome | pep2genome | 蛋白质到基因组比对工具|Protein to genome alignment tool |
| phyto_effector | phyto-effector | Phytophthora效应子鉴定(rxlr/crn)|Phytophthora effector identification (rxlr/crn) |
| promoter_extractor | promoter-extractor | 启动子提取工具|Promoter extraction tool |
| resistify | resistify | Resistify NLR分析工具|Resistify NLR analysis tool |
| rnaseq_val | rnaseq-val | 转录组验证注释|Transcriptome validation for genome annotation |
| rxlr_scanner | rxlr-scanner | RxLR效应蛋白扫描工具|RxLR effector protein scanner |

## 转录组|Transcriptome(6)|6 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| dual_rnaseq | dual-rnaseq | 互作转录组分析|Dual RNA-seq Analysis |
| longrnaseq | longrnaseq | 三代转录组比对工具|Long RNA-seq alignment tool |
| rnabloom | rnabloom | RNA-Bloom转录组从头组装工具|RNA-Bloom de novo transcriptome assembly tool |
| rnaseq | rnaseq | RNA-seq表达定量流程|RNA-seq expression quantification pipeline |
| rnaseq2vcf | rnaseq2vcf | 转录组变异检测(到VCF)|RNA-seq variant calling (to VCF) |
| transcript_assembly | transcript-assembly | 转录本组装(FASTQ/BAM→GFF3,支持长读)|Transcript assembly (FASTQ/BAM→GFF3, long-read support) |

## 变异检测|Variants(23)|23 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| annovar | annovar | ANNOVAR变异注释|ANNOVAR Variant Annotation |
| bwa_gatk | bwa-gatk | 全基因组比对和变异检测|Whole genome alignment and variant detection |
| fastq2vcf_gtx | fastq2vcf-gtx | Fastq到VCF (GTX) 全流程分析|Fastq to VCF (GTX) pipeline |
| fastq2vcf_parabricks | fastq2vcf-parabricks | Fastq到VCF (Parabricks) 全流程分析|Fastq to VCF (Parabricks) pipeline |
| filter_snp_indel | filter-snp-indel | SNP和INDEL过滤工具|SNP and INDEL filtering tool |
| gatk_joint | gatk-joint | GATK Joint Genotyping工具|GATK Joint Genotyping tool |
| genome2sv | genome2sv | assembly-to-assembly SV calling (minimap2+svim-asm+SURVIVOR)|组装间结构变异检测 |
| gtx | gtx | 运行GTX WGS流程|Run GTX WGS pipeline |
| gtx_joint | gtx-joint | GTX Joint Calling命令生成工具|GTX Joint Calling command generator |
| indel_marker | indel-marker | 抗病/感病INDEL共显性标记开发|R/S INDEL codominant marker |
| indelpav | indelpav | INDEL PAV分析工具|INDEL PAV analysis tool |
| parabricks | parabricks | 基于GPU的全基因组流程|GPU-based whole genome pipeline |
| snp_region_gene | snp-region-gene | SNP区域基因提取工具|SNP Region Gene Extractor |
| swave | swave | Swave结构变异检测工具|Swave structural variant detection tool |
| vcf2gene | vcf2gene | VCF变异基因注释工具|VCF variant gene annotation tool |
| vcf2genotype | vcf2genotype | VCF基因型提取|VCF genotype extraction |
| vcf2pav | vcf2pav | VCF转PAV(Presence/Absence)矩阵|VCF to PAV (Presence/Absence) matrix |
| vcf_filter | vcf-filter | VCF文件筛选|VCF file filtering |
| vcf_merger | vcf-merger | VCF按染色体合并工具|VCF chromosome merge tool |
| vcf_renamer | vcf-renamer | VCF样品名称重命名工具|VCF sample name renaming tool |
| vcf_sample_hete | vcf-sample-hete | VCF样本基因型统计|VCF sample genotype statistics |
| vcf_sampler | vcf-sampler | VCF文件SNP抽样工具|VCF file SNP sampling tool |
| vcf_sequence | vcf-sequence | 从基因组和VCF提取序列|Extract sequences from genome and VCF |

## 泛基因组|Pan-genome(17)|17 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| cactus | cactus | Cactus泛基因组构建和分析工具|Cactus pangenome construction and analysis tool |
| genomesyn | genomesyn | 基因组共线性分析|Genome collinearity analysis |
| genomesyn2 | genomesyn2 | GenomeSyn2比较基因组学可视化工具|GenomeSyn2 comparative genomics visualization |
| microsynteny | microsynteny | 微观共线性分析工具|Microsynteny analysis tool |
| minigraph | minigraph | Minigraph泛基因组图构建和分析工具|Minigraph pangenome graph construction and analysis tool |
| ngenomesyn | ngenomesyn | NGenomeSyn可视化工具|NGenomeSyn visualization tool |
| orthofinder | orthofinder | OrthoFinder泛基因组分析工具包|OrthoFinder pan-genome analysis toolkit |
| pan_blocks | pan-blocks | 泛基因组Block构建工具|Pan-genome Block Construction Tool |
| pandepth | pandepth | PanDepth覆盖度计算工具|PanDepth coverage calculation tool |
| panedta | panedta | PanEDTA泛基因组转座子注释|PanEDTA Pan-genome TE annotation |
| panhite | panhite | panHiTE群体基因组TE分析|panHiTE pan-genome TE analysis |
| panman | panman | Panman泛基因组构建和分析工具|Panman pan-genome construction and analysis tool |
| panvar | panvar | 泛基因组变异分析|Pan-genome Variant Analysis |
| pggb | pggb | PGGB泛基因组图构建工具|PGGB Pangenome Graph Builder |
| psvcp | psvcp | PSVCP线性泛基因组构建(MUMmer+Assemblytics)|PSVCP linear pangenome construction |
| vg | vg | VG变异图分析工具（construct/index/giraffe/deconstruct）|VG variation graph analysis tool |
| wgdi | wgdi | WGDI比较基因组学分析工具|WGDI comparative genomics analysis tool |

## 三维基因组|3D&Chromatin(7)|7 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| allhic | allhic | 使用ALLHiC进行染色体挂载|Use ALLHiC for chromosome scaffolding |
| bismark | bismark | 全基因组甲基化分析|Whole genome methylation analysis |
| cphasing | cphasing | CPhasing基因组分相和挂载工具|CPhasing genome phasing and scaffolding tool |
| haphic | haphic | HapHiC基因组scaffolding工具|HapHiC genome scaffolding tool |
| hic_qc | hic-qc | Hi-C数据质量评估工具|Hi-C data quality assessment tool |
| subgenome_assign | subgenome-assign | 基于亲本比对的亚基因组归属|Subgenome assignment via parental alignment |
| subphaser | subphaser | SubPhaser异源多倍体亚基因组分离|SubPhaser subgenome phasing of allopolyploids |

## 系统发育|Phylogeny(17)|17 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| easyhap | easyhap | 区域单倍型分析(EasyHap)|Regional haplotype analysis (EasyHap) |
| genome2tree | genome2tree | 基因组目录免比对物种树(waster)|Alignment-free species tree from genome dir (waster) |
| reads2tree | reads2tree | fastq目录免组装物种树(WASTER)|De novo species tree from reads (WASTER) |
| splitstree6 | splitstree6 | SplitsTree6免比对建网/建树(VCF默认)|SplitsTree6 network/tree (VCF default) |
| iqtree | iqtree | IQ-TREE系统发育树分析工具|IQ-TREE phylogenetic tree analysis tool |
| kaks | kaks | Ka/Ks计算|Ka/Ks calculation |
| mafft_fasttree | mafft-fasttree | 系统发育树构建工具|Phylogenetic tree construction tool |
| msa | msa | 多序列比对分析工具|Multiple sequence alignment analysis tool |
| ncbi_taxo | ncbi-taxo | NCBI分类学注释工具|NCBI Taxonomy Annotation Tool |
| phylo_selector | phylo-selector | 系统发育树样品选择工具|Phylogenetic tree sample selection tool |
| phylo_trim | phylo-trim | 整合mafft-fasttree+trimal,出trimal前后两棵树|Integrate mafft-fasttree+trimal, before/after trees |
| raxml | raxml | RAxML系统发育树|RAxML phylogenetic tree |
| raxml_ng | raxml-ng | RAxML-NG系统发育树(maximum likelihood)|RAxML-NG phylogenetic tree (ML) |
| trimal | trimal | 多序列比对自动修剪(trimAl)|MSA automated trimming with trimAl |
| vcf2nj | vcf2nj | VCF构建NJ进化树|VCF NJ phylogenetic tree construction |
| vcf2phylip | vcf2phylip | VCF转phylip格式|VCF to phylip format conversion |
| vcf2tree | vcf2tree | VCF转系统发育树|VCF to phylogenetic tree |

## 蛋白质分析|Proteins(9)|9 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| deeploc | deeploc | DeepLoc 2.1蛋白质亚细胞定位预测工具|DeepLoc 2.1 protein subcellular localization prediction tool |
| deeptmhmm | deeptmhmm | DeepTMHMM跨膜螺旋/信号肽预测|DeepTMHMM TM helix & signal peptide prediction |
| hmmsearch | hmmsearch | HMMsearch结果处理工具|HMMsearch result processing tool |
| meme_parser | meme-parser | MEME Motif发现和解析工具|MEME Motif discovery and parser tool |
| phobius | phobius | Phobius跨膜拓扑+信号肽预测|Phobius TM topology & signal peptide prediction |
| predgpi | predgpi | PredGPI GPI锚定蛋白预测|PredGPI GPI-anchor prediction |
| protein_stats | protein-stats | Protein Stats理化性质分析工具|Protein Stats physicochemical properties analysis tool |
| signalp | signalp | SignalP 6.0信号肽预测工具|SignalP 6.0 signal peptide prediction tool |
| tmhmm | tmhmm | TMHMM跨膜螺旋预测|TMHMM transmembrane helix prediction |

## 重复序列|Repeats(5)|5 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| edta | edta | EDTA转座子注释|EDTA TE annotation |
| hite | hite | HiTE转座子检测与注释|HiTE transposon detection and annotation |
| lai | lai | LAI组装质量指数计算工具|LAI Assembly Index calculator |
| repeat_analyzer | repeat-analyzer | 重复序列分析模块|Repeat sequence analysis module |
| repeatmask | repeatmask | 重复序列屏蔽工具|Repeat masking tool |

## 可视化|Visualization(8)|8 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| aliner | aliner | a-liner共线性可视化pipeline(FASTA→minimap2→图)|a-liner synteny pipeline (FASTA->minimap2->plot) |
| bam_view | bam-view | BAM比对可视化工具|BAM Alignment Visualization Tool |
| hap_type | hap-type | 单倍型可视化工具|Haplotype visualization tool |
| hic_heatmap | hic-heatmap | Hi-C全基因组热图分析|Hi-C whole genome heatmap analysis |
| jcvi | jcvi | JCVI共线性分析工具集|JCVI Synteny Analysis Toolkit |
| msaviz | msaviz | MSA可视化工具（自动比对+可视化）|MSA Visualization Tool (Auto-align + Visualize) |
| plotsr | plotsr | 多基因组共线性可视化工具|Multi-genome synteny visualization tool |
| samplot | samplot | Samplot结构变异可视化工具|Samplot SV visualization tool |

## 微生物组|Microbiome(5)|5 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| faprotaxtax | faprotaxtax | FAPROTAX功能注释|FAPROTAX functional annotation |
| kraken2 | kraken2 | Kraken2宏基因组分类工具|Kraken2 metagenomic classification tool |
| mcyc | mcyc | 甲烷循环基因丰度分析工具|Methane cycle gene abundance analysis tool |
| picrust2 | picrust2 | PICRUSt2微生物群落功能丰度预测|PICRUSt2 functional abundance prediction |
| qiime2 | qiime2 | QIIME2微生物组多样性分析|QIIME2 Microbiome Diversity Analysis |

## 工具类|Utilities(34)|34 modules

| 模块 | 命令 | 说明 |
|------|------|------|
| allelic_genes | allelic-genes | [已废弃,请使用jcvi allelic]|[DEPRECATED, use jcvi allelic] |
| bam2fastq | bam2fastq | BAM to FASTQ批量转换工具|BAM to FASTQ batch conversion tool |
| bam_cov | bam-cov | BAM覆盖度统计|BAM Coverage Statistics |
| bam_stats | bam-stats | BAM文件批量统计分析|BAM File Batch Statistics Analysis |
| check_reads | check-reads | fastq完整性检查(gz/0字节/配对)|FASTQ integrity check (gzip/empty/pairing) |
| blast | blast | BLAST序列比对分析|BLAST Sequence Alignment Analysis |
| bwa | bwa | 全基因组比对工具|Whole genome alignment tool |
| coverage | coverage | BAM覆盖度分析工具|BAM coverage analysis tool |
| coverage_filter | coverage-filter | 基于覆盖度的序列质量过滤|Sequence quality filtering based on coverage |
| ncbi_datasets | ncbi-datasets | NCBI taxon基因组批量下载(datasets CLI)|NCBI taxon genome batch download (datasets CLI) |
| ena_downloader | ena-downloader | ENA数据下载工具|ENA data download tool |
| extract_reads | extract-reads | 基于contig-reads对应关系提取fastq reads|Extract fastq reads by contig-reads mapping |
| fastp | fastp | FASTQ数据质量控制|FASTQ data quality control |
| fastq_gc_filter | fastq-gc-filter | FASTQ文件GC含量和序列长度过滤|FASTQ file GC content and sequence length filtering |
| fastq_stats | fq-stats | FASTQ文件统计工具|FASTQ file statistics tool |
| fof | fof | FOF文件生成工具(样品名→路径映射表)|FOF file generator (sample→path mapping table) |
| get_link_from_CNCB | get-link-from-CNCB | 从CNCB批量获取测序数据下载链接|Batch download links from CNCB |
| iseq | iseq | 公共测序数据下载工具|Public sequencing data download tool |
| kmc | kmc | KMC k-mer统计和分析工具|KMC k-mer counting and analysis tool |
| kmertools | kmertools | K-mer工具集 - 建库、查询和分析|K-mer Tools - Build, Query and Analysis |
| minibwa | minibwa | Minibwa短读长比对（标准/Hi-C/BS-seq/长读）|Minibwa short-read alignment (standard/Hi-C/BS-seq/long-read) |
| minimap2 | minimap2 | Minimap2比对与区域提取|Minimap2 alignment and region extraction |
| needle_identity | needle-identity | 序列两两identity计算(EMBOSS needle)|Pairwise sequence identity (EMBOSS needle) |
| pair_fastq | pair-fastq | FASTQ配对修复工具|FASTQ pair fixing tool |
| parse_seq | parse-seq | 核酸或蛋白序列提取工具|Nucleotide or protein sequence extraction tool |
| primer3 | primer3 | Primer3引物设计工具|Primer3 Primer Design Tool |
| seq2genome | seq2genome | 序列到基因组比对工具（支持DNA/蛋白质自动检测）|Sequence to genome alignment tool (DNA/protein with auto-detection) |
| seq_extract | seq-extract | 序列提取(seqkit封装,自动识别ID/ID文件/BED)|Sequence extraction (seqkit wrapper, auto-detect ID/ID file/BED) |
| seq_len | seq-len | FASTA序列长度统计(文件/文件夹,合并+汇总)|FASTA sequence length statistics (file/folder, merged + summary) |
| smudgescope | smudgescope | GenomeScope2+Smudgeplot基因组评估工具|GenomeScope2 and Smudgeplot genome evaluation tool |
| split_fasta_id | split-fasta-id | 分割FASTA文件ID|Split FASTA file ID |
| sra2fastq | sra2fastq | SRA转FASTQ转换工具|SRA to FASTQ conversion tool |
| subseq | subseq | 序列子集提取工具|Sequence subset extraction tool |
| wgsim | wgsim | Wgsim基因组测序数据模拟|Wgsim genome sequencing simulation |

总覆盖|Total: 215/215;缺失|Missing: 无;重复|Duplicated: 无
