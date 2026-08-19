"""工具到功能域环境映射表|Tool → domain environment mapping

功能域环境合并方案(2026-08-16)的代码侧实现|Code-side implementation of the
domain consolidation plan. 每个工具映射到其所属功能域环境|Each tool maps to
its domain environment (~/miniforge3/envs/<domain>/bin/<tool>).

注意|Notes:
- 键为**二进制名**(与代码中默认路径 bin/<tool> 一致), 不是配置键名
  |Keys are binary names (matching bin/<tool> in default paths), not config keys
- Tier2/Tier3 例外环境(braker/fanc/deeptmhmm/EDTA/cphasing/jcvi/qiime/picrust/
  EGAPx/telocomp/rnaseq_val/singularity 等)不在此表, 继续用旧环境
  |Tier2/Tier3 exceptions are not here; they keep using legacy envs
- 域环境尚未安装的工具(存在性检查失败)自动回退旧默认路径, 无破坏
  |Tools not yet in the domain env fall back to legacy defaults automatically
"""

# 工具 -> 功能域|tool -> domain
TOOL_DOMAIN_MAP = {
    # ===== align 比对与变异核心 =====
    "gatk": "align",
    "bcftools": "align",
    "bgzip": "align",
    "tabix": "align",
    "samtools": "align",
    "bwa": "align",
    "freebayes": "align",
    "bedtools": "align",
    "wgsim": "align",
    "minimap2": "align",  # Genome_dedup/telocomp 归 align; cphasing 例外保留
    # ===== pop 群体遗传 =====
    "plink": "pop",
    "vcftools": "pop",
    "admixture": "pop",
    "treemix": "pop",
    "pixy": "pop",
    "PopLDdecay": "pop",
    "poplddecay": "pop",
    "RAiSD": "pop",
    "xpclr": "pop",
    "easyhap": "pop",
    # ===== asm 组装 =====
    "canu": "asm",
    "hifiasm": "asm",
    "kmc": "asm",
    "kmc_tools": "asm",
    "jellyfish": "asm",
    "merqury": "asm",
    "meryl": "asm",
    "purge_dups": "asm",
    "genomescope2": "asm",
    "tidk": "asm",
    "get_organelle_from_reads.py": "asm",
    "spades.py": "asm",
    "spades": "asm",
    # ===== hic =====
    "haphic": "hic",
    "pairtools": "hic",
    "yahs": "hic",
    "juicer": "hic",
    "matlock": "hic",
    "samblaster": "hic",
    "filter_bam": "hic",
    # ===== annot 注释 =====
    "augustus": "annot",
    "bam2hints": "annot",
    "etraining": "annot",
    "gff2gbSmallDNA.pl": "annot",
    "gtf2gff.pl": "annot",
    "agat_convert_sp_gxf2gxf.pl": "annot",
    "gffcompare": "annot",
    "miniprot": "annot",
    "TransDecoder.LongOrfs": "annot",
    "TransDecoder.Predict": "annot",
    "emapper.py": "annot",
    "orthofinder": "annot",
    "blastn": "annot",
    "blastp": "annot",
    "makeblastdb": "annot",
    "diamond": "annot",
    "mmseqs": "annot",
    "gt": "annot",
    "eviann": "annot",
    # ===== repeat 重复序列 =====
    "RepeatMasker": "repeat",
    "BuildDatabase": "repeat",
    "RepeatModeler": "repeat",
    "LTR_retriever": "repeat",
    "ltr_harvest_parallel": "repeat",
    "ltr_finder_parallel": "repeat",
    # ===== rna =====
    "hisat2": "rna",
    "hisat2-build": "rna",
    "extract_exons.py": "rna",
    "extract_splice_sites.py": "rna",
    "stringtie": "rna",
    "fastp": "rna",
    "gffread": "rna",
    "infer_experiment.py": "rna",
    # ===== protein 蛋白 =====
    "signalp6": "protein",
    "signalp": "protein",
    "resistify": "protein",
    "hmmsearch": "protein",
    "needle": "protein",
    "needleall": "protein",
    "embossversion": "protein",
    "meme": "protein",
    "phobius.pl": "protein",
    "tmhmm": "protein",
    "tmhmm2": "protein",
    # ===== phylo 系统发育 =====
    "iqtree": "phylo",
    "iqtree2": "phylo",
    "mafft": "phylo",
    "trimal": "phylo",
    "nw_reroot": "phylo",
    "nw_display": "phylo",
    "wgdi": "phylo",
    "KaKs_Calculator": "phylo",
    "raxml-ng": "phylo",
    # ===== pan 泛基因组 =====
    "pggb": "pan",
    "vg": "pan",
    "kmtricks": "pan",
    "kmindex": "pan",
    "panman": "pan",
    "nucmer": "pan",
    "delta-filter": "pan",
    "show-coords": "pan",
    "minigraph": "pan",
    # ===== viz 可视化 =====
    "samplot": "viz",
    "pycirclize": "viz",
    # ===== misc 杂项 =====
    "iseq": "misc",
    "primer3_core": "misc",
    "bbmap.sh": "misc",
    "repair.sh": "misc",
    "seqkit": "misc",
    "sra-tools": "misc",
    "fasterq-dump": "misc",
    "prefetch": "misc",
    "axel": "misc",
    "pigz": "misc",
    # ===== r R生态 =====
    "mstmap": "r",
    # ===== busco 独立 =====
    "busco": "busco",
}
