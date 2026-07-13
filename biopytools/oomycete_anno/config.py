"""疫霉菌基因组注释配置类|Oomycete Genome Annotation Config Class"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path, get_samtools_path, get_tool_path


@dataclass
class OomyceteAnnoConfig:
    """疫霉菌基因组注释配置|Oomycete genome annotation config.

    复刻 T2T 证据驱动 Augustus 流程, 支持基因组 + (可选)RNA-seq/同源蛋白/三代转录本。
    |Replicates the T2T evidence-driven Augustus pipeline; supports genome +
    (optional) RNA-seq / homologous proteins / long-read transcripts.
    """

    # ===== 必需参数|Required =====
    genome: str  # 基因组 FASTA|Genome FASTA
    species: str  # Augustus 物种名(etraining/config 用)|Augustus species name

    # ===== 可选证据(graceful degradation)|Optional evidence =====
    rnaseq_dirs: Optional[List[str]] = None  # 二代 RNA-seq 目录列表|short RNA-seq dirs
    prot_seq: Optional[str] = None  # 同源蛋白(Phase2)|homologous proteins (P2)
    isoseq: Optional[str] = None  # 三代转录本(Phase2)|long-read transcripts (P2)
    # 已知效应子蛋白(Phase3 救援): miniprot 全长比对当基因模型直接替换 Augustus 错注/漏注位点
    # |Known effectors (Phase3 rescue): use full-length miniprot alignments as gene models
    # to replace Augustus mis/missed annotations at effector loci
    effectors: Optional[str] = None

    # ===== 文件识别|File patterns =====
    read1_pattern: str = "_1.clean.fq.gz"  # R1 后缀|R1 suffix
    read2_pattern: str = "_2.clean.fq.gz"  # R2 后缀|R2 suffix
    # 链特异性: '' = 非链特异性(用户数据), 'RF'/'FR' = 链特异性
    # |Strandness: '' = unstranded (user data), 'RF'/'FR' = stranded
    rna_strandness: str = ""

    # ===== 输出/流程参数|Output & pipeline params =====
    output_dir: str = "./oomycete_anno_output"
    threads: int = 12
    soft_masking: bool = True  # RepeatMasker 软屏蔽|soft masking

    # ===== 重复屏蔽工具|Repeat masking tools =====
    repeatmodeler_bin: str = field(
        default_factory=lambda: get_tool_path(
            "repeatmodeler",
            "~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler",
            "REPEATMODELER_PATH",
        )
    )
    repeatmasker_bin: str = field(
        default_factory=lambda: get_tool_path(
            "repeatmasker",
            "~/miniforge3/envs/repeat_identiy/bin/RepeatMasker",
            "REPEATMASKER_PATH",
        )
    )
    build_database_bin: str = field(
        default_factory=lambda: get_tool_path(
            "build_database",
            "~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/BuildDatabase",
            "BUILDDATABASE_PATH",
        )
    )

    # ===== 二代 RNA-seq 比对工具|Short RNA-seq alignment tools =====
    hisat2_bin: str = field(
        default_factory=lambda: get_tool_path(
            "hisat2", "~/miniforge3/envs/RNA_Seq/bin/hisat2", "HISAT2_PATH"
        )
    )
    hisat2_build_bin: str = field(
        default_factory=lambda: get_tool_path(
            "hisat2_build",
            "~/miniforge3/envs/RNA_Seq/bin/hisat2-build",
            "HISAT2_BUILD_PATH",
        )
    )
    samtools_bin: str = field(default_factory=get_samtools_path)

    # ===== Augustus 工具|Augustus tools (host-side, Augustus_v.3.5.0) =====
    bam2hints_bin: str = field(
        default_factory=lambda: get_tool_path(
            "bam2hints",
            "~/miniforge3/envs/Augustus_v.3.5.0/bin/bam2hints",
            "BAM2HINTS_PATH",
        )
    )
    augustus_bin: str = field(
        default_factory=lambda: get_tool_path(
            "augustus",
            "~/miniforge3/envs/Augustus_v.3.5.0/bin/augustus",
            "AUGUSTUS_PATH",
        )
    )
    etraining_bin: str = field(
        default_factory=lambda: get_tool_path(
            "etraining",
            "~/miniforge3/envs/Augustus_v.3.5.0/bin/etraining",
            "ETRAINING_PATH",
        )
    )
    # Augustus 训练集转换脚本(GTF->GFF->GenBank; etraining 要 GenBank 位置参数)
    # |Augustus training-set converters (GTF->GFF->GenBank; etraining takes GenBank positional)
    gtf2gff_bin: str = field(
        default_factory=lambda: get_tool_path(
            "gtf2gff",
            "~/miniforge3/envs/Augustus_v.3.5.0/bin/gtf2gff.pl",
            "GTF2GFF_PATH",
        )
    )
    gff2gb_bin: str = field(
        default_factory=lambda: get_tool_path(
            "gff2gb",
            "~/miniforge3/envs/Augustus_v.3.5.0/bin/gff2gbSmallDNA.pl",
            "GFF2GB_PATH",
        )
    )
    # Augustus config 源(只读, 会被拷贝到输出目录供 etraining 写 species 参数)
    # |Augustus config source (read-only; copied to output so etraining can write species params)
    augustus_config_src: str = field(
        default_factory=lambda: get_tool_path(
            "augustus_config_src",
            "~/miniforge3/envs/Augustus_v.3.5.0/config",
            "AUGUSTUS_CONFIG_SRC",
        )
    )
    # hints 证据惩罚配置(intron+protein 用 MPE)|extrinsic cfg for intron+protein hints
    extrinsic_cfg: str = field(
        default_factory=lambda: expand_path(
            "~/miniforge3/envs/Augustus_v.3.5.0/config/extrinsic/extrinsic.MPE.cfg"
        )
    )

    # ===== GeneMark-ES/ET/EP+ (特殊: perl 环境与路径解耦)|GeneMark (special) =====
    gmes_petap_path: str = field(
        default_factory=lambda: get_tool_path(
            "gmes_petap",
            "~/software/GeneMark/gmes_linux_64_4/gmes_petap.pl",
            "GMES_PETAP_PATH",
        )
    )
    # GeneMark 的 perl 提供环境(系统 perl 缺 CPAN 模块)|perl provider env for GeneMark
    genemark_perl_env: str = "braker_v.3.0.8"

    # ===== Phase2 工具(首版未用, 留字段)|P2 tools (unused in MVP, fields reserved) =====
    gmap_bin: str = field(
        default_factory=lambda: get_tool_path(
            "gmap", "~/miniforge3/envs/pasa_v.2.5.3/bin/gmap", "GMAP_PATH"
        )
    )
    gmap_build_bin: str = field(
        default_factory=lambda: get_tool_path(
            "gmap_build", "~/miniforge3/envs/pasa_v.2.5.3/bin/gmap_build", "GMAP_BUILD_PATH"
        )
    )
    transdecoder_predict_bin: str = field(
        default_factory=lambda: get_tool_path(
            "transdecoder_predict",
            "~/miniforge3/envs/transdecoder_v.5.5.0/bin/TransDecoder.Predict",
            "TRANSDECODER_PREDICT_PATH",
        )
    )
    transdecoder_longorfs_bin: str = field(
        default_factory=lambda: get_tool_path(
            "transdecoder_longorfs",
            "~/miniforge3/envs/transdecoder_v.5.5.0/bin/TransDecoder.LongOrfs",
            "TRANSDECODER_LONGORFS_PATH",
        )
    )
    miniprot_bin: str = field(
        default_factory=lambda: get_tool_path(
            "miniprot", "~/miniforge3/envs/miniprot_v.0.18/bin/miniprot", "MINIPROT_PATH"
        )
    )
    # TransDecoder 坐标映射脚本(三代 ORF -> 基因组坐标)|TransDecoder ORF->genome mapper
    cdna_orf_to_genome_bin: str = field(
        default_factory=lambda: get_tool_path(
            "cdna_orf_to_genome",
            "~/miniforge3/envs/transdecoder_v.5.5.0/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl",
            "CDNA_ORF_TO_GENOME_PATH",
        )
    )

    # ===== LTR 注释工具(正交 TE)|LTR tools (orthogonal TE) =====
    gt_bin: str = field(
        default_factory=lambda: get_tool_path(
            "gt", "~/miniforge3/envs/genometools_v.1.6.5/bin/gt", "GT_PATH"
        )
    )
    ltr_retriever_bin: str = field(
        default_factory=lambda: get_tool_path(
            "ltr_retriever",
            "~/miniforge3/envs/ltr_retriever_v.3.0.1/bin/LTR_retriever",
            "LTR_RETRIEVER_PATH",
        )
    )

    # ===== 步骤跳过|Step skip flags =====
    skip_repeat: bool = False
    skip_rna: bool = False
    skip_iso: bool = False  # Phase2
    skip_protein: bool = False  # Phase2
    skip_ltr: bool = False  # Phase2
    skip_rescue: bool = False  # Phase3 效应子救援|effector rescue

    # ===== Phase3 效应子救援参数|Effector rescue params =====
    # miniprot 比对身份阈值(注: 全长判定靠 Target 起始=1 + stop_codon, 非高 identity)
    # |miniprot identity cutoff (full-length judged by Target start=1 + stop_codon, not high identity)
    rescue_min_identity: float = 0.85
    # Augustus 基因与效应子模型重叠占效应子模型长度 > 此比例 -> 视为冲突, 替换
    # |Augustus gene overlapping > this fraction of effector model length -> conflict, replace
    rescue_conflict_overlap: float = 0.50

    def __post_init__(self):
        """初始化后处理: 展开路径 + 建子目录|Post-init: expand paths, make subdirs."""
        self.working_dir = os.getcwd()

        # 展开所有含 ~ 的路径(规范§11.3.1)|expand all ~ paths (spec §11.3.1)
        self.genome = expand_path(self.genome)
        self.output_dir = expand_path(self.output_dir)
        self.repeatmodeler_bin = expand_path(self.repeatmodeler_bin)
        self.repeatmasker_bin = expand_path(self.repeatmasker_bin)
        self.build_database_bin = expand_path(self.build_database_bin)
        self.hisat2_bin = expand_path(self.hisat2_bin)
        self.hisat2_build_bin = expand_path(self.hisat2_build_bin)
        self.samtools_bin = expand_path(self.samtools_bin)
        self.bam2hints_bin = expand_path(self.bam2hints_bin)
        self.augustus_bin = expand_path(self.augustus_bin)
        self.etraining_bin = expand_path(self.etraining_bin)
        self.gtf2gff_bin = expand_path(self.gtf2gff_bin)
        self.gff2gb_bin = expand_path(self.gff2gb_bin)
        self.augustus_config_src = expand_path(self.augustus_config_src)
        self.extrinsic_cfg = expand_path(self.extrinsic_cfg)
        self.gmes_petap_path = expand_path(self.gmes_petap_path)
        self.gmap_bin = expand_path(self.gmap_bin)
        self.gmap_build_bin = expand_path(self.gmap_build_bin)
        self.transdecoder_predict_bin = expand_path(self.transdecoder_predict_bin)
        self.transdecoder_longorfs_bin = expand_path(self.transdecoder_longorfs_bin)
        self.miniprot_bin = expand_path(self.miniprot_bin)
        self.cdna_orf_to_genome_bin = expand_path(self.cdna_orf_to_genome_bin)
        self.gt_bin = expand_path(self.gt_bin)
        self.ltr_retriever_bin = expand_path(self.ltr_retriever_bin)

        # 展开可选证据路径|expand optional evidence paths
        if self.prot_seq:
            self.prot_seq = expand_path(self.prot_seq)
        if self.isoseq:
            self.isoseq = expand_path(self.isoseq)
        if self.effectors:
            self.effectors = expand_path(self.effectors)
        if self.rnaseq_dirs:
            self.rnaseq_dirs = [expand_path(d) for d in self.rnaseq_dirs]

        # 创建输出目录 + 子目录(by-step, 规范§12)|create output + subdirs
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.pipeline_info_dir = self.output_path / "00_pipeline_info"
        self.repeat_dir = self.output_path / "01_repeat_masking"
        self.rna_align_dir = self.output_path / "02_rna_align"
        self.iso_align_dir = self.output_path / "03_iso_align"
        self.protein_align_dir = self.output_path / "04_protein_align"
        self.hints_dir = self.output_path / "05_hints"
        self.training_dir = self.output_path / "06_training"
        self.augustus_dir = self.output_path / "07_augustus"
        self.ltr_dir = self.output_path / "08_ltr"
        self.rescue_dir = self.output_path / "09_effector_rescue"  # Phase3
        self.log_dir = self.output_path / "99_logs"
        for d in (
            self.pipeline_info_dir, self.repeat_dir, self.rna_align_dir,
            self.iso_align_dir, self.protein_align_dir, self.hints_dir,
            self.training_dir, self.augustus_dir, self.ltr_dir, self.rescue_dir,
            self.log_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def has_evidence(self) -> bool:
        """是否有任何证据数据|Whether any evidence is provided."""
        return bool(self.rnaseq_dirs or self.prot_seq or self.isoseq)

    def validate(self):
        """验证配置参数|Validate config."""
        errors = []

        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome not found: {self.genome}")

        if self.prot_seq and not os.path.exists(self.prot_seq):
            errors.append(f"蛋白文件不存在|Protein file not found: {self.prot_seq}")
        if self.isoseq and not os.path.exists(self.isoseq):
            errors.append(f"三代转录本不存在|Iso-seq file not found: {self.isoseq}")
        if self.effectors and not os.path.exists(self.effectors):
            errors.append(f"效应子蛋白文件不存在|Effectors file not found: {self.effectors}")
        if self.rnaseq_dirs:
            for i, d in enumerate(self.rnaseq_dirs):
                if not os.path.exists(d):
                    errors.append(
                        f"RNA-seq 目录不存在|RNA-seq dir not found [{i}]: {d}"
                    )

        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")

        if self.rna_strandness and self.rna_strandness not in ("FR", "RF"):
            errors.append(
                f"链特异性只能为 '' / 'FR' / 'RF'|Strandness must be ''/'FR'/'RF': "
                f"{self.rna_strandness}"
            )

        if errors:
            raise ValueError("\n".join(errors))
        return True
