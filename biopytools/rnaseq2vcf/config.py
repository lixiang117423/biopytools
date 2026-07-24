"""rnaseq2vcf 配置管理|Rnaseq2vcf configuration management (VCF-only, no annotation)"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class Rnaseq2vcfConfig:
    """转录组变异检测配置(到 VCF 为止,不含注释)|
    RNA-seq variant calling configuration (up to VCF, no annotation)"""

    # 必需|Required
    ref_genome_fa: str
    output_dir: str

    # 可选: HISAT2 剪接位点(提供则用,否则 HISAT2 de novo 发现 junction)|
    # Optional: HISAT2 splice sites (used if provided; else HISAT2 discovers junctions de novo)
    gff3_file: Optional[str] = None

    # 输入目录二选一|Input dir (one of two)
    raw_fastq_dir: Optional[str] = None
    clean_fastq_dir: Optional[str] = None

    # 共享目录(默认在 output_dir 下)|Shared dirs (default under output_dir)
    genome_index_dir: Optional[str] = None

    # 工具路径(支持~与环境变量)|Tool paths (~ and env-var overridable)
    hisat2_path: str = field(
        default_factory=lambda: get_tool_path('hisat2', '~/miniforge3/envs/RNA_Seq/bin/hisat2', 'HISAT2_PATH'))
    fastp_path: str = field(
        default_factory=lambda: get_tool_path('fastp', '~/miniforge3/envs/RNA_Seq/bin/fastp', 'FASTP_PATH'))
    samtools_path: str = field(
        default_factory=lambda: get_tool_path('samtools', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools', 'SAMTOOLS_PATH'))
    bcftools_path: str = field(
        default_factory=lambda: get_tool_path('bcftools', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/bcftools', 'BCFTOOLS_PATH'))
    gatk_path: str = field(
        default_factory=lambda: get_tool_path('gatk', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/gatk', 'GATK_PATH'))

    # 参数|Params
    threads: int = 12
    min_conf: int = 20
    fs_threshold: float = 30.0
    qd_threshold: float = 2.0
    cluster_window: int = 35
    cluster_size: int = 3
    read1_pattern: str = "_1.fq.gz"
    read2_pattern: str = "_2.fq.gz"

    # 流程控制|Control
    enable_checkpoint: bool = True
    dry_run: bool = False
    force: bool = False
    skip_qc: bool = False
    step: Optional[int] = None  # 0=仅建索引|index only; None/1-4=全流程(断点续传)|full pipeline
    log_file: Optional[str] = None
    log_level: str = "INFO"

    def __post_init__(self):
        """展开~、规范化路径、预计算共享目录与索引名|Expand ~, normalize, precompute shared dirs + index name"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        self.ref_genome_fa = os.path.normpath(os.path.abspath(expand_path(self.ref_genome_fa)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        if self.gff3_file:
            self.gff3_file = os.path.normpath(os.path.abspath(expand_path(self.gff3_file)))

        for attr in ('raw_fastq_dir', 'clean_fastq_dir', 'genome_index_dir', 'log_file'):
            v = getattr(self, attr)
            if v:
                setattr(self, attr, os.path.normpath(os.path.abspath(expand_path(v))))

        for attr in ('hisat2_path', 'fastp_path', 'samtools_path', 'bcftools_path', 'gatk_path'):
            setattr(self, attr, expand_path(getattr(self, attr)))

        # 索引命名(取基因组 basename)|index naming from genome basename
        self.genome_name = os.path.splitext(os.path.basename(self.ref_genome_fa))[0]

        if not self.genome_index_dir:
            self.genome_index_dir = os.path.join(self.output_dir, "genome_index")
        os.makedirs(self.genome_index_dir, exist_ok=True)

        # tmp 目录(§12.4.1)|tmp dir under output_dir
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        os.makedirs(self.tmp_dir, exist_ok=True)

    def validate(self):
        """校验配置|Validate configuration"""
        errors = []
        if not os.path.exists(self.ref_genome_fa):
            errors.append(f"参考基因组不存在|Reference genome does not exist: {self.ref_genome_fa}")
        if self.gff3_file and not os.path.exists(self.gff3_file):
            errors.append(f"GFF3 文件不存在|GFF3 file does not exist: {self.gff3_file}")
        if not self.raw_fastq_dir and not self.clean_fastq_dir:
            errors.append("必须提供 -i/--input(原始) 或 --clean-fastq-dir 之一|"
                          "Must provide -i/--input (raw) or --clean-fastq-dir")
        if self.threads <= 0:
            errors.append("线程数必须为正|Thread count must be positive")
        if self.step is not None and self.step not in (0, 1, 2, 3, 4):
            errors.append(f"无效步骤|Invalid step: {self.step} (应为 0-4|should be 0-4)")
        if errors:
            raise ValueError("\n".join(errors))
        return True
