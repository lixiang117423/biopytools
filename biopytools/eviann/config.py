"""
EviAnn配置模块|EviAnn Configuration Module
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import os


@dataclass
class EviAnnConfig:
    """EviAnn配置类|EviAnn Configuration Class"""

    # 必需参数|Required parameters
    genome: str

    # 数据输入参数（至少需要一个）|Data input parameters (at least one required)
    rnaseq: Optional[str] = None  # RNA-seq描述文件（已弃用，建议使用short_reads和long_reads）|RNA-seq description file (deprecated, use short_reads and long_reads)
    short_reads: Optional[str] = None  # 二代转录组数据文件或目录|Short-read RNA-seq file or directory
    long_reads: Optional[str] = None  # 三代转录组数据文件或目录|Long-read RNA-seq file or directory
    transcripts: Optional[str] = None
    proteins: Optional[str] = None

    # 可选参数|Optional parameters
    uniprot: Optional[str] = None
    threads: int = 12
    max_intron: Optional[int] = None
    ploidy: int = 2
    cds_gff: Optional[str] = None
    lncrna_tpm: float = 1.0
    partial: bool = False
    functional: bool = False
    mito_contigs: Optional[str] = None
    extra_gff: Optional[str] = None
    debug: bool = False
    verbose: bool = False

    # 软件路径|Software path
    eviann_path: str = '~/miniforge3/envs/eviann_v.2.0.5'

    # 输出目录|Output directory (使用基因组文件名作为前缀)
    output_prefix: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        from .utils import expand_path
        self.eviann_path = expand_path(self.eviann_path)
        self.genome = os.path.abspath(self.genome)

        # 设置输出前缀|Set output prefix
        if self.output_prefix is None:
            # 使用基因组文件名作为前缀|Use genome filename as prefix
            genome_name = Path(self.genome).name
            self.output_prefix = genome_name

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 检查至少有一个数据输入|Check at least one data input
        has_rnaseq = self.rnaseq or self.short_reads or self.long_reads
        if not has_rnaseq and not self.transcripts:
            errors.append("必须提供RNA-seq数据(-r/--short-reads/--long-reads)或转录本数据(-e)|Must provide RNA-seq data (-r/--short-reads/--long-reads) or transcripts (-e)")

        # 检查可选文件是否存在|Check optional files if provided
        if self.rnaseq and not os.path.exists(self.rnaseq):
            errors.append(f"RNA-seq文件不存在|RNA-seq file not found: {self.rnaseq}")

        if self.short_reads and not os.path.exists(self.short_reads):
            errors.append(f"二代转录组文件不存在|Short-reads file not found: {self.short_reads}")

        if self.long_reads and not os.path.exists(self.long_reads):
            errors.append(f"三代转录组文件不存在|Long-reads file not found: {self.long_reads}")

        if self.transcripts and not os.path.exists(self.transcripts):
            errors.append(f"转录本文件不存在|Transcripts file not found: {self.transcripts}")

        if self.proteins and not os.path.exists(self.proteins):
            errors.append(f"蛋白质文件不存在|Proteins file not found: {self.proteins}")

        if self.uniprot and not os.path.exists(self.uniprot):
            errors.append(f"UniProt文件不存在|UniProt file not found: {self.uniprot}")

        if self.cds_gff and not os.path.exists(self.cds_gff):
            errors.append(f"CDS GFF文件不存在|CDS GFF file not found: {self.cds_gff}")

        if self.mito_contigs and not os.path.exists(self.mito_contigs):
            errors.append(f"线粒体contig文件不存在|Mito contigs file not found: {self.mito_contigs}")

        if self.extra_gff and not os.path.exists(self.extra_gff):
            errors.append(f"额外GFF文件不存在|Extra GFF file not found: {self.extra_gff}")

        # 检查EviAnn路径|Check EviAnn path
        eviann_sh = os.path.join(self.eviann_path, 'bin', 'eviann.sh')
        if not os.path.exists(eviann_sh):
            errors.append(f"EviAnn未找到|EviAnn not found at: {eviann_sh}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
