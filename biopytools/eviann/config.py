"""
EviAnn配置模块|EviAnn Configuration Module
"""

import os
from dataclasses import dataclass
from typing import Optional

from biopytools.common.paths import expand_path, get_db_path


@dataclass
class EviAnnConfig:
    """EviAnn配置类|EviAnn Configuration Class"""

    # 必需参数|Required parameters
    genome: str
    output_dir: str

    # 数据输入(三种模式互斥,至少一种或 -e)|Input modes (exclusive, at least one or -e)
    rnaseq_data: Optional[str] = None  # 文件/目录(逗号分隔),自动识别|files/dirs, auto-classified
    sample_sheet: Optional[str] = None  # 样本清单TSV|Sample sheet TSV
    rnaseq: Optional[str] = None  # EviAnn原生-r描述文件透传|Passthrough -r file
    transcripts: Optional[str] = None
    proteins: Optional[str] = None

    # 可选参数|Optional parameters
    uniprot: Optional[str] = None
    threads: int = 12
    max_intron: Optional[int] = None
    ploidy: int = 2
    cds_gff: Optional[str] = None
    lncrna_tpm: float = 1.0
    min_prot: Optional[int] = None
    partial: bool = False
    functional: bool = False
    mito_contigs: Optional[str] = None
    extra_gff: Optional[str] = None
    debug: bool = False
    verbose: bool = False

    # 软件路径|Software path
    eviann_path: str = '~/miniforge3/envs/eviann_v.2.0.5'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # UniProt 默认:显式参数 > 环境变量/配置文件 databases 段 > 默认 ~/database
        # |UniProt default: explicit > env/config databases > ~/database default
        if not self.uniprot:
            self.uniprot = get_db_path(
                'uniprot_sprot', '~/database/uniprot/uniprot_sprot.fasta',
                'UNIPROT_SPROT_PATH')
        # 展开所有~路径|Expand all ~ paths
        for attr in ('genome', 'output_dir', 'sample_sheet', 'rnaseq',
                     'transcripts', 'proteins', 'uniprot', 'cds_gff',
                     'mito_contigs', 'extra_gff', 'eviann_path'):
            val = getattr(self, attr)
            if val:
                setattr(self, attr, expand_path(val))
        # rnaseq_data 逗号分隔逐条展开|expand each comma-separated entry
        if self.rnaseq_data:
            self.rnaseq_data = ','.join(
                expand_path(p.strip())
                for p in self.rnaseq_data.split(',') if p.strip())
        os.makedirs(self.output_dir, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 输入模式互斥且至少一种(或 -e)|Input modes exclusive, at least one (or -e)
        n_input = sum(
            1 for v in (self.rnaseq_data, self.sample_sheet, self.rnaseq)
            if v)
        if n_input > 1:
            errors.append(
                "输入模式互斥|Input modes are mutually exclusive: "
                "--rnaseq-data/--sample-sheet/-r 只能用一个|only one allowed")
        if n_input == 0 and not self.transcripts:
            errors.append(
                "必须提供转录组输入(--rnaseq-data/--sample-sheet/-r)或转录本(-e)"
                "|Must provide RNA-seq input or transcripts (-e)")

        # 检查可选文件|Check optional files
        for label, attr in [("转录本|Transcripts", "transcripts"),
                            ("蛋白质|Proteins", "proteins"),
                            ("UniProt|UniProt", "uniprot"),
                            ("CDS GFF|CDS GFF", "cds_gff"),
                            ("线粒体contig|Mito contigs", "mito_contigs"),
                            ("额外GFF|Extra GFF", "extra_gff")]:
            val = getattr(self, attr)
            if val and not os.path.exists(val):
                errors.append(f"{label}文件不存在|File not found: {val}")

        if self.sample_sheet and not os.path.exists(self.sample_sheet):
            errors.append(
                f"样本清单不存在|Sample sheet not found: {self.sample_sheet}")
        if self.rnaseq and not os.path.exists(self.rnaseq):
            errors.append(
                f"描述文件不存在|Description file not found: {self.rnaseq}")
        if self.rnaseq_data:
            for p in self.rnaseq_data.split(','):
                if not os.path.exists(p):
                    errors.append(f"输入不存在|Input not found: {p}")

        # 检查EviAnn|Check EviAnn
        eviann_sh = os.path.join(self.eviann_path, 'bin', 'eviann.sh')
        if not os.path.exists(eviann_sh):
            errors.append(f"EviAnn未找到|EviAnn not found at: {eviann_sh}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
