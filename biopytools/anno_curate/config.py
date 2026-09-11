"""anno_curate配置模块|anno_curate Configuration Module"""

import os
from dataclasses import dataclass
from typing import Optional, Tuple

from biopytools.common.paths import expand_path, get_domain_tool_path

VALID_LOG_LEVELS = ('DEBUG', 'INFO', 'WARNING', 'ERROR')
GXF_SUFFIXES = ('.gff', '.gff3', '.gtf')
FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')


@dataclass
class AnnoCurateConfig:
    """anno_curate配置类|anno_curate Configuration Class"""

    input: str = ''
    output_dir: str = './anno_curate_output'
    genome: Optional[str] = None
    rna_bams: Optional[Tuple[str, ...]] = None
    skip_diagnosis: bool = False
    check_utr: bool = True
    min_cds_nt: int = 90
    max_cds_nt: int = 10500
    utr5_frac_max: float = 0.15
    utr3_frac_max: float = 0.30
    relax: float = 0.5
    cpc2_path: str = ''
    samtools_path: str = ''
    force: bool = False
    log_level: str = 'INFO'

    def __post_init__(self):
        """~展开+工具路径解析+负relax截0|Expand ~, resolve tools, clamp relax"""
        for attr in ('input', 'output_dir', 'genome'):
            if getattr(self, attr):
                setattr(self, attr, expand_path(getattr(self, attr)))
        if self.rna_bams:
            self.rna_bams = tuple(expand_path(p) for p in self.rna_bams)
        if not self.cpc2_path:
            # 入口名以 Task 1 实测为准|Entry name verified in Task 1
            self.cpc2_path = get_domain_tool_path(
                'CPC2.py', '~/miniforge3/envs/annot/bin/CPC2.py', 'CPC2_PATH')
        if not self.samtools_path:
            self.samtools_path = get_domain_tool_path(
                'samtools', '~/miniforge3/envs/align/bin/samtools',
                'SAMTOOLS_PATH')
        self.cpc2_path = expand_path(self.cpc2_path)
        self.samtools_path = expand_path(self.samtools_path)
        # 负 relax 截 0(对齐 GSAman CLI)|Negative relax clamps to 0
        if self.relax < 0:
            self.relax = 0.0
        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def sample(self) -> str:
        """样品名=输入basename剥全部后缀层|Sample = basename, all layers stripped"""
        name = os.path.basename(self.input)
        while True:
            root, ext = os.path.splitext(name)
            if ext.lower() in ('.gz',) + GXF_SUFFIXES:
                name = root
            else:
                break
        return name

    def validate(self):
        """收集全部错误一次抛出|Collect all errors, raise once"""
        errors = []
        if not self.input:
            errors.append("必须提供输入|-i/--input is required")
        elif not os.path.exists(self.input):
            errors.append(f"输入不存在|Input not found: {self.input}")
        elif not self.input.lower().endswith(GXF_SUFFIXES):
            errors.append(f"输入应为 {'/'.join(GXF_SUFFIXES)}|Input must be "
                          f"one of {GXF_SUFFIXES}: {self.input}")
        if self.genome and not os.path.exists(self.genome):
            errors.append(f"基因组不存在|Genome not found: {self.genome}")
        for bam in (self.rna_bams or ()):
            if not os.path.exists(bam):
                errors.append(f"--rna-bam 不存在|not found: {bam}")
        if self.min_cds_nt < 0 or self.min_cds_nt >= self.max_cds_nt:
            errors.append(f"阈值非法|min-cds-nt 必须非负且小于 max-cds-nt "
                          f"({self.min_cds_nt} vs {self.max_cds_nt})")
        for label, v in (('--utr5-frac-max', self.utr5_frac_max),
                         ('--utr3-frac-max', self.utr3_frac_max)):
            if not 0 < v <= 1:
                errors.append(f"{label} 必须在 (0,1]|must be in (0,1]: {v}")
        if self.log_level.upper() not in VALID_LOG_LEVELS:
            errors.append(f"日志级别无效|Invalid log level: {self.log_level}")
        if self.skip_diagnosis and (self.genome or self.rna_bams):
            errors.append("--skip-diagnosis 与 --genome/--rna-bam 互斥"
                          "|--skip-diagnosis conflicts with --genome/--rna-bam")
        if errors:
            raise ValueError("\n".join(errors))
        return True
