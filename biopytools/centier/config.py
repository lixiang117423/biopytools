"""
CentIER着丝粒鉴定配置管理模块|CentIER Centromere Identification Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class CentIERConfig:
    """CentIER着丝粒鉴定配置类|CentIER Centromere Identification Configuration Class"""

    # 必需文件|Required files
    genome_fasta: str

    # 路径配置|Path configuration
    centier_path: str = field(
        default_factory=lambda: get_tool_path(
            'centier',
            '~/software/CentIER/CentIER-main',
            'CENTIER_PATH'
        )
    )
    output_dir: str = './centier_output'

    # 可选文件|Optional files
    gff_annotation: Optional[str] = None

    # Hi-C数据文件(可选)|Hi-C data files (optional)
    matrix1: Optional[str] = None
    matrix2: Optional[str] = None
    bed1: Optional[str] = None
    bed2: Optional[str] = None

    # 分析参数|Analysis parameters
    threads: int = 12
    kmer_size: int = 21
    center_tolerance: int = 15
    step_len: int = 10000
    mul_cents: bool = False
    mingap: int = 2
    signal_threshold: float = 0.7

    # 步骤控制|Step control
    step: Optional[int] = None  # 1-6, None表示运行全部步骤|None means run all steps

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.centier_path = expand_path(self.centier_path)
        self.genome_fasta = expand_path(self.genome_fasta)
        self.output_dir = expand_path(self.output_dir)

        # 展开可选路径|Expand optional paths
        if self.gff_annotation:
            self.gff_annotation = expand_path(self.gff_annotation)
        if self.matrix1:
            self.matrix1 = expand_path(self.matrix1)
        if self.matrix2:
            self.matrix2 = expand_path(self.matrix2)
        if self.bed1:
            self.bed1 = expand_path(self.bed1)
        if self.bed2:
            self.bed2 = expand_path(self.bed2)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome_fasta):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_fasta}")

        # 检查CentIER路径|Check CentIER path
        if not os.path.exists(self.centier_path):
            errors.append(f"CentIER路径不存在|CentIER path does not exist: {self.centier_path}")

        # 检查可选文件|Check optional files
        if self.gff_annotation and not os.path.exists(self.gff_annotation):
            errors.append(f"GFF注释文件不存在|GFF annotation file not found: {self.gff_annotation}")

        # 检查Hi-C数据完整性|Check Hi-C data completeness
        hic_files = [self.matrix1, self.matrix2, self.bed1, self.bed2]
        if any(hic_files) and not all(hic_files):
            errors.append("Hi-C分析需要所有4个文件(medtrx1, matrix2, bed1, bed2)|"
                         "Hi-C analysis requires all 4 files (matrix1, matrix2, bed1, bed2)")

        if any(hic_files):
            for i, f in enumerate(hic_files):
                if f and not os.path.exists(f):
                    errors.append(f"Hi-C文件不存在|Hi-C file not found: {f}")

        # 检查步骤参数|Check step parameter
        if self.step is not None and self.step not in [1, 2, 3, 4, 5, 6]:
            errors.append(f"无效的步骤编号|Invalid step number: {self.step} (应为1-6|should be 1-6)")

        # 检查参数范围|Check parameter ranges
        if self.kmer_size <= 0:
            errors.append("kmer_size必须为正数|kmer_size must be positive")

        if self.step_len <= 0:
            errors.append("step_len必须为正数|step_len must be positive")

        if not (0 <= self.signal_threshold <= 1):
            errors.append("signal_threshold必须在0-1之间|signal_threshold must be between 0 and 1")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_centier_script_path(self) -> str:
        """获取centIER.py脚本路径|Get centIER.py script path"""
        return os.path.join(self.centier_path, 'centIER.py')

    def get_bin_path(self) -> str:
        """获取bin目录路径|Get bin directory path"""
        return os.path.join(self.centier_path, 'bin')
