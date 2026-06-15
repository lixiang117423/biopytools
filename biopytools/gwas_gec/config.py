"""
GWAS GEC分析配置管理模块|GWAS GEC Analysis Configuration Management Module
功能: 支持PLINK binary和VCF两种格式|Features: Support both PLINK binary and VCF formats
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class GECConfig:
    """GEC分析配置类|GEC Analysis Configuration Class"""

    # 必需参数|Required parameters
    pfile: str  # GWAS P值文件|GWAS P-value file
    reference: str  # 参考文件（VCF或PLINK binary前缀）|Reference file (VCF or PLINK binary prefix)

    # 路径配置|Path configuration
    kggsee_jar: str = '~/software/kmmsee/kggsee.jar'
    output_dir: str = './gec_output'

    # 处理参数|Processing parameters
    threads: int = 12
    memory: str = '100G'
    filter_maf_le: float = 0.05
    p_value_cutoff: float = 0.05
    keep_ref: bool = True

    # 染色体范围|Chromosome range
    chromosomes: Optional[str] = None

    # P值文件列名配置|P-value file column names
    chrom_col: str = 'CHR'
    pos_col: str = 'BP'
    p_col: str = 'P'

    # Alpha水平|Alpha level
    alpha: float = 0.05

    # 染色体格式转换|Chromosome format conversion
    convert_chrom_format: bool = True  # 自动转换染色体格式|Auto-convert chromosome format

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.pfile = os.path.normpath(os.path.abspath(self.pfile))
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        self.kggsee_jar = os.path.normpath(os.path.abspath(expand_path(self.kggsee_jar)))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 设置日志文件路径|Set log file path
        self.log_file = os.path.join(self.output_dir, 'gwas_gec.log')

        # KGGSee只支持VCF格式
        # KGGSee only supports VCF format

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查P值文件|Check P-value file
        if not os.path.exists(self.pfile):
            errors.append(f"P值文件不存在|P-value file does not exist: {self.pfile}")

        # 检查参考VCF文件|Check reference VCF file
        # KGGSee只支持VCF格式
        if not os.path.exists(self.reference):
            errors.append(f"参考VCF文件不存在|Reference VCF not found: {self.reference}")

        # 检查KGGSee JAR|Check KGGSee JAR
        if not os.path.exists(self.kggsee_jar):
            errors.append(f"KGGSee JAR文件不存在|KGGSee JAR not found: {self.kggsee_jar}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive")

        if not (0 <= self.filter_maf_le <= 1):
            errors.append(f"MAF过滤阈值必须在0-1之间|MAF filter must be between 0-1")

        if not (0 <= self.p_value_cutoff <= 1):
            errors.append(f"P值阈值必须在0-1之间|P-value cutoff must be between 0-1")

        if not (0 <= self.alpha <= 1):
            errors.append(f"Alpha水平必须在0-1之间|Alpha must be between 0-1")

        if errors:
            raise ValueError("\n".join(errors))

        return True
