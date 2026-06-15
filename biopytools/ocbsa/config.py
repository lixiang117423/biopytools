"""
OcBSA - 配置模块|OcBSA - Config Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class OcbsaConfig:
    """OcBSA配置类|OcBSA Configuration Class"""

    # 必需参数|Required parameters
    output_dir: str = "./output"

    # VCF输入|VCF input
    input_vcf: Optional[str] = None
    parent1: int = 0
    parent2: int = 0
    pool1: int = 0
    pool2: int = 0

    # 深度过滤|Depth filter
    parent_min_dep: int = 10
    parent_max_dep: int = 100
    pool_min_dep: int = 20
    pool_max_dep: int = 500

    # 滑窗|Sliding window
    window_size: int = 1000000

    # F1特有参数|F1-specific
    pvalue: float = 99

    # F2特有参数|F2-specific
    method: str = "snpindex"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if self.input_vcf:
            self.input_vcf = os.path.normpath(expand_path(self.input_vcf))
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if self.input_vcf and not os.path.isfile(self.input_vcf):
            errors.append(f"VCF文件不存在|VCF file not found: {self.input_vcf}")

        if self.parent1 < 1 or self.parent2 < 1 or self.pool1 < 1 or self.pool2 < 1:
            errors.append("亲本和混池列号必须大于0|Parent and pool column numbers must be > 0")

        if self.method not in ("snpindex", "ED"):
            errors.append(f"分析方法无效|Invalid method: {self.method} (可选|valid: snpindex, ED)")

        if self.window_size <= 0:
            errors.append(f"滑窗大小必须大于0|Window size must be > 0: {self.window_size}")

        if errors:
            raise ValueError("\n".join(errors))


@dataclass
class BsaFigConfig:
    """BSA绘图配置类|BSA Figure Configuration Class"""

    input_file: str = ""
    output_file: str = "output.png"
    plot_type: str = "ocvalue"
    position: Optional[str] = None
    color: str = "plasma_r"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_file = os.path.normpath(expand_path(self.input_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.isfile(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if self.plot_type not in ("ocvalue", "snpindex", "ed"):
            errors.append(f"图表类型无效|Invalid plot type: {self.plot_type} (可选|valid: ocvalue, snpindex, ed)")

        if not (self.output_file.endswith('.png') or self.output_file.endswith('.pdf')):
            errors.append(f"输出文件必须为.png或.pdf格式|Output must be .png or .pdf: {self.output_file}")

        if errors:
            raise ValueError("\n".join(errors))


@dataclass
class BsaPrimerConfig:
    """BSA引物设计配置类|BSA Primer Design Configuration Class"""

    genome: str = ""
    ocvalue_file: str = ""
    region: str = ""
    output_dir: str = "./output"

    # 引物参数|Primer parameters
    primer_num: int = 10
    flank_length: int = 200
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 24
    product_min: int = 70
    product_max: int = 200
    min_tm: float = 50.0
    max_tm: float = 65.0
    min_gc: float = 35.0
    max_gc: float = 65.0
    tm_diff: float = 0.5

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.genome = os.path.normpath(expand_path(self.genome))
        self.ocvalue_file = os.path.normpath(expand_path(self.ocvalue_file))
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.isfile(self.genome):
            errors.append(f"参考基因组不存在|Reference genome not found: {self.genome}")

        if not os.path.isfile(self.ocvalue_file):
            errors.append(f"OcValue文件不存在|OcValue file not found: {self.ocvalue_file}")

        if not self.region:
            errors.append("目标区间不能为空|Region cannot be empty (format: chr,start,end)")

        if errors:
            raise ValueError("\n".join(errors))
