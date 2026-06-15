"""
VCF抽样配置管理模块|VCF Sampling Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class VCFSamplerConfig:
    """VCF抽样配置类|VCF Sampling Configuration Class"""

    # 必需文件|Required files
    input_vcf: str
    output_vcf: str

    # 抽样参数|Sampling parameters
    sample_rate: float = 0.25  # 抽样比例 (0.0-1.0)|Sampling rate (0.0-1.0)

    # 随机种子|Random seed
    random_seed: Optional[int] = 1288  # 随机种子，用于可重复性|Random seed for reproducibility

    # 输出参数|Output parameters
    keep_header: bool = True  # 保留VCF头信息|Keep VCF header

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_vcf = os.path.normpath(os.path.abspath(self.input_vcf))
        self.output_vcf = os.path.normpath(os.path.abspath(self.output_vcf))

        # 创建输出目录|Create output directory
        output_dir = Path(self.output_vcf).parent
        output_dir.mkdir(parents=True, exist_ok=True)

        # 验证抽样比例|Validate sample rate
        if not 0.0 < self.sample_rate <= 1.0:
            raise ValueError(
                f"抽样比例必须在0.0到1.0之间|Sample rate must be between 0.0 and 1.0, "
                f"got: {self.sample_rate}"
            )

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_vcf):
            errors.append(
                f"输入VCF文件不存在|Input VCF file does not exist: {self.input_vcf}"
            )

        # 检查文件是否为gz格式|Check if file is gz format
        if not self.input_vcf.endswith('.gz'):
            errors.append(
                f"输入VCF文件必须是gzip压缩格式 (.vcf.gz)|"
                f"Input VCF file must be gzip compressed (.vcf.gz): {self.input_vcf}"
            )

        # 检查输出文件扩展名|Check output file extension
        if not self.output_vcf.endswith('.vcf.gz'):
            errors.append(
                f"输出VCF文件必须是.vcf.gz格式|"
                f"Output VCF file must be .vcf.gz format: {self.output_vcf}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
