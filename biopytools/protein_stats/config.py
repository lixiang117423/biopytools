"""
Protein Stats配置管理模块|Protein Stats Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ProteinStatsConfig:
    """Protein Stats配置类|Protein Stats Configuration Class"""

    # 必需参数|Required parameters
    protein_fasta: str

    # 输出配置|Output configuration
    output_file: str = 'protein_stats.tsv'
    output_format: str = 'tsv'  # tsv, csv, excel

    # 计算选项|Calculation options
    calculate_length: bool = True
    calculate_mw: bool = True
    calculate_pi: bool = True
    calculate_aa_composition: bool = False
    calculate_instability_index: bool = False
    calculate_gravy: bool = False  # 脂肪指数|Gravy (hydropathy)
    calculate_aromacity: bool = False  # 芳香性|Aromaticity

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.protein_fasta = os.path.normpath(os.path.abspath(self.protein_fasta))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.protein_fasta):
            errors.append(f"蛋白序列文件不存在|Protein FASTA file not found: {self.protein_fasta}")

        # 检查输出格式|Check output format
        valid_formats = ['tsv', 'csv', 'excel']
        if self.output_format not in valid_formats:
            errors.append(f"输出格式必须是以下之一|Output format must be one of: {', '.join(valid_formats)}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
