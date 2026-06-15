"""
GWAS2Gene配置管理模块|GWAS2Gene Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def expand_path(path):
    """展开路径|Expand path"""
    return os.path.expandvars(os.path.expanduser(path))


@dataclass
class GWAS2GeneConfig:
    """GWAS2Gene配置类|GWAS2Gene Configuration Class"""

    # 必需参数|Required parameters
    gwas_file: str
    pval_col: str
    gff_file: str
    output_file: str

    # 可选参数|Optional parameters
    threshold: float = 1e-5
    window: int = 200000
    func_file: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.gwas_file = expand_path(self.gwas_file)
        self.gff_file = expand_path(self.gff_file)
        self.output_file = expand_path(self.output_file)

        if self.func_file:
            self.func_file = expand_path(self.func_file)

        # 确保输出目录存在|Ensure output directory exists
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查GWAS文件|Check GWAS file
        if not os.path.exists(self.gwas_file):
            errors.append(f"GWAS文件不存在|GWAS file does not exist: {self.gwas_file}")

        # 检查GFF文件|Check GFF file
        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file does not exist: {self.gff_file}")

        # 检查功能注释文件|Check function annotation file
        if self.func_file and not os.path.exists(self.func_file):
            errors.append(
                f"功能注释文件不存在|Function annotation file does not exist: "
                f"{self.func_file}"
            )

        # 检查参数范围|Check parameter ranges
        if self.threshold <= 0 or self.threshold >= 1:
            errors.append(
                f"P值阈值必须在0-1之间|P-value threshold must be between 0-1: "
                f"{self.threshold}"
            )

        if self.window < 0:
            errors.append(f"窗口大小必须为非负数|Window size must be non-negative: {self.window}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
