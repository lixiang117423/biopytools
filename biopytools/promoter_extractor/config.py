"""
启动子提取器配置管理模块|Promoter Extractor Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class PromoterExtractorConfig:
    """启动子提取器配置类|Promoter Extractor Configuration Class"""

    # 必需文件|Required files
    gff_file: str
    genome_file: str
    output_prefix: str = "promoters"

    # 启动子参数|Promoter parameters
    promoter_length: int = 2000
    min_length: int = 0  # 最小接受长度|Minimum acceptable length

    # 基因选择|Gene selection
    gene_list: Optional[str] = None  # 基因ID列表文件|Gene ID list file

    # 输出选项|Output options
    output_bed: bool = True
    output_stats: bool = True

    # 质量控制|Quality control
    allow_partial: bool = True  # 允许边界截断|Allow boundary truncation

    # 处理选项|Processing options
    threads: int = 1
    keep_intermediate: bool = False

    # 日志配置|Logging configuration
    log_level: str = "INFO"
    quiet: bool = False
    verbose: int = 0

    # 执行控制|Execution control
    force: bool = False
    dry_run: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_prefix).parent
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.gff_file = os.path.abspath(self.gff_file)
        self.genome_file = os.path.abspath(self.genome_file)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file does not exist: {self.gff_file}")

        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_file}")

        # 检查基因列表文件|Check gene list file
        if self.gene_list and not os.path.exists(self.gene_list):
            errors.append(f"基因列表文件不存在|Gene list file does not exist: {self.gene_list}")

        # 检查启动子长度|Check promoter length
        if self.promoter_length <= 0:
            errors.append(f"启动子长度必须为正数|Promoter length must be positive: {self.promoter_length}")

        if self.min_length < 0:
            errors.append(f"最小长度不能为负数|Min length cannot be negative: {self.min_length}")

        if self.min_length > self.promoter_length:
            errors.append(f"最小长度不能大于启动子长度|Min length cannot exceed promoter length: {self.min_length} > {self.promoter_length}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("配置错误|Configuration error:\n" + "\n".join(errors))

        return True
