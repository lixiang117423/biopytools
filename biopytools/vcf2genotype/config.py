"""
VCF基因型提取配置管理模块|VCF Genotype Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Union


@dataclass
class VCFConfig:
    """VCF基因型提取配置类|VCF Genotype Extraction Configuration Class"""

    # 必需输入|Required inputs
    vcf_file: str

    # 路径配置|Path configuration
    output_prefix: str = "vcf_genotype"
    output_dir: str = "./"

    # 样本选择|Sample selection
    samples: Union[str, List[str]] = "all"  # "all" 或样本名称列表|"all" or list of sample names

    # 过滤选项|Filtering options
    biallelic_only: bool = False  # 是否只要双等位位点|Whether to keep only biallelic sites
    min_length: Optional[int] = None  # 最小变异长度|Minimum variant length
    max_length: Optional[int] = None  # 最大变异长度|Maximum variant length

    # 输出配置|Output configuration
    split_by_chromosome: bool = False  # 是否按染色体拆分输出|Whether to split output by chromosome
    output_type: str = "txt"  # 输出格式：txt, csv, excel|Output format: txt, csv, excel

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    # 高级选项|Advanced options
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""

        # 标准化路径|Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 处理样本参数|Handle samples parameter
        if isinstance(self.samples, str) and self.samples != "all":
            self.samples = [s.strip() for s in self.samples.split(",")]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查VCF文件|Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        # 检查输出格式|Check output format
        if self.output_type not in ["txt", "csv", "excel"]:
            errors.append(f"不支持的输出格式|Unsupported output format: {self.output_type}")

        # 检查长度过滤参数|Check length filter parameters
        if self.min_length is not None and self.min_length < 0:
            errors.append(f"最小长度必须为非负数|Minimum length must be non-negative: {self.min_length}")
        if self.max_length is not None and self.max_length < 0:
            errors.append(f"最大长度必须为非负数|Maximum length must be non-negative: {self.max_length}")
        if (self.min_length is not None and self.max_length is not None and
            self.min_length > self.max_length):
            errors.append(f"最小长度不能大于最大长度|Minimum length cannot be greater than maximum length: {self.min_length} > {self.max_length}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
