"""
Resistify Parser配置管理模块|Resistify Parser Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ResistifyParserConfig:
    """Resistify Parser配置类|Resistify Parser Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str

    # 输出配置|Output configuration
    output_prefix: str = 'resistify_results'
    output_dir: str = '.'

    # 序列提取选项|Sequence extraction options
    extract_nlr_sequences: bool = False
    extract_nbarc_sequences: bool = False

    # 筛选选项|Filtering options
    filter_classification: str = None  # 按分类筛选，如"TN", "CNL", "NL"等
    min_length: int = None
    max_length: int = None
    min_lrr_length: int = None

    # 输出格式|Output format
    output_tsv: bool = True
    output_csv: bool = True
    output_excel: bool = True

    # 是否包含motifs详情|Whether to include motifs details
    include_motifs: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 检查必需文件|Check required files
        required_files = ['results.tsv', 'domains.tsv']
        for file in required_files:
            file_path = os.path.join(self.input_dir, file)
            if not os.path.exists(file_path):
                errors.append(f"必需文件不存在|Required file not found: {file}")

        # 检查可选文件|Check optional files
        if self.extract_nlr_sequences:
            nlr_fasta = os.path.join(self.input_dir, 'nlr.fasta')
            if not os.path.exists(nlr_fasta):
                errors.append(f"NLR序列文件不存在|NLR fasta file not found: {nlr_fasta}")

        if self.extract_nbarc_sequences:
            nbarc_fasta = os.path.join(self.input_dir, 'nbarc.fasta')
            if not os.path.exists(nbarc_fasta):
                errors.append(f"NB-ARC序列文件不存在|NB-ARC fasta file not found: {nbarc_fasta}")

        # 检查长度参数|Check length parameters
        if self.min_length is not None and self.min_length <= 0:
            errors.append(f"最小长度必须为正数|Min length must be positive: {self.min_length}")

        if self.max_length is not None and self.max_length <= 0:
            errors.append(f"最大长度必须为正数|Max length must be positive: {self.max_length}")

        if self.min_length and self.max_length and self.min_length > self.max_length:
            errors.append(f"最小长度不能大于最大长度|Min length cannot be greater than max length")

        if errors:
            raise ValueError("\n".join(errors))

        return True
