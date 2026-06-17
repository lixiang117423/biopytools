"""
FASTA序列ID重命名配置管理模块|FASTA ID Renamer Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastaIDRenamerConfig:
    """FASTA序列ID重命名配置类|FASTA ID Renamer Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_file: str

    # 重命名规则配置|Renaming rule configuration
    prefix: str = "Chr"                # 序列前缀|Sequence prefix
    use_zero_padding: bool = True      # 是否使用零填充(如Chr01)|Use zero padding (e.g., Chr01)
    padding_width: int = 2             # 填充宽度|Padding width (01, 001, etc.)

    # 染色体提取选项|Chromosome extraction options
    chr_count: int = 0                 # 染色体数量(0表示不提取)|Chromosome count (0 means don't extract)

    # ID映射文件|ID mapping file
    save_mapping: bool = True          # 是否保存ID映射文件|Save ID mapping file
    mapping_file: Optional[str] = None # ID映射文件路径|ID mapping file path

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 如果没有指定映射文件路径，自动生成|Auto-generate mapping file path if not specified
        if self.save_mapping and self.mapping_file is None:
            output_dir = Path(self.output_file).parent
            base_name = Path(self.output_file).stem
            self.mapping_file = str(output_dir / f"{base_name}_id_mapping.txt")

        # 验证格式选项|Validate format options
        if self.use_zero_padding and self.padding_width < 1:
            raise ValueError(f"填充宽度必须>=1|Padding width must be >= 1: {self.padding_width}")

        # 验证染色体数量|Validate chromosome count
        if self.chr_count < 0:
            raise ValueError(f"染色体数量必须>=0|Chromosome count must be >= 0: {self.chr_count}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        # 检查输入文件格式|Check input file format
        if os.path.exists(self.input_file):
            if not self._is_fasta_file(self.input_file):
                errors.append(f"输入文件不是有效的FASTA格式|Input file is not valid FASTA format: {self.input_file}")

        # 检查输出目录|Check output directory
        output_dir = Path(self.output_file).parent
        if not output_dir.exists():
            errors.append(f"输出目录不存在|Output directory not found: {output_dir}")

        # 检查前缀名称合法性|Check prefix name validity
        if not self.prefix or len(self.prefix.strip()) == 0:
            errors.append(f"前缀不能为空|Prefix cannot be empty")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def _is_fasta_file(self, file_path: str) -> bool:
        """检查是否为FASTA文件|Check if file is FASTA format"""
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except Exception:
            return False
