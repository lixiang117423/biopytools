"""
Chromosome Rename Configuration
染色体重命名配置类
"""

import os


class ChromosomeRenameConfig:
    """染色体重命名配置类 | Chromosome Rename Configuration Class"""

    def __init__(
        self,
        input_file: str,
        output_file: str,
        chromosome_number: int
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            input_file: 输入FASTA文件路径 | Input FASTA file path
            output_file: 输出FASTA文件路径 | Output FASTA file path
            chromosome_number: 染色体数量 | Number of chromosomes
        """
        self.input_file = input_file
        self.output_file = output_file
        self.chromosome_number = chromosome_number

    def validate(self):
        """
        验证配置参数 | Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时 | When configuration parameters are invalid
        """
        # 检查输入文件是否存在 | Check if input file exists
        if not self.input_file or not os.path.exists(self.input_file):
            raise ValueError(f"输入文件不存在 | Input file not found: {self.input_file}")

        # 检查参数范围 | Check parameter ranges
        if self.chromosome_number <= 0:
            raise ValueError(f"染色体数量必须大于0 | Chromosome number must be > 0: {self.chromosome_number}")

        # 检查输出目录是否存在，如不存在则创建 | Check if output directory exists
        output_dir = os.path.dirname(self.output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        return True

    def __repr__(self):
        """配置的字符串表示 | String representation of configuration"""
        return (
            f"ChromosomeRenameConfig(\n"
            f"  input_file={self.input_file!r},\n"
            f"  output_file={self.output_file!r},\n"
            f"  chromosome_number={self.chromosome_number!r}\n"
            f")"
        )
