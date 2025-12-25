"""
VCF重命名配置类 | VCF Renamer Configuration Class
"""

import os


class VCFRenamerConfig:
    """VCF样品名称重命名配置类 | VCF Sample Name Renamer Configuration Class"""

    def __init__(
        self,
        input_vcf: str,
        output_vcf: str,
        prefix: str = "S",
        keep_mapping: bool = True,
        mapping_file: str = None
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            input_vcf: 输入VCF文件路径 | Input VCF file path
            output_vcf: 输出VCF文件路径 | Output VCF file path
            prefix: 新样品名前缀 | New sample name prefix
            keep_mapping: 是否保留映射文件 | Whether to keep mapping file
            mapping_file: 映射文件路径 | Mapping file path
        """
        self.input_vcf = input_vcf
        self.output_vcf = output_vcf
        self.prefix = prefix
        self.keep_mapping = keep_mapping

        # 设置映射文件路径 | Set mapping file path
        if mapping_file:
            self.mapping_file = mapping_file
        else:
            # 从输出文件名生成映射文件名 | Generate mapping filename from output filename
            base = os.path.basename(output_vcf).replace('.vcf.gz', '').replace('.vcf', '')
            self.mapping_file = f"{base}_name_mapping.txt"

        # 运行时变量 | Runtime variables
        self.sample_count = 0
        self.old_samples = []

    def validate(self):
        """
        验证配置参数 | Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时 | When configuration parameters are invalid
        """
        # 检查输入文件 | Check input file
        if not self.input_vcf:
            raise ValueError("输入VCF文件路径不能为空 | Input VCF path cannot be empty")
        if not os.path.exists(self.input_vcf):
            raise ValueError(f"输入VCF文件不存在 | Input VCF not found: {self.input_vcf}")

        # 检查输出文件 | Check output file
        if not self.output_vcf:
            raise ValueError("输出VCF文件路径不能为空 | Output VCF path cannot be empty")

        # 检查前缀 | Check prefix
        if not self.prefix:
            raise ValueError("样品名前缀不能为空 | Sample name prefix cannot be empty")

        # 检查文件扩展名 | Check file extensions
        if not self.input_vcf.endswith(('.vcf.gz', '.vcf')):
            raise ValueError(f"输入文件必须是VCF格式 | Input must be VCF format: {self.input_vcf}")

        return True

    def __repr__(self):
        """配置的字符串表示 | String representation of configuration"""
        return (
            f"VCFRenamerConfig(\n"
            f"  input_vcf={self.input_vcf!r},\n"
            f"  output_vcf={self.output_vcf!r},\n"
            f"  prefix={self.prefix!r},\n"
            f"  mapping_file={self.mapping_file!r}\n"
            f")"
        )
