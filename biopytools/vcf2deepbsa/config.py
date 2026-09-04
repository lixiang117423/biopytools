"""vcf2deepbsa 配置|vcf2deepbsa configuration"""
from dataclasses import dataclass
from pathlib import Path

from ..common.paths import expand_path


@dataclass
class Vcf2DeepBsaConfig:
    """VCF→DeepBSA CSV转换配置|VCF to DeepBSA CSV conversion config"""

    # 必需|Required
    input_vcf: str       # 输入 VCF 文件(支持 .gz)|Input VCF file (.gz supported)
    # 可选|Optional
    output_dir: str = "./vcf2deepbsa_output"  # 输出目录|Output directory
    log_level: str = "INFO"

    def __post_init__(self):
        """展开路径,创建输出目录|Expand paths, create output dir"""
        self.input_vcf = expand_path(self.input_vcf)
        self.output_dir = expand_path(self.output_dir)
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not Path(self.input_vcf).is_file():
            errors.append(f"输入文件不存在|Input file not found: {self.input_vcf}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
