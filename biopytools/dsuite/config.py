"""
Dsuite Configuration
Dsuite配置类
"""

import os


class DsuiteConfig:
    """Dsuite配置类 | Dsuite Configuration Class"""

    def __init__(
        self,
        vcf_file: str,
        sets_file: str,
        output_dir: str,
        output_prefix: str = "dsuite",
        dsuite_bin: str = "/share/org/YZWL/yzwl_lixg/software/Dsuite/Build/Dsuite",
        min_alleles: int = 2,
        max_alleles: int = 2,
        variant_type: str = "snps",
        bcftools: str = "bcftools"
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            vcf_file: 输入VCF文件路径 | Input VCF file path
            sets_file: SETS分组文件路径 | SETS file path
            output_dir: 输出目录 | Output directory
            output_prefix: 输出文件前缀 | Output file prefix
            dsuite_bin: Dsuite可执行文件路径 | Dsuite binary path
            min_alleles: 最小等位基因数 | Min number of alleles
            max_alleles: 最大等位基因数 | Max number of alleles
            variant_type: 变异类型 (snps/indels/both) | Variant type
            bcftools: bcftools命令路径 | bcftools command path
        """
        self.vcf_file = vcf_file
        self.sets_file = sets_file
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.dsuite_bin = dsuite_bin
        self.min_alleles = min_alleles
        self.max_alleles = max_alleles
        self.variant_type = variant_type
        self.bcftools = bcftools

    def validate(self):
        """
        验证配置参数 | Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时 | When configuration parameters are invalid
        """
        # 检查VCF文件是否存在 | Check if VCF file exists
        if not self.vcf_file or not os.path.exists(self.vcf_file):
            raise ValueError(f"VCF文件不存在 | VCF file not found: {self.vcf_file}")

        # 检查SETS文件是否存在 | Check if SETS file exists
        if not self.sets_file or not os.path.exists(self.sets_file):
            raise ValueError(f"SETS文件不存在 | SETS file not found: {self.sets_file}")

        # 检查Dsuite是否存在 | Check if Dsuite exists
        if not self.dsuite_bin or not os.path.exists(self.dsuite_bin):
            raise ValueError(f"Dsuite不存在 | Dsuite not found: {self.dsuite_bin}")

        # 检查参数范围 | Check parameter ranges
        if self.min_alleles < 1:
            raise ValueError(f"最小等位基因数必须>=1 | Min alleles must be >= 1: {self.min_alleles}")

        if self.max_alleles < self.min_alleles:
            raise ValueError(f"最大等位基因数必须>=最小等位基因数 | Max alleles must be >= min alleles")

        # 检查变异类型 | Check variant type
        valid_types = ['snps', 'indels', 'both', 'none']
        if self.variant_type not in valid_types:
            raise ValueError(f"无效的变异类型 | Invalid variant type: {self.variant_type}")

        # 创建输出目录 | Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        return True

    def __repr__(self):
        """配置的字符串表示 | String representation of configuration"""
        return (
            f"DsuiteConfig(\n"
            f"  vcf_file={self.vcf_file!r},\n"
            f"  sets_file={self.sets_file!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  output_prefix={self.output_prefix!r},\n"
            f"  min_alleles={self.min_alleles!r},\n"
            f"  max_alleles={self.max_alleles!r},\n"
            f"  variant_type={self.variant_type!r}\n"
            f")"
        )
