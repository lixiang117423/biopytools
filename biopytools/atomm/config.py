"""ATOMM配置模块|ATOMM Configuration Module"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class ATOMMConfig:
    """ATOMM分析配置类|ATOMM Analysis Configuration Class

    Attributes:
        host_genotype: 宿主基因型文件路径(ATOMM格式)|Host genotype file path (ATOMM format)
        pathogen_genotype: 病原基因型文件路径(ATOMM格式)|Pathogen genotype file path (ATOMM format)
        phenotype: 表型文件路径(ATOMM格式)|Phenotype file path (ATOMM format)
        output_dir: 输出目录|Output directory
        maf_threshold: 最小等位基因频率阈值|Minor allele frequency threshold
        host_snp_range: 宿主SNP测试范围(start, end)|Host SNP test range
        pathogen_snp_range: 病原SNP测试范围(start, end)|Pathogen SNP test range
        interaction_host_range: 交互检验宿主SNP范围|Interaction test host SNP range
        interaction_pathogen_range: 交互检验病原SNP范围|Interaction test pathogen SNP range
        optimize_tol: 优化器容忍度|Optimizer tolerance
        optimize_maxiter: 优化器最大迭代次数|Optimizer max iterations
        host_vcf: 宿主VCF文件路径|Host VCF file path
        pathogen_vcf: 病原VCF文件路径|Pathogen VCF file path
        phenotype_matrix: 交叉感染矩阵文件(行=宿主,列=病原)|Cross-infection matrix file
        encoding: 基因型编码方式 "auto"|"haploid"|"dosage"|Genotype encoding
        convert_maf_threshold: VCF转换时的MAF阈值|MAF threshold for VCF conversion
        missing_value: 表型缺失值标记|Missing value marker in phenotype matrix
    """

    host_genotype: str = ""
    pathogen_genotype: str = ""
    phenotype: str = ""
    output_dir: str = "./output"

    maf_threshold: float = 0.05

    host_snp_range: Optional[tuple] = None
    pathogen_snp_range: Optional[tuple] = None
    interaction_host_range: Optional[tuple] = None
    interaction_pathogen_range: Optional[tuple] = None

    optimize_tol: float = 1e-6
    optimize_maxiter: int = 10000

    host_vcf: str = ""
    pathogen_vcf: str = ""
    phenotype_matrix: str = ""
    encoding: str = "auto"
    convert_maf_threshold: float = 0.05
    missing_value: str = "NA"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    @property
    def has_convert_inputs(self):
        """是否使用VCF/表型矩阵转换模式|Whether VCF/matrix conversion mode is active"""
        return bool(self.host_vcf or self.pathogen_vcf or self.phenotype_matrix)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if self.maf_threshold <= 0 or self.maf_threshold >= 0.5:
            errors.append(f"MAF阈值必须在(0, 0.5)之间|MAF threshold must be in (0, 0.5): {self.maf_threshold}")

        if self.has_convert_inputs:
            if self.host_vcf and not Path(self.host_vcf).exists():
                errors.append(f"宿主VCF文件不存在|Host VCF file not found: {self.host_vcf}")
            if self.pathogen_vcf and not Path(self.pathogen_vcf).exists():
                errors.append(f"病原VCF文件不存在|Pathogen VCF file not found: {self.pathogen_vcf}")
            if self.phenotype_matrix and not Path(self.phenotype_matrix).exists():
                errors.append(f"表型矩阵文件不存在|Phenotype matrix file not found: {self.phenotype_matrix}")
            if self.encoding not in ("auto", "haploid", "dosage"):
                errors.append(f"编码方式必须是auto/haploid/dosage之一|Encoding must be auto/haploid/dosage: {self.encoding}")
        else:
            if not self.host_genotype or not Path(self.host_genotype).exists():
                errors.append(f"宿主基因型文件不存在|Host genotype file not found: {self.host_genotype}")
            if not self.pathogen_genotype or not Path(self.pathogen_genotype).exists():
                errors.append(f"病原基因型文件不存在|Pathogen genotype file not found: {self.pathogen_genotype}")
            if not self.phenotype or not Path(self.phenotype).exists():
                errors.append(f"表型文件不存在|Phenotype file not found: {self.phenotype}")

        if errors:
            raise ValueError("\n".join(errors))
