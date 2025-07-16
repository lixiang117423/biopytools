"""
PLINK GWAS分析配置管理模块 | PLINK GWAS Analysis Configuration Management Module
"""

from dataclasses import dataclass
from pathlib import Path

@dataclass
class PlinkGWASConfig:
    """PLINK GWAS分析配置类 | PLINK GWAS Analysis Configuration Class"""
    
    # 输入文件 | Input files
    vcf_file: str
    phenotype_file: str
    output_dir: str = "gwas_results"
    
    # 质量控制参数 | Quality control parameters
    mind: float = 0.05  # 个体缺失率阈值 | Individual missing rate threshold
    geno: float = 0.05  # SNP缺失率阈值 | SNP missing rate threshold
    maf: float = 0.01   # 最小等位基因频率 | Minor allele frequency
    hwe: float = 1e-6   # Hardy-Weinberg平衡检验P值阈值 | HWE p-value threshold
    
    # LD剪枝参数 | LD pruning parameters
    ld_window_size: int = 50     # LD剪枝窗口大小(kb) | LD window size (kb)
    ld_step_size: int = 5        # LD剪枝步长(SNP数) | LD step size (SNPs)
    ld_r2_threshold: float = 0.2 # LD剪枝r²阈值 | LD r² threshold
    
    # 主成分参数 | PCA parameters
    pca_components: int = 10  # 计算的主成分数量 | Number of PCA components to compute
    pca_use: int = 5          # 关联分析中使用的主成分数量 | Number of PCs to use in association
    
    # 显著性阈值 | Significance thresholds
    significance_threshold: str = "5e-8"  # 全基因组显著性阈值 | Genome-wide significance
    suggestive_threshold: str = "1e-5"    # 提示性关联阈值 | Suggestive threshold
    
    # 计算参数 | Computing parameters
    threads: int = 1  # 使用的线程数 | Number of threads
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        if not Path(self.vcf_file).exists():
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        if not Path(self.phenotype_file).exists():
            errors.append(f"表型文件不存在 | Phenotype file does not exist: {self.phenotype_file}")
        
        if self.mind < 0 or self.mind > 1:
            errors.append(f"个体缺失率阈值必须在0-1之间 | Individual missing rate must be 0-1: {self.mind}")
        
        if self.geno < 0 or self.geno > 1:
            errors.append(f"SNP缺失率阈值必须在0-1之间 | SNP missing rate must be 0-1: {self.geno}")
        
        if self.maf < 0 or self.maf > 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间 | MAF threshold must be 0-0.5: {self.maf}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
