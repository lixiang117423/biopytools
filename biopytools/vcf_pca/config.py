"""
VCF PCA分析配置管理模块|VCF PCA Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PCAConfig:
    """PCA分析配置类|PCA Analysis Configuration Class"""
    
    # 输入文件|Input files
    vcf_file: str
    output_dir: str = './pca_output'
    sample_info_file: Optional[str] = None
    
    # PCA参数|PCA parameters
    components: int = 10
    
    # 质控参数|Quality control parameters
    maf: float = 0.05  # Minor allele frequency threshold
    missing_rate: float = 0.1  # Missing genotype rate threshold
    hwe_pvalue: float = 1e-6  # Hardy-Weinberg equilibrium p-value
    skip_qc: bool = False  # 跳过质控过滤|Skip quality control filtering
    
    # 可视化参数|Visualization parameters
    plot: bool = False
    group_column: Optional[str] = None
    
    # 工具路径|Tool paths
    plink_path: str = 'plink'
    bcftools_path: str = 'bcftools'
    
    # 内部属性|Internal attributes
    base_name: str = 'vcf_pca'
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.sample_info_file:
            self.sample_info_file = os.path.normpath(os.path.abspath(self.sample_info_file))
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件|Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")
        
        # 检查样本信息文件|Check sample info file
        if self.sample_info_file and not os.path.exists(self.sample_info_file):
            errors.append(f"样本信息文件不存在|Sample info file does not exist: {self.sample_info_file}")
        
        # 检查参数范围|Check parameter ranges
        if self.components <= 0:
            errors.append(f"主成分数量必须为正整数|Number of components must be positive: {self.components}")
        
        if not 0 < self.maf < 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间|MAF threshold must be between 0-0.5: {self.maf}")
        
        if not 0 < self.missing_rate < 1:
            errors.append(f"缺失率阈值必须在0-1之间|Missing rate threshold must be between 0-1: {self.missing_rate}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
