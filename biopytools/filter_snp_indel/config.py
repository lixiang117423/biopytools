"""
VCF过滤配置管理模块 | VCF Filtering Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class FilterConfig:
    """过滤配置类 | Filtering Configuration Class"""
    
    # 输入输出文件 | Input/Output files
    vcf_file: str
    output_dir: str = './filtered_vcf'
    
    # 通用参数 | General parameters
    threads: int = 88  # 默认线程数 | Default thread number
    bcftools_path: str = 'bcftools'
    
    # SNP过滤参数 (根据文献常用值设置) | SNP filtering parameters (based on literature)
    snp_qual: float = 30.0      # 最小质量值 | Minimum quality score
    snp_dp: int = 10            # 最小测序深度 | Minimum depth
    snp_mq: float = 40.0        # 最小比对质量 | Minimum mapping quality
    snp_qd: float = 2.0         # 质量/深度比 | Quality by depth
    snp_fs: float = 60.0        # FisherStrand偏倚 | FisherStrand bias
    snp_sor: float = 3.0        # Strand Odds Ratio
    snp_mqrs: float = -12.5     # MappingQualityRankSum
    snp_rprs: float = -8.0      # ReadPosRankSum
    snp_maf: float = 0.05       # [新增] 最小次等位基因频率 | Minimum Minor Allele Frequency
    
    # INDEL过滤参数 (根据文献常用值设置) | INDEL filtering parameters (based on literature)
    indel_qual: float = 30.0    # 最小质量值 | Minimum quality score
    indel_dp: int = 10          # 最小测序深度 | Minimum depth
    indel_mq: float = 40.0      # 最小比对质量 | Minimum mapping quality
    indel_qd: float = 2.0       # 质量/深度比 | Quality by depth
    indel_fs: float = 200.0     # FisherStrand偏倚 | FisherStrand bias
    indel_sor: float = 10.0     # Strand Odds Ratio
    indel_rprs: float = -20.0   # ReadPosRankSum
    
    # 内部属性 | Internal attributes
    base_name: str = 'variation'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"❌ VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if self.snp_qual < 0 or self.indel_qual < 0:
            errors.append(f"❌ 质量值不能为负数 | Quality score cannot be negative")
        
        if self.snp_dp < 0 or self.indel_dp < 0:
            errors.append(f"❌ 测序深度不能为负数 | Depth cannot be negative")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
