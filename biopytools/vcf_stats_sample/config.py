"""
VCF基因型统计配置管理模块 | VCF Genotype Statistics Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class VCFStatsConfig:
    """VCF基因型统计配置类 | VCF Genotype Statistics Configuration Class"""
    
    # 必需文件 | Required files
    vcf_file: str
    output_dir: str = './vcf_stats_output'
    
    # 处理参数 | Processing parameters
    min_depth: int = 0  # 最小深度过滤 | Minimum depth filter
    min_qual: float = 0.0  # 最小质量分数过滤 | Minimum quality score filter
    exclude_missing: bool = False  # 是否排除缺失基因型 | Whether to exclude missing genotypes
    
    # 输出选项 | Output options
    output_detailed: bool = True  # 是否输出详细统计 | Whether to output detailed statistics
    output_summary: bool = True  # 是否输出总结统计 | Whether to output summary statistics
    
    # 内部属性 | Internal attributes
    sample_names: List[str] = None  # 样本名称列表 | Sample names list
    total_snps: int = 0  # 总SNP数 | Total number of SNPs
    
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
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.min_depth < 0:
            errors.append(f"最小深度必须非负 | Minimum depth must be non-negative: {self.min_depth}")
        
        if self.min_qual < 0:
            errors.append(f"最小质量分数必须非负 | Minimum quality score must be non-negative: {self.min_qual}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
