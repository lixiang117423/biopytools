"""
INDEL PAV分析配置管理模块 | INDEL PAV Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PAVConfig:
    """PAV分析配置类 | PAV Analysis Configuration Class"""
    
    # 输入文件 | Input files
    vcf_file: str
    output_file: str = './indel_pav.txt'
    
    # 分析参数 | Analysis parameters
    threads: int = 88
    min_length: int = 1  # 最小INDEL长度 | Minimum INDEL length
    max_length: Optional[int] = None  # 最大INDEL长度 | Maximum INDEL length
    
    # 过滤参数 | Filtering parameters
    min_quality: float = 20.0  # 最小质量分数 | Minimum quality score
    min_depth: int = 5  # 最小深度 | Minimum depth
    max_missing_rate: float = 0.8  # 最大缺失率 | Maximum missing rate
    
    # 输出参数 | Output parameters
    include_complex: bool = False  # 是否包含复杂变异 | Include complex variants
    compress_output: bool = False  # 是否压缩输出 | Compress output
    
    # 工具路径 | Tool paths
    bcftools_path: str = 'bcftools'
    
    # 内部属性 | Internal attributes
    base_name: str = 'indel_pav'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_file).parent
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if self.min_length < 1:
            errors.append(f"最小INDEL长度必须>=1 | Minimum INDEL length must be >=1: {self.min_length}")
        
        if self.max_length is not None and self.max_length < self.min_length:
            errors.append(f"最大INDEL长度必须>=最小长度 | Maximum INDEL length must be >= minimum length")
        
        if not 0 < self.max_missing_rate <= 1:
            errors.append(f"最大缺失率必须在0-1之间 | Maximum missing rate must be between 0-1: {self.max_missing_rate}")
        
        if self.min_quality < 0:
            errors.append(f"最小质量分数必须>=0 | Minimum quality score must be >=0: {self.min_quality}")
        
        if self.min_depth < 1:
            errors.append(f"最小深度必须>=1 | Minimum depth must be >=1: {self.min_depth}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
