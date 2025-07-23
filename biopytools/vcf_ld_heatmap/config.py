"""
VCF LD热图分析配置管理模块 | VCF LD Heatmap Analysis Configuration Management Module
"""

import os
from pathlib import Path
from typing import Optional, List

class LDHeatmapConfig:
    """LD热图分析配置类 | LD Heatmap Analysis Configuration Class"""
    
    def __init__(self, vcf_file: str, output_file: str = 'ld_heatmap.png', 
                 region: Optional[str] = None, maf: float = 0.01, 
                 max_snps: int = 1000, samples: Optional[List[str]] = None,
                 figsize: List[int] = [10, 8], dpi: int = 300, 
                 colormap: str = 'RdYlGn_r', title: Optional[str] = None,
                 save_matrix: Optional[str] = None, ld_threshold: float = 0.0,
                 triangle_only: bool = False, verbose: bool = False):
        
        # 输入输出文件 | Input/Output files
        self.vcf_file = os.path.normpath(os.path.abspath(vcf_file))
        self.output_file = os.path.normpath(os.path.abspath(output_file))
        self.output_dir = Path(self.output_file).parent
        self.base_name = Path(self.output_file).stem
        self.output_path = Path(self.output_file)
        
        # 分析参数 | Analysis parameters
        self.region = region
        self.maf = maf
        self.max_snps = max_snps
        self.samples = samples
        self.ld_threshold = ld_threshold
        
        # 图形参数 | Graphics parameters
        self.figsize = tuple(figsize)
        self.dpi = dpi
        self.colormap = colormap
        self.title = title
        self.triangle_only = triangle_only
        
        # 输出参数 | Output parameters
        self.save_matrix = save_matrix
        if self.save_matrix:
            self.save_matrix = os.path.normpath(os.path.abspath(self.save_matrix))
        
        # 其他参数 | Other parameters
        self.verbose = verbose
        
        # 创建输出目录 | Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 验证配置 | Validate configuration
        self.validate()
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查参数范围 | Check parameter ranges
        if not 0 <= self.maf <= 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间 | MAF threshold must be between 0-0.5: {self.maf}")
        
        if self.max_snps <= 0:
            errors.append(f"最大SNP数量必须为正整数 | Maximum SNPs must be positive: {self.max_snps}")
        
        if not 0 <= self.ld_threshold <= 1:
            errors.append(f"LD阈值必须在0-1之间 | LD threshold must be between 0-1: {self.ld_threshold}")
        
        if self.dpi <= 0:
            errors.append(f"DPI必须为正整数 | DPI must be positive: {self.dpi}")
        
        if len(self.figsize) != 2 or any(x <= 0 for x in self.figsize):
            errors.append(f"图形尺寸必须为两个正数 | Figure size must be two positive numbers: {self.figsize}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
