"""
群体遗传分析配置管理模块 | Population Genetics Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class PopGenConfig:
    """群体遗传分析配置类 | Population Genetics Analysis Configuration Class"""
    
    # 输入文件 | Input files
    vcf_file: str
    output_dir: str = './popgen_output'
    group_file: Optional[str] = None

    # 质控参数 | Quality control parameters
    skip_qc: bool = True  # 新增：跳过质控过滤
    maf: float = 0.01
    missing_rate: float = 0.1
    hwe_pvalue: float = 1e-6
    min_dp: int = 10
    max_dp: int = 100
    
    # # 分析参数选择 | Analysis parameter selection
    # calculate_all: bool = False
    # calculate_fst: bool = True
    # calculate_pi: bool = True
    # calculate_theta_w: bool = True
    # calculate_tajima_d: bool = True
    # calculate_ibd: bool = True
    # calculate_ld: bool = True
    # calculate_ne: bool = False
    # 分析参数选择 | Analysis parameter selection
    calculate_all: bool = False
    calculate_fst: bool = False
    calculate_pi: bool = False
    calculate_theta_w: bool = False
    calculate_tajima_d: bool = False
    calculate_ibd: bool = False
    calculate_ld: bool = False
    calculate_ne: bool = False
    
    # 滑动窗口参数 | Sliding window parameters
    window_sizes: List[int] = None
    window_overlap: float = 0.5
    display_window_size: int = 500000
    
    # 质控参数 | Quality control parameters
    maf: float = 0.05
    missing_rate: float = 0.1
    hwe_pvalue: float = 1e-6
    min_dp: int = 10
    max_dp: int = 100
    
    # 输出格式 | Output format
    output_format: str = 'txt'
    
    # 工具路径 | Tool paths
    vcftools_path: str = 'vcftools'
    plink_path: str = 'plink'
    bcftools_path: str = 'bcftools'
    smcpp_path: str = 'smc++'
    
    # 线程数 | Number of threads
    threads: int = 4
    
    # 内部属性 | Internal attributes
    base_name: str = 'popgen'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        if self.window_sizes is None:
            self.window_sizes = [10000, 100000, 500000]
            
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.group_file:
            self.group_file = os.path.normpath(os.path.abspath(self.group_file))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        if self.group_file and not os.path.exists(self.group_file):
            errors.append(f"分组文件不存在 | Group file does not exist: {self.group_file}")
        
        if not 0 < self.maf < 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间 | MAF threshold must be between 0-0.5: {self.maf}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
