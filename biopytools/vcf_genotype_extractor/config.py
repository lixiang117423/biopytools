"""
VCF基因型提取配置管理模块 | VCF Genotype Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Union

@dataclass
class VCFConfig:
    """VCF基因型提取配置类 | VCF Genotype Extraction Configuration Class"""
    
    # 必需参数 | Required parameters
    vcf_file: str
    output_prefix: str = "vcf_genotype"
    
    # 可选参数 | Optional parameters
    samples: Union[str, List[str]] = "all"  # "all" 或样本名称列表 | "all" or list of sample names
    biallelic_only: bool = False  # 是否只要双等位位点 | Whether to keep only biallelic sites
    split_by_chromosome: bool = False  # 是否按染色体拆分输出 | Whether to split output by chromosome
    output_type: str = "txt"  # 输出格式：txt, csv, excel | Output format: txt, csv, excel
    output_dir: str = "./"  # 输出目录 | Output directory
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 确保输出目录存在 | Ensure output directory exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        # 处理样本参数 | Handle samples parameter
        if isinstance(self.samples, str) and self.samples != "all":
            self.samples = [s.strip() for s in self.samples.split(",")]
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查输出格式 | Check output format
        if self.output_type not in ["txt", "csv", "excel"]:
            errors.append(f"不支持的输出格式 | Unsupported output format: {self.output_type}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
