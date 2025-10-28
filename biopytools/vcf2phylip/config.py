"""
VCF转换工具配置管理模块 | VCF Converter Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class ConverterConfig:
    """VCF转换配置类 | VCF Conversion Configuration Class"""
    
    # 输入输出文件 | Input/Output files
    input_file: str
    output_dir: str = './converted_output'
    output_prefix: Optional[str] = None
    
    # 转换参数 | Conversion parameters
    min_samples_locus: int = 4  # 位点最少样本数 | Minimum samples per locus
    outgroup: str = ""  # 外群样本名 | Outgroup sample name
    
    # 输出格式控制 | Output format control
    phylip_disable: bool = False  # 禁用PHYLIP输出 | Disable PHYLIP output
    fasta: bool = False  # 启用FASTA输出 | Enable FASTA output
    nexus: bool = False  # 启用NEXUS输出 | Enable NEXUS output
    nexus_binary: bool = False  # 启用二进制NEXUS输出 | Enable binary NEXUS output
    
    # 处理选项 | Processing options
    resolve_IUPAC: bool = False  # 随机解析杂合子基因型 | Randomly resolve heterozygous genotypes
    write_used_sites: bool = False  # 保存使用的位点坐标 | Save coordinates of used sites
    threads: int = 88  # 线程数 | Number of threads
    
    # 内部属性 | Internal attributes
    base_name: str = 'vcf_conversion'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入VCF文件 | Check input VCF file
        if not os.path.exists(self.input_file):
            errors.append(f"🚨 输入VCF文件不存在 | Input VCF file does not exist: {self.input_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.min_samples_locus <= 0:
            errors.append(f"🚨 最小样本数必须为正整数 | Minimum samples must be positive: {self.min_samples_locus}")
        
        if self.threads <= 0:
            errors.append(f"🚨 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        # 检查是否至少启用一种输出格式 | Check if at least one output format is enabled
        if self.phylip_disable and not self.fasta and not self.nexus and not self.nexus_binary:
            errors.append("🚨 必须至少启用一种输出格式 | At least one output format must be enabled")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
