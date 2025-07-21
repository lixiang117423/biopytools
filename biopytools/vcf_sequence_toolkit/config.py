"""
序列提取配置管理模块 | Sequence Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class SequenceConfig:
    """序列提取配置类 | Sequence Extraction Configuration Class"""
    
    # 必需文件 | Required files
    vcf_file: str
    genome_file: str
    chrom: str
    start: int
    end: int
    
    # 输出配置 | Output configuration
    output_dir: str = './sequence_output'
    
    # 处理参数 | Processing parameters
    use_first_allele: bool = True  # 是否只使用第一个等位基因
    include_reference: bool = True  # 是否包含参考序列
    export_format: str = "tab"  # 输出格式: tab, fasta, csv
    
    # 过滤参数 | Filtering parameters
    min_qual: Optional[int] = None  # 最小质量值过滤
    sample_list: Optional[List[str]] = None  # 指定样品列表
    exclude_samples: Optional[List[str]] = None  # 排除样品列表
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 确保坐标正确 | Ensure coordinates are correct
        if self.start >= self.end:
            raise ValueError(f"起始位置必须小于结束位置 | Start position must be less than end position: {self.start} >= {self.end}")
        
        if self.start < 1:
            raise ValueError(f"起始位置必须为正整数 | Start position must be positive: {self.start}")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        required_files = [
            ('VCF文件 | VCF file', self.vcf_file),
            ('基因组文件 | Genome file', self.genome_file),
        ]
        
        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在 | does not exist: {file_path}")
        
        # 检查输出格式 | Check output format
        if self.export_format not in ["tab", "fasta", "csv"]:
            errors.append(f"不支持的输出格式 | Unsupported output format: {self.export_format}")
        
        # 检查坐标范围 | Check coordinate range
        if self.end - self.start > 10000:
            errors.append(f"区间长度过大 | Region too large: {self.end - self.start} bp (建议<10kb | recommend <10kb)")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
