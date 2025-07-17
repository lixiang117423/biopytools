"""
最长转录本提取配置管理模块 | Longest mRNA Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class LongestMRNAConfig:
    """最长转录本提取配置类 | Longest mRNA Extraction Configuration Class"""
    
    # 必需文件 | Required files
    genome_file: str
    gff3_file: str
    output_file: str
    
    # 可选参数 | Optional parameters
    gene_info_file: Optional[str] = None  # 如果不指定，会自动生成 | Auto-generated if not specified
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.gff3_file = os.path.normpath(os.path.abspath(self.gff3_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        # 自动生成基因信息文件路径 | Auto-generate gene info file path
        if self.gene_info_file is None:
            genome_dir = os.path.dirname(self.genome_file)
            genome_basename = os.path.basename(self.genome_file).split('.')[0]
            self.gene_info_file = os.path.join(genome_dir, f'{genome_basename}.gene.info.txt')
        else:
            self.gene_info_file = os.path.normpath(os.path.abspath(self.gene_info_file))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        required_files = [
            ('基因组文件 | Genome file', self.genome_file),
            ('GFF3文件 | GFF3 file', self.gff3_file),
        ]
        
        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在 | does not exist: {file_path}")
        
        # 检查输出目录是否可写 | Check if output directory is writable
        output_dir = os.path.dirname(self.output_file)
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
            except Exception as e:
                errors.append(f"无法创建输出目录 | Cannot create output directory: {output_dir} - {e}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
