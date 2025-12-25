"""
基因序列提取配置管理模块 🧬 | Gene Sequence Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class ExtractionConfig:
    """基因序列提取配置类 🧬 | Gene Sequence Extraction Configuration Class"""
    
    # 输入文件 | Input files
    genome_file: str
    gff_file: str
    output_file: str
    
    # 提取参数 | Extraction parameters
    feature_type: str = 'gene'  # 要提取的特征类型 | Feature type to extract
    min_length: int = 0  # 最小基因长度过滤 | Minimum gene length filter
    threads: int = 88  # 线程数 | Number of threads
    
    # 输出选项 | Output options
    verbose: bool = False  # 详细输出 | Verbose output
    line_width: int = 60  # FASTA行宽度 | FASTA line width
    
    # 内部属性 | Internal attributes
    base_name: str = 'gene_extraction'
    
    def __post_init__(self):
        """初始化后处理 🔧 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        # 创建输出目录 | Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 ✅ | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genome file
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查GFF文件 | Check GFF file
        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在 | GFF file does not exist: {self.gff_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.min_length < 0:
            errors.append(f"最小长度必须为非负整数 | Minimum length must be non-negative: {self.min_length}")
        
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if self.line_width <= 0:
            errors.append(f"行宽度必须为正整数 | Line width must be positive: {self.line_width}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
