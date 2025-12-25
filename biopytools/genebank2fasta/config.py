"""
GenBank序列提取配置管理模块 🔧 | GenBank Sequence Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class ExtractorConfig:
    """序列提取配置类 ⚙️ | Sequence Extraction Configuration Class"""
    
    # 输入输出路径 | Input/Output paths
    input_dir: str
    output_dir: str = './genbank_output'
    
    # 处理参数 | Processing parameters
    threads: int = 88
    min_protein_length: int = 10  # 最小蛋白质长度 | Minimum protein length
    skip_unknown_genes: bool = True  # 跳过unknown基因 | Skip unknown genes
    
    # 输出选项 | Output options
    create_phylogenetic_matrix: bool = False  # 创建系统发育矩阵 | Create phylogenetic matrix
    separate_by_sample: bool = True  # 按样品分离 | Separate by sample
    separate_by_gene: bool = True  # 按基因分离 | Separate by gene
    
    # 内部属性 | Internal attributes
    gb_dir: str = ""
    cds_dir: str = ""
    pep_dir: str = ""
    
    def __post_init__(self):
        """初始化后处理 📋 | Post-initialization processing"""
        # 创建输出路径 | Create output path
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 设置子目录 | Set subdirectories
        self.gb_dir = self.input_dir
        self.cds_dir = os.path.join(self.output_dir, "cds")
        self.pep_dir = os.path.join(self.output_dir, "pep")
    
    def validate(self):
        """验证配置参数 ✅ | Validate configuration parameters"""
        errors = []
        
        # 检查输入目录 | Check input directory
        if not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if self.min_protein_length <= 0:
            errors.append(f"最小蛋白质长度必须为正整数 | Minimum protein length must be positive: {self.min_protein_length}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
