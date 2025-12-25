"""
⚙️ 系统发育树构建配置管理模块 | Phylogenetic Tree Builder Configuration Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PhyloConfig:
    """系统发育树构建配置类 | Phylogenetic Tree Builder Configuration Class"""
    
    # 输入文件 | Input files
    input_file: str
    output_dir: str = './phylo_output'
    
    # 序列参数 | Sequence parameters
    seq_type: Optional[str] = None  # 'protein' or 'nucleotide', auto-detect if None
    
    # 运行参数 | Run parameters  
    threads: int = 88
    
    # MAFFT参数 | MAFFT parameters
    mafft_params: str = '--auto'  # Default MAFFT parameters
    
    # FastTree参数 | FastTree parameters
    fasttree_params: str = ''  # Default FastTree parameters
    
    # 工具路径 | Tool paths
    mafft_path: str = 'mafft'
    fasttree_path: str = 'fasttree'
    
    # 内部属性 | Internal attributes
    base_name: str = 'sequences'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 从输入文件名提取基础名称 | Extract base name from input file
        input_basename = os.path.basename(self.input_file)
        self.base_name = os.path.splitext(input_basename)[0]
        
        # 设置输出文件路径 | Set output file paths
        self.cleaned_file = self.output_path / f"{self.base_name}.cleaned.fa"
        self.mafft_file = self.output_path / f"{self.base_name}.mafft.fa"
        self.tree_file = self.output_path / f"{self.base_name}.fasttree.nwk"
        self.id_mapping_file = self.output_path / f"{self.base_name}.id_mapping.txt"
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"❌ 输入文件不存在 | Input file does not exist: {self.input_file}")
        
        # 检查线程数 | Check threads
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Threads must be positive: {self.threads}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
