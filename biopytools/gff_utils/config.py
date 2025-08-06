"""
GFF3工具配置管理模块 | GFF3 Utilities Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Set

@dataclass
class GFFConfig:
    """GFF3处理配置类 | GFF3 Processing Configuration Class"""
    
    # 必需参数 | Required parameters
    gff3_file: str
    output_file: str
    
    # 可选参数 | Optional parameters
    transcript_types: Set[str] = None
    gene_type: str = 'gene'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 设置默认转录本类型 | Set default transcript types
        if self.transcript_types is None:
            self.transcript_types = {'mRNA', 'transcript'}
        
        # 标准化路径 | Normalize paths
        self.gff3_file = os.path.normpath(os.path.abspath(self.gff3_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        # 创建输出目录 | Create output directory
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input file
        if not os.path.exists(self.gff3_file):
            errors.append(f"GFF3文件不存在 | GFF3 file does not exist: {self.gff3_file}")
        
        # 检查转录本类型 | Check transcript types
        if not self.transcript_types:
            errors.append("转录本类型不能为空 | Transcript types cannot be empty")
        
        # 检查基因类型 | Check gene type
        if not self.gene_type:
            errors.append("基因类型不能为空 | Gene type cannot be empty")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
