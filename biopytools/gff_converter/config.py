"""
GFF格式转换配置管理模块 | GFF Format Conversion Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class GFFConfig:
    """GFF转换配置类 | GFF Conversion Configuration Class"""
    
    # 输入输出文件 | Input/Output files
    input_file: str
    output_file: str
    
    # 物种信息 | Species information
    species_name: str  # e.g., "OV53"
    species_prefix: str  # e.g., "Ov"
    
    # 编号参数 | Numbering parameters
    start_number: int = 10  # 起始编号 | Starting number
    step: int = 10  # 编号步长 | Step size
    
    # 处理参数 | Processing parameters
    threads: int = 88  # 线程数 | Thread count
    
    # 输出选项 | Output options
    verbose: bool = True  # 详细输出 | Verbose output
    keep_intermediate: bool = False  # 保留中间文件 | Keep intermediate files
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        # 创建输出目录 | Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入GFF文件不存在 | Input GFF file does not exist: {self.input_file}")
        
        # 检查文件扩展名 | Check file extension
        if not self.input_file.lower().endswith(('.gff', '.gff3')):
            errors.append(f"输入文件不是GFF格式 | Input file is not GFF format: {self.input_file}")
        
        # 检查物种信息 | Check species information
        if not self.species_name.strip():
            errors.append("物种名称不能为空 | Species name cannot be empty")
        
        if not self.species_prefix.strip():
            errors.append("物种缩写不能为空 | Species prefix cannot be empty")
        
        # 检查数值参数 | Check numeric parameters
        if self.start_number < 1:
            errors.append(f"起始编号必须大于0 | Start number must be greater than 0: {self.start_number}")
        
        if self.step < 1:
            errors.append(f"编号步长必须大于0 | Step size must be greater than 0: {self.step}")
        
        if self.threads < 1:
            errors.append(f"线程数必须大于0 | Thread count must be greater than 0: {self.threads}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
