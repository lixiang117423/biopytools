"""
FASTA ID分割配置管理模块 | FASTA ID Splitting Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

@dataclass
class SplitConfig:
    """FASTA ID分割配置类 | FASTA ID Splitting Configuration Class"""
    
    # 输入输出文件 | Input and output files
    input_file: str
    output_file: str = "output.fasta"
    
    # 分割参数 | Splitting parameters
    position: int = 0  # 提取位置，0表示第一个元素 | Extract position, 0 means first element
    delimiter: str = "auto"  # 分隔符: "auto"(自动), "space"(空格), "tab"(制表符), "both"(两者), 或任意字符如","、"|"等 | Delimiter: "auto", "space", "tab", "both", or any character like ",", "|", etc.
    keep_original: bool = False  # 是否保留原始格式作为备份 | Whether to keep original format as backup
    
    # 处理选项 | Processing options
    skip_empty: bool = True  # 跳过空的序列名称行 | Skip empty sequence name lines
    preserve_comments: bool = False  # 是否保留序列名称行中的注释部分 | Whether to preserve comments in sequence name lines
    
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
            errors.append(f"❌ 输入文件不存在 | Input file does not exist: {self.input_file}")
        
        # 检查文件格式 | Check file format
        if not (self.input_file.lower().endswith('.fasta') or 
                self.input_file.lower().endswith('.fa') or
                self.input_file.lower().endswith('.fas')):
            errors.append(f"⚠️ 输入文件不是FASTA格式 | Input file is not in FASTA format: {self.input_file}")
        
        # 检查位置参数 | Check position parameter
        if self.position < 0:
            errors.append(f"❌ 提取位置必须为非负整数 | Extract position must be non-negative integer: {self.position}")
        
        # 检查分隔符选项 | Check delimiter options
        # 现在支持任意字符，只需要检查是否为空字符串
        if self.delimiter == "":
            errors.append(f"❌ 分隔符不能为空字符串 | Delimiter cannot be empty string")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
