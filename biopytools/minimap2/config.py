"""
Minimap2分析配置管理模块 | Minimap2 Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class Minimap2Config:
    """Minimap2分析配置类 | Minimap2 Analysis Configuration Class"""
    
    # 必需文件 | Required files
    target_genome: str
    query_genome: str
    output_dir: str = './minimap2_output'
    
    # Minimap2参数 | Minimap2 parameters
    preset: str = 'asm5'  # asm5, asm10, asm20等
    threads: int = 8
    minimap2_path: str = 'minimap2'  # minimap2可执行文件路径
    
    # 筛选参数 | Filtering parameters
    min_match_length: int = 1000  # number_match最小值
    min_unmapped_length: int = 1000  # 未比对区间最小长度
    tp_type: str = 'P'  # 保留的tp类型：'S'=secondary, 'P'=primary, 'SP'=both
    
    # 工具路径 | Tool paths
    seqkit_path: str = 'seqkit'  # seqkit可执行文件路径
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.target_genome = os.path.normpath(os.path.abspath(self.target_genome))
        self.query_genome = os.path.normpath(os.path.abspath(self.query_genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 设置输出文件名 | Set output file names
        query_basename = os.path.splitext(os.path.basename(self.query_genome))[0]
        self.paf_file = os.path.join(self.output_dir, f"{query_basename}_alignment.paf")
        self.bed_file = os.path.join(self.output_dir, f"{query_basename}_unmapped.bed")
        self.unmapped_fasta = os.path.join(self.output_dir, f"{query_basename}_unmapped.fa")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        required_files = [
            ('目标基因组文件 | Target genome file', self.target_genome),
            ('查询基因组文件 | Query genome file', self.query_genome),
        ]
        
        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在 | does not exist: {file_path}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")
        
        if self.min_match_length <= 0:
            errors.append(f"最小匹配长度必须为正整数 | Minimum match length must be positive: {self.min_match_length}")
        
        if self.min_unmapped_length <= 0:
            errors.append(f"最小未比对长度必须为正整数 | Minimum unmapped length must be positive: {self.min_unmapped_length}")
        
        if self.tp_type.upper() not in ['S', 'P', 'SP']:
            errors.append(f"无效的tp_type参数 | Invalid tp_type parameter: {self.tp_type} (应为 S, P, 或 SP | should be S, P, or SP)")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
