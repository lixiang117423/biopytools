"""
Augustus流水线配置管理模块 | Augustus Pipeline Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PipelineConfig:
    """Augustus流水线配置类 | Augustus Pipeline Configuration Class"""
    
    # 必需参数 | Required parameters
    species_name: str
    genome_file: str
    gff_file: str
    
    # 输出设置 | Output settings
    output_dir: str = './augustus_output'
    
    # Augustus设置 | Augustus settings
    augustus_path: str = '/share/org/YZWL/yzwl_lixg/miniforge3/envs/Augustus_v.3.5.0/bin'
    
    # 训练参数 | Training parameters
    train_ratio: float = 0.8
    flank_length: int = 1000
    
    # 内部属性 | Internal attributes
    training_file: str = ""
    train_file: str = ""
    test_file: str = ""
    prediction_file: str = ""
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置 | Validate configuration"""
        # 检查必需文件 | Check required files
        if not os.path.exists(self.genome_file):
            raise FileNotFoundError(f"Genome file not found: {self.genome_file}")
        
        if not os.path.exists(self.gff_file):
            raise FileNotFoundError(f"GFF file not found: {self.gff_file}")
        
        # 检查Augustus路径 | Check Augustus path
        if not os.path.exists(self.augustus_path):
            raise FileNotFoundError(f"Augustus path not found: {self.augustus_path}")
        
        # 检查训练比例 | Check training ratio
        if not 0.1 <= self.train_ratio <= 0.9:
            raise ValueError(f"Training ratio must be between 0.1 and 0.9, got: {self.train_ratio}")
        
        # 检查物种名称 | Check species name
        if not self.species_name or not self.species_name.replace('_', '').replace('-', '').isalnum():
            raise ValueError(f"Invalid species name: {self.species_name}")
