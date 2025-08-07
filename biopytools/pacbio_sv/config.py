"""
PacBio HiFi结构变异检测配置模块 | PacBio HiFi Structural Variant Detection Configuration Module
"""

import os
from pathlib import Path
from typing import Optional
from dataclasses import dataclass


@dataclass
class PacBioSVConfig:
    """PacBio HiFi SV检测配置类 | PacBio HiFi SV Detection Configuration Class"""
    
    # 输入文件 | Input files
    bam_file: str
    ref_genome: str
    sample_name: str
    
    # 输出目录 | Output directory
    output_dir: str = "./SV_analysis"
    
    # 分析参数 | Analysis parameters
    threads: int = 16
    min_sv_length: int = 50
    min_support: int = 3
    quality_threshold: float = 20.0
    
    # 大片段SV阈值 | Large SV thresholds
    large_sv_threshold: int = 1000
    very_large_sv_threshold: int = 10000
    
    # SURVIVOR合并参数 | SURVIVOR merge parameters
    survivor_distance: int = 1000
    survivor_min_callers: int = 2
    
    # 工具路径 | Tool paths
    pbsv_path: str = "pbsv"
    sniffles_path: str = "sniffles"
    cutesv_path: str = "cuteSV"
    samtools_path: str = "samtools"
    bcftools_path: str = "bcftools"
    bgzip_path: str = "bgzip"
    tabix_path: str = "tabix"
    survivor_path: str = "SURVIVOR"
    
    def __post_init__(self):
        """配置验证 | Configuration validation"""
        # 验证输入文件 | Validate input files
        if not Path(self.bam_file).exists():
            raise ValueError(f"BAM文件不存在 | BAM file not found: {self.bam_file}")
        
        if not Path(self.ref_genome).exists():
            raise ValueError(f"参考基因组文件不存在 | Reference genome not found: {self.ref_genome}")
        
        # 创建输出目录结构 | Create output directory structure
        self.output_dir = Path(self.output_dir).absolute()
        self.results_dir = self.output_dir / "results"
        self.temp_dir = self.output_dir / "temp"
        self.logs_dir = self.output_dir / "logs"
        
        for dir_path in [self.results_dir, self.temp_dir, self.logs_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # 验证参数范围 | Validate parameter ranges
        if self.threads <= 0:
            raise ValueError("线程数必须大于0 | Thread count must be greater than 0")
        
        if self.min_sv_length <= 0:
            raise ValueError("最小SV长度必须大于0 | Minimum SV length must be greater than 0")
