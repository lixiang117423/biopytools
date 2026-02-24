"""
BUSCO分析配置管理模块|BUSCO Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class BUSCOConfig:
    """BUSCO分析配置类|BUSCO Analysis Configuration Class"""
    
    # 输入文件|Input files
    input_path: str
    lineage: str
    output_dir: str = './busco_output'
    
    # BUSCO参数|BUSCO parameters
    mode: str = 'genome'  # genome, transcriptome, proteins
    threads: int = 88
    sample_suffix: str = '*.fa'
    
    # 高级参数|Advanced parameters
    force: bool = False
    augustus: bool = False
    augustus_parameters: Optional[str] = None
    augustus_species: Optional[str] = None
    auto_lineage: bool = False
    auto_lineage_euk: bool = False
    auto_lineage_prok: bool = False
    contig_break: int = 10
    datasets_version: str = 'odb12'
    download_path: Optional[str] = None
    evalue: float = 1e-3
    limit: int = 3
    long: bool = False
    metaeuk: bool = False
    metaeuk_parameters: Optional[str] = None
    metaeuk_rerun_parameters: Optional[str] = None
    miniprot: bool = False
    skip_bbtools: bool = False
    offline: bool = False
    restart: bool = False
    quiet: bool = False
    scaffold_composition: bool = False
    tar: bool = False
    
    # 输出格式|Output format
    output_format: str = 'txt'
    
    # 工具路径|Tool paths
    busco_path: str = 'busco'
    
    # 内部属性|Internal attributes
    base_name: str = 'busco_analysis'
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.download_path:
            self.download_path = os.path.normpath(os.path.abspath(self.download_path))
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查输入路径|Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")
        
        # 检查模式|Check mode
        valid_modes = ['genome', 'geno', 'transcriptome', 'tran', 'proteins', 'prot']
        if self.mode not in valid_modes:
            errors.append(f"无效的分析模式|Invalid analysis mode: {self.mode}. Valid modes: {valid_modes}")
        
        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"CPU线程数必须为正整数|CPU threads must be positive: {self.threads}")
        
        if self.evalue <= 0:
            errors.append(f"E值必须为正数|E-value must be positive: {self.evalue}")
        
        if self.limit <= 0:
            errors.append(f"候选区域限制必须为正整数|Limit must be positive: {self.limit}")
        
        if self.contig_break < 0:
            errors.append(f"Contig断点数必须为非负整数|Contig break must be non-negative: {self.contig_break}")
        
        # 检查样本后缀格式|Check sample suffix format
        if '*' not in self.sample_suffix:
            errors.append(f"样本后缀必须包含通配符*|Sample suffix must contain wildcard *: {self.sample_suffix}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
