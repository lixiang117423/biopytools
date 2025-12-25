"""
基因组组装配置管理模块 ⚙️ | Genome Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path

@dataclass
class AssemblyConfig:
    """组装配置类 | Assembly Configuration Class"""
    
    # 必需参数 | Required parameters
    hifi_data: str  # HiFi数据文件 | HiFi data file
    hic_r1: str     # Hi-C R1文件 | Hi-C R1 file
    hic_r2: str     # Hi-C R2文件 | Hi-C R2 file
    
    # 基本参数 | Basic parameters
    prefix: str = "genome_sample"  # 样本前缀 | Sample prefix
    threads: int = 88             # 线程数 | Number of threads
    genome_size: str = "1.45g"    # 基因组大小 | Genome size estimate
    n_hap: int = 2                # 倍性 | Ploidy
    
    # 路径参数 | Path parameters
    base_dir: str = "./assembly_output"  # 基础输出目录 | Base output directory
    
    # 内部属性 | Internal attributes
    work_dir: str = None
    raw_dir: str = None
    fasta_dir: str = None
    log_dir: str = None
    stat_dir: str = None
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.base_dir = os.path.normpath(os.path.abspath(self.base_dir))
        self.work_dir = os.path.join(self.base_dir, self.prefix)
        
        # 定义子目录 | Define subdirectories
        self.raw_dir = os.path.join(self.work_dir, "01.raw_output")
        self.fasta_dir = os.path.join(self.work_dir, "02.fasta")
        self.log_dir = os.path.join(self.work_dir, "03.logs")
        self.stat_dir = os.path.join(self.work_dir, "04.statistics")
        
        # 创建目录结构 | Create directory structure
        for dir_path in [self.work_dir, self.raw_dir, self.fasta_dir, self.log_dir, self.stat_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input files
        if not os.path.exists(self.hifi_data):
            errors.append(f"❌ HiFi文件不存在 | HiFi file not found: {self.hifi_data}")
        
        if not os.path.exists(self.hic_r1):
            errors.append(f"❌ Hi-C R1文件不存在 | Hi-C R1 file not found: {self.hic_r1}")
            
        if not os.path.exists(self.hic_r2):
            errors.append(f"❌ Hi-C R2文件不存在 | Hi-C R2 file not found: {self.hic_r2}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Threads must be positive: {self.threads}")
            
        if self.n_hap <= 0:
            errors.append(f"❌ 倍性必须为正整数 | Ploidy must be positive: {self.n_hap}")
        
        if not self.genome_size.lower().endswith(('g', 'm', 'k')):
            errors.append(f"❌ 基因组大小格式错误 | Genome size format error: {self.genome_size}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
