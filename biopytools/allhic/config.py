"""
ALLHiC流水线配置管理模块 ⚙️ | ALLHiC Pipeline Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

@dataclass
class PipelineConfig:
    """流水线配置类 | Pipeline Configuration Class"""
    
    # 必需参数 | Required parameters
    reference: str    # 参考基因组文件 | Reference genome file
    read1: str       # Hi-C读段1文件 | Hi-C read 1 file
    read2: str       # Hi-C读段2文件 | Hi-C read 2 file
    chr_num: int     # 染色体数量 | Number of chromosomes
    
    # 基本参数 | Basic parameters
    enzyme: str = "GATC"                    # 酶切位点 | Restriction enzyme motif
    threads: int = 88                       # 线程数 | Number of threads
    workdir: str = "./allhic_output"        # 工作目录 | Working directory
    
    # MapQ过滤参数 | MapQ filtering parameters
    mapq_step1: int = 1                     # Step 1 MapQ阈值 | Step 1 MapQ threshold
    
    # 绘图参数 | Plotting parameters
    bin_size: str = "500k"                  # 二进制大小 | Bin size
    min_bin_size: str = "50k"               # 最小二进制大小 | Minimum bin size
    
    # 软件路径 | Software paths (支持环境变量覆盖 | Support environment variable override)
    bwa_path: str = field(default_factory=lambda: os.getenv("ALLHIC_BWA_PATH", "/share/org/YZWL/yzwl_lixg/miniforge3/envs/Population_genetics/bin"))
    allhic_software_path: str = field(default_factory=lambda: os.getenv("ALLHIC_SOFTWARE_PATH", "/share/org/YZWL/yzwl_lixg/software/ALLHiC/bin"))
    allhic_bin_path: str = field(default_factory=lambda: os.getenv("ALLHIC_BIN_PATH", "/share/org/YZWL/yzwl_lixg/software/ALLHiC/bin"))  # ALLHiC工具在bin目录
    asmkit_path: str = field(default_factory=lambda: os.getenv("ALLHIC_ASMKIT_PATH", "/share/org/YZWL/yzwl_lixg/software/asmkit/asmkit"))
    assembly_visualizer: str = field(default_factory=lambda: os.getenv("ALLHIC_VISUALIZER_PATH", "/share/org/YZWL/yzwl_lixg/software/3d-dna/visualize/run-assembly-visualizer.sh"))
    
    # 跳过步骤标志 | Skip step flags
    skip_steps: Dict[str, bool] = field(default_factory=lambda: {
        "mapping": False,
        "allele": False,
        "prune": False,
        "partition": False,
        "extract": False,
        "rescue": False,
        "optimize": False,
        "build": False,
        "plot": False,
        "asmkit": False
    })
    
    # 诊断模式 | Diagnostic mode
    diagnose_mode: bool = False
    
    # 内部属性 | Internal attributes
    work_dir_abs: str = field(init=False)
    directories: Dict[str, str] = field(default_factory=dict)
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.work_dir_abs = os.path.normpath(os.path.abspath(self.workdir))
        
        # 定义目录结构 | Define directory structure
        self.directories = {
            "data": os.path.join(self.work_dir_abs, "00_data"),
            "mapping": os.path.join(self.work_dir_abs, "01_mapping"),
            "allele": os.path.join(self.work_dir_abs, "02_allele_table"),
            "prune": os.path.join(self.work_dir_abs, "03_pruning"),
            "partition": os.path.join(self.work_dir_abs, "04_partition"),
            "extract": os.path.join(self.work_dir_abs, "05_extract_matrix"),
            "rescue": os.path.join(self.work_dir_abs, "06_rescue"),
            "optimize": os.path.join(self.work_dir_abs, "07_optimize"),
            "build": os.path.join(self.work_dir_abs, "08_build"),
            "plot": os.path.join(self.work_dir_abs, "09_plot"),
            "jbat": os.path.join(self.work_dir_abs, "10_jbat_asmkit")
        }
        
        # 创建目录结构 | Create directory structure
        self._create_directories()
        
        # 标准化输入路径 | Normalize input paths
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        self.read1 = os.path.normpath(os.path.abspath(self.read1))
        self.read2 = os.path.normpath(os.path.abspath(self.read2))
    
    def _create_directories(self):
        """创建目录结构 | Create directory structure"""
        for dir_path in self.directories.values():
            Path(dir_path).mkdir(parents=True, exist_ok=True)
        
        # 创建日志目录 | Create log directory
        log_dir = os.path.join(self.work_dir_abs, "logs")
        Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        warnings = []

        # 检查输入文件 | Check input files
        for file_path, file_name in [
            (self.reference, "Reference genome"),
            (self.read1, "Hi-C read 1"),
            (self.read2, "Hi-C read 2")
        ]:
            if not os.path.exists(file_path):
                errors.append(f"❌ {file_name}不存在 | {file_name} not found: {file_path}")

        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Threads must be positive: {self.threads}")

        if self.chr_num <= 0:
            errors.append(f"❌ 染色体数必须为正整数 | Chromosome number must be positive: {self.chr_num}")

        if self.mapq_step1 < 0:
            errors.append(f"❌ MapQ阈值不能为负数 | MapQ threshold cannot be negative: {self.mapq_step1}")

        # 检查软件路径 | Check software paths
        software_paths = [
            (self.bwa_path, "BWA"),
            (self.allhic_software_path, "ALLHiC"),
            (self.allhic_bin_path, "ALLHiC bin"),
            (self.asmkit_path, "asmkit"),
            (self.assembly_visualizer, "Assembly Visualizer")
        ]

        for path, name in software_paths:
            if not os.path.exists(path):
                warnings.append(f"⚠️ {name}路径不存在 | {name} path not found: {path}")

        if errors:
            raise ValueError("\n".join(errors))

        if warnings:
            # 在诊断模式下，只显示警告而不中断
            if not self.diagnose_mode:
                raise ValueError("\n".join(warnings))

        return True
