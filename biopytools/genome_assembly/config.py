"""
基因组组装配置管理模块 | Genome Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class AssemblyConfig:
    """基因组组装配置类 | Genome Assembly Configuration Class"""
    
    # 必需文件 | Required files
    hifi_reads: str
    output_dir: str
    
    # 可选输入文件 | Optional input files
    ont_reads: Optional[str] = None
    reference_genome: Optional[str] = None
    
    # 组装参数 | Assembly parameters
    genome_size: str = '3.2g'
    ploidy: int = 2
    
    # Verkko参数 | Verkko parameters
    verkko_version: str = '1.4.1'
    verkko_grid: bool = False
    verkko_local_cpus: int = 32
    verkko_local_memory: int = 128
    verkko_screen: bool = True
    
    # hifiasm参数 | hifiasm parameters
    hifiasm_version: str = '0.19.6'
    hifiasm_ultra_long: bool = True
    hifiasm_purge_level: int = 3
    hifiasm_similarity: float = 0.8
    
    # Graphasing参数 | Graphasing parameters
    graphasing_version: str = '0.3.1-alpha'
    graphasing_kmer_size: int = 21
    
    # 质控参数 | Quality control parameters
    run_contamination_screen: bool = True
    fcs_version: str = '0.4.0'
    run_flagger: bool = True
    flagger_version: str = '0.3.3'
    run_merqury: bool = True
    merqury_version: str = '1.0'
    run_inspector: bool = True
    inspector_version: str = '1.2'
    
    # 质量评估参数 | Quality assessment parameters
    run_deepvariant: bool = True
    deepvariant_version: str = '1.6'
    run_compleasm: bool = True
    compleasm_version: str = '0.2.5'
    orthodb_version: str = '10'
    
    # 比对参数 | Alignment parameters
    minimap2_version: str = '2.26'
    mashmap_version: str = '3.1.3'
    
    # 系统参数 | System parameters
    threads: int = 32
    memory: int = 128
    tmp_dir: str = '/tmp'
    keep_temp: bool = False
    
    # 工具路径 | Tool paths
    verkko_path: str = 'verkko'
    hifiasm_path: str = 'hifiasm'
    graphasing_path: str = 'graphasing'
    fcs_path: str = 'fcs.py'
    flagger_path: str = 'flagger'
    merqury_path: str = 'merqury.sh'
    nucfreq_path: str = 'nucfreq'
    inspector_path: str = 'inspector.py'
    deepvariant_path: str = 'run_deepvariant'
    compleasm_path: str = 'compleasm'
    minimap2_path: str = 'minimap2'
    mashmap_path: str = 'mashmap'
    
    # 家系分析参数 | Family trio parameters
    trio_mode: bool = False
    parent1_reads: Optional[str] = None
    parent2_reads: Optional[str] = None
    child_reads: Optional[str] = None
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.hifi_reads = os.path.normpath(os.path.abspath(self.hifi_reads))
        if self.ont_reads:
            self.ont_reads = os.path.normpath(os.path.abspath(self.ont_reads))
        if self.reference_genome:
            self.reference_genome = os.path.normpath(os.path.abspath(self.reference_genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        if not os.path.exists(self.hifi_reads):
            errors.append(f"HiFi reads文件不存在 | HiFi reads file does not exist: {self.hifi_reads}")
        
        if self.ont_reads and not os.path.exists(self.ont_reads):
            errors.append(f"ONT reads文件不存在 | ONT reads file does not exist: {self.ont_reads}")
            
        if self.reference_genome and not os.path.exists(self.reference_genome):
            errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {self.reference_genome}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.memory <= 0:
            errors.append(f"内存大小必须为正数 | Memory size must be positive: {self.memory}")
            
        if self.ploidy not in [1, 2]:
            errors.append(f"倍性必须为1或2 | Ploidy must be 1 or 2: {self.ploidy}")
        
        # 家系模式检查 | Trio mode validation
        if self.trio_mode:
            trio_files = [self.parent1_reads, self.parent2_reads, self.child_reads]
            if not all(trio_files):
                errors.append("家系模式需要提供所有三个样本的reads文件 | Trio mode requires all three sample read files")
            else:
                for i, reads_file in enumerate(trio_files, 1):
                    if not os.path.exists(reads_file):
                        errors.append(f"家系样本{i}的reads文件不存在 | Trio sample {i} reads file does not exist: {reads_file}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
