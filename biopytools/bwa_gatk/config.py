"""
⚙️ BWA-GATK流程配置管理模块 | BWA-GATK Pipeline Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

@dataclass
class PipelineConfig:
    """流程配置类 | Pipeline Configuration Class"""
    
    # 📁 输入输出 | Input/Output
    reference: str
    input_path: str
    output_dir: str = './bwa_gatk_output'
    
    # 🧬 分析参数 | Analysis parameters
    threads: int = 88
    ploidy: int = 2
    mem_per_thread: int = 10  # GB
    
    # 🎯 可选功能 | Optional features
    intervals: Optional[str] = None  # BED文件 | BED file
    force_restart: bool = False
    dry_run: bool = False
    verbose: bool = False
    
    # 🔧 工具路径 | Tool paths
    bwa_path: str = 'bwa'
    samtools_path: str = 'samtools'
    gatk_path: str = 'gatk'
    
    # 📊 过滤参数 | Filtering parameters (GATK Best Practices)
    # InDel过滤 | InDel filtering
    indel_qd: float = 2.0
    indel_fs: float = 200.0
    indel_sor: float = 10.0
    indel_inbreeding: float = -0.8
    indel_readpos: float = -20.0
    
    # SNP过滤 | SNP filtering
    snp_qd: float = 2.0
    snp_fs: float = 60.0
    snp_sor: float = 5.5
    snp_mq: float = 40.0
    snp_mqranksum: float = -12.5
    snp_readpos: float = -8.0
    
    # 🗂️ 内部属性 | Internal attributes
    output_path: Path = field(init=False)
    bam_dir: Path = field(init=False)
    gvcf_dir: Path = field(init=False)
    vcf_dir: Path = field(init=False)
    stats_dir: Path = field(init=False)
    logs_dir: Path = field(init=False)
    temp_dir: Path = field(init=False)
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.reference = os.path.abspath(self.reference)
        self.input_path = os.path.abspath(self.input_path)
        self.output_dir = os.path.abspath(self.output_dir)
        
        if self.intervals:
            self.intervals = os.path.abspath(self.intervals)
        
        # 创建输出目录 | Create output directories
        self.output_path = Path(self.output_dir)
        self.bam_dir = self.output_path / 'bam'
        self.gvcf_dir = self.output_path / 'gvcf'
        self.vcf_dir = self.output_path / 'vcf'
        self.stats_dir = self.output_path / 'stats'
        self.logs_dir = self.output_path / 'logs'
        self.temp_dir = self.output_path / 'temp'
        
        # 创建所有必需目录 | Create all required directories
        for dir_path in [self.bam_dir, self.gvcf_dir, self.vcf_dir, 
                         self.stats_dir, self.logs_dir, self.temp_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # 计算总内存 | Calculate total memory
        self.total_memory = self.threads * self.mem_per_thread
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查参考基因组 | Check reference genome
        if not os.path.exists(self.reference):
            errors.append(f"❌ 参考基因组不存在 | Reference genome does not exist: {self.reference}")
        
        # 检查输入路径 | Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"❌ 输入路径不存在 | Input path does not exist: {self.input_path}")
        
        # 检查区间文件 | Check intervals file
        if self.intervals and not os.path.exists(self.intervals):
            errors.append(f"❌ 区间文件不存在 | Intervals file does not exist: {self.intervals}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Threads must be positive: {self.threads}")
        
        if self.ploidy <= 0:
            errors.append(f"❌ 倍性必须为正整数 | Ploidy must be positive: {self.ploidy}")
        
        if self.mem_per_thread <= 0:
            errors.append(f"❌ 每线程内存必须为正数 | Memory per thread must be positive: {self.mem_per_thread}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def get_filter_expression(self, variant_type: str) -> str:
        """获取过滤表达式 | Get filter expression"""
        if variant_type.upper() == 'SNP':
            return (f"QD < {self.snp_qd} || FS > {self.snp_fs} || "
                   f"SOR > {self.snp_sor} || MQ < {self.snp_mq} || "
                   f"MQRankSum < {self.snp_mqranksum} || "
                   f"ReadPosRankSum < {self.snp_readpos}")
        elif variant_type.upper() == 'INDEL':
            return (f"QD < {self.indel_qd} || FS > {self.indel_fs} || "
                   f"SOR > {self.indel_sor} || "
                   f"InbreedingCoeff < {self.indel_inbreeding} || "
                   f"ReadPosRankSum < {self.indel_readpos}")
        else:
            raise ValueError(f"Unknown variant type: {variant_type}")
