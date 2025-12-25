"""
基因组组装配置管理模块 | Genome Assembly Configuration Management Module 🔧
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List

@dataclass
@dataclass
class AssemblyConfig:
    """基因组组装配置类 | Genome Assembly Configuration Class"""
    
    # ===== 必需参数（无默认值）=====
    input_dir: str  # 必须放在最前面，因为没有默认值
    
    # ===== 可选参数（有默认值）=====
    # 输入输出 | Input/Output
    output_dir: str = './assembly_output'
    project_name: str = 'genome_assembly'
    
    # 数据文件 | Data files (auto-detected from input_dir)
    hifi_reads: Optional[str] = None
    hic_r1: Optional[str] = None
    hic_r2: Optional[str] = None
    ont_reads: Optional[str] = None
    ngs_r1: Optional[str] = None
    ngs_r2: Optional[str] = None
    
    # 组装策略 | Assembly strategy
    hic_strategy: str = "complete_juicer"  # complete_juicer, standard_3ddna, simplified_salsa2
    
    # Hifiasm参数 | Hifiasm parameters
    threads: int = 88
    purge_level: int = 1
    purge_max: int = 80
    similarity_threshold: float = 0.75
    n_haplotypes: int = 2
    
    # Hi-C处理参数 | Hi-C processing parameters
    restriction_enzyme: str = "MboI"  # MboI, DpnII, HindIII
    min_contig_size: int = 15000
    edit_rounds: int = 2
    
    # 物种特异参数 | Species-specific parameters
    genome_size: str = "3g"  # auto, or specific size like 3g, 500m
    telomere_motif: str = "CCCTAA"  # human default
    species_type: str = "diploid"  # diploid, haploid, polyploid
    
    # 质量控制参数 | Quality control parameters
    skip_fastqc: bool = True  # 默认跳过FastQC，节省时间
    min_hifi_coverage: int = 30
    min_hic_coverage: int = 50
    min_mapping_rate: float = 0.7
    busco_lineage: str = "auto"  # auto, mammalia_odb10, etc.
    
    # 工具路径 | Tool paths
    hifiasm_path: str = 'hifiasm'
    bwa_path: str = 'bwa'
    samtools_path: str = 'samtools'
    juicer_path: str = 'juicer.sh'
    pipeline_3ddna: str = '3d-dna/run-asm-pipeline.sh'
    juicer_tools: str = 'juicer_tools.jar'
    salsa2_path: str = 'run_pipeline.py'
    skip_dependency_check: bool = False
    
    # 内部属性 | Internal attributes (必须放在最后，因为使用了field)
    detected_data_types: List[str] = field(default_factory=list)
    assembly_strategy: str = "unknown"
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 确保detected_data_types是列表
        if not isinstance(self.detected_data_types, list):
            self.detected_data_types = []
        
        # 检测输入数据 | Detect input data
        try:
            self._detect_input_files()
            self._determine_assembly_strategy()
        except Exception as e:
            # 如果检测失败，至少确保有基本的数据类型
            print(f"⚠️ 数据检测出现问题: {e}")
            if not self.detected_data_types:
                self.detected_data_types = ['HiFi']  # 至少假设有HiFi数据
            self.assembly_strategy = "hifi_only"
    
    def _detect_input_files(self):
        """检测输入目录中的数据文件 | Detect data files in input directory"""
        input_path = Path(self.input_dir)
        if not input_path.exists():
            raise ValueError(f"输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
        # 重置列表
        self.detected_data_types = []
        
        # 检测HiFi数据 | Detect HiFi data
        hifi_patterns = ['*hifi*.f*q*', '*pacbio*.f*q*', '*ccs*.f*q*']
        for pattern in hifi_patterns:
            hifi_files = list(input_path.glob(pattern))
            if hifi_files:
                self.hifi_reads = str(hifi_files[0])  # 取第一个文件
                self.detected_data_types.append('HiFi')
                break
        
        # 检测Hi-C数据 | Detect Hi-C data
        hic_r1_patterns = ['*hic*R1*.f*q*', '*hic*_1.f*q*', '*R1*.f*q*']
        hic_r2_patterns = ['*hic*R2*.f*q*', '*hic*_2.f*q*', '*R2*.f*q*']
        
        for pattern in hic_r1_patterns:
            r1_files = list(input_path.glob(pattern))
            if r1_files:
                self.hic_r1 = str(r1_files[0])
                break
                
        for pattern in hic_r2_patterns:
            r2_files = list(input_path.glob(pattern))
            if r2_files:
                self.hic_r2 = str(r2_files[0])
                break
        
        if self.hic_r1 and self.hic_r2:
            self.detected_data_types.append('Hi-C')
        
        # 检测ONT数据 | Detect ONT data
        ont_patterns = ['*ont*.f*q*', '*nanopore*.f*q*', '*long*.f*q*']
        for pattern in ont_patterns:
            ont_files = list(input_path.glob(pattern))
            if ont_files:
                self.ont_reads = str(ont_files[0])
                self.detected_data_types.append('ONT')
                break
        
        # 检测NGS数据 | Detect NGS data
        ngs_r1_patterns = ['*ngs*R1*.f*q*', '*illumina*_1.f*q*', '*short*R1*.f*q*']
        ngs_r2_patterns = ['*ngs*R2*.f*q*', '*illumina*_2.f*q*', '*short*R2*.f*q*']
        
        for pattern in ngs_r1_patterns:
            r1_files = list(input_path.glob(pattern))
            if r1_files:
                self.ngs_r1 = str(r1_files[0])
                break
                
        for pattern in ngs_r2_patterns:
            r2_files = list(input_path.glob(pattern))
            if r2_files:
                self.ngs_r2 = str(r2_files[0])
                break
        
        if self.ngs_r1 and self.ngs_r2:
            self.detected_data_types.append('NGS')
        
        # 确保至少检测到HiFi数据
        if 'HiFi' not in self.detected_data_types:
            raise ValueError("未检测到HiFi数据！HiFi数据是必需的 | No HiFi data detected! HiFi data is required")
    
    def _determine_assembly_strategy(self):
        """确定组装策略 | Determine assembly strategy"""
        # 确保detected_data_types不为空
        if not self.detected_data_types:
            self.assembly_strategy = "unknown"
            return
            
        data_combo = tuple(sorted(self.detected_data_types))
        
        strategy_map = {
            ('HiFi',): "hifi_only",
            ('Hi-C', 'HiFi'): "hifi_hic",
            ('HiFi', 'ONT'): "hifi_ont", 
            ('HiFi', 'NGS'): "hifi_ngs",
            ('Hi-C', 'HiFi', 'ONT'): "hifi_hic_ont",
            ('Hi-C', 'HiFi', 'NGS'): "hifi_hic_ngs",
            ('HiFi', 'NGS', 'ONT'): "hifi_ont_ngs",
            ('Hi-C', 'HiFi', 'NGS', 'ONT'): "all_data"
        }
        
        self.assembly_strategy = strategy_map.get(data_combo, "unknown")