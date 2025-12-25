"""
⚙️ 转录组预测分析配置管理模块 | Transcriptome Prediction Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class TranscriptomeConfig:
    """🔧 转录组预测分析配置类 | Transcriptome Prediction Analysis Configuration Class"""
    
    # 输入文件 | Input files
    genome_file: str
    rna_seq_files: List[str]
    output_dir: str = './transcriptome_output'
    samples_file: Optional[str] = None  # Trinity样本文件 | Trinity samples file
    
    # 🧵 分析参数 | Analysis parameters
    threads: int = 88  # 默认线程数 | Default thread count
    
    # 🧬 HISAT2比对参数 | HISAT2 alignment parameters
    hisat2_min_intron: int = 20
    hisat2_max_intron: int = 500000
    hisat2_novel_splicesite: bool = True
    hisat2_dta: bool = True  # 针对转录本组装器优化 | Optimize for transcript assemblers
    
    # 🧩 StringTie重构参数 | StringTie reconstruction parameters
    stringtie_min_length: int = 50
    stringtie_min_coverage: float = 0
    stringtie_min_fpkm: float = 1
    stringtie_min_iso: float = 0.01
    stringtie_conservative: bool = False
    
    # 🔗 Trinity组装参数 | Trinity assembly parameters
    trinity_min_contig_length: int = 200
    trinity_max_memory: str = "200G"
    trinity_cpu: int = 88  # Trinity CPU数量 | Trinity CPU count
    trinity_ss_lib_type: Optional[str] = None  # 链特异性 | Strand specificity
    
    # 📍 PASA参数 | PASA parameters
    pasa_max_intron_length: int = 100000
    pasa_min_percent_aligned: int = 90
    pasa_min_avg_per_id: int = 95
    pasa_aligners: str = "gmap,blat"
    pasa_cpu: int = 88  # PASA CPU数量 | PASA CPU count
    
    # 🔍 TransDecoder参数 | TransDecoder parameters
    transdecoder_min_protein_len: int = 100
    transdecoder_genetic_code: str = "universal"
    transdecoder_complete_orfs_only: bool = False
    
    # 🛠️ 工具路径 | Tool paths
    hisat2_path: str = 'hisat2'
    hisat2_build_path: str = 'hisat2-build'
    stringtie_path: str = 'stringtie'
    trinity_path: str = 'Trinity'
    pasa_path: str = 'Launch_PASA_pipeline.pl'
    transdecoder_longorfs_path: str = 'TransDecoder.LongOrfs'
    transdecoder_predict_path: str = 'TransDecoder.Predict'
    samtools_path: str = 'samtools'
    
    # 📝 内部属性 | Internal attributes
    base_name: str = 'transcriptome_prediction'
    
    def __post_init__(self):
        """🚀 初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.rna_seq_files = [os.path.normpath(os.path.abspath(f)) for f in self.rna_seq_files]
        
        if self.samples_file:
            self.samples_file = os.path.normpath(os.path.abspath(self.samples_file))
        
        # 同步线程数设置 | Sync thread settings
        if hasattr(self, 'threads') and self.trinity_cpu == 88:
            self.trinity_cpu = self.threads
        if hasattr(self, 'threads') and self.pasa_cpu == 88:
            self.pasa_cpu = self.threads
    
    def validate(self):
        """✅ 验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genome file
        if not os.path.exists(self.genome_file):
            errors.append(f"❌ 基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查RNA-seq文件 | Check RNA-seq files
        for rna_file in self.rna_seq_files:
            if not os.path.exists(rna_file):
                errors.append(f"❌ RNA-seq文件不存在 | RNA-seq file does not exist: {rna_file}")
        
        # 检查样本文件 | Check samples file
        if self.samples_file and not os.path.exists(self.samples_file):
            errors.append(f"❌ 样本文件不存在 | Samples file does not exist: {self.samples_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if self.stringtie_min_length <= 0:
            errors.append(f"❌ StringTie最小长度必须为正整数 | StringTie min length must be positive: {self.stringtie_min_length}")
        
        if self.trinity_min_contig_length <= 0:
            errors.append(f"❌ Trinity最小contig长度必须为正整数 | Trinity min contig length must be positive: {self.trinity_min_contig_length}")
        
        if self.transdecoder_min_protein_len <= 0:
            errors.append(f"❌ TransDecoder最小蛋白质长度必须为正整数 | TransDecoder min protein length must be positive: {self.transdecoder_min_protein_len}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
@dataclass
class TranscriptomeConfig:
    """🔧 转录组预测分析配置类 | Transcriptome Prediction Analysis Configuration Class"""
    
    # 输入文件 | Input files
    genome_file: str
    rna_seq_files: List[str]
    output_dir: str = './transcriptome_output'
    samples_file: Optional[str] = None  # Trinity样本文件 | Trinity samples file
    
    # 🧵 分析参数 | Analysis parameters
    threads: int = 88  # 默认线程数 | Default thread count
    
    # 🔄 工作流控制参数 | Workflow control parameters
    skip_trinity: bool = False  # 跳过Trinity de novo组装 | Skip Trinity de novo assembly
    
    # ... 保持其他参数不变 ...
    # 🧬 HISAT2比对参数 | HISAT2 alignment parameters
    hisat2_min_intron: int = 20
    hisat2_max_intron: int = 500000
    hisat2_novel_splicesite: bool = True
    hisat2_dta: bool = True  # 针对转录本组装器优化 | Optimize for transcript assemblers
    
    # 🧩 StringTie重构参数 | StringTie reconstruction parameters
    stringtie_min_length: int = 50
    stringtie_min_coverage: float = 0
    stringtie_min_fpkm: float = 1
    stringtie_min_iso: float = 0.01
    stringtie_conservative: bool = False
    
    # 🔗 Trinity组装参数 | Trinity assembly parameters
    trinity_min_contig_length: int = 200
    trinity_max_memory: str = "200G"
    trinity_cpu: int = 88  # Trinity CPU数量 | Trinity CPU count
    trinity_ss_lib_type: Optional[str] = None  # 链特异性 | Strand specificity
    
    # 📍 PASA参数 | PASA parameters
    pasa_max_intron_length: int = 100000
    pasa_min_percent_aligned: int = 90
    pasa_min_avg_per_id: int = 95
    pasa_aligners: str = "gmap,blat"
    pasa_cpu: int = 88  # PASA CPU数量 | PASA CPU count
    
    # 🔍 TransDecoder参数 | TransDecoder parameters
    transdecoder_min_protein_len: int = 100
    transdecoder_genetic_code: str = "universal"
    transdecoder_complete_orfs_only: bool = False
    
    # 🛠️ 工具路径 | Tool paths
    hisat2_path: str = 'hisat2'
    hisat2_build_path: str = 'hisat2-build'
    stringtie_path: str = 'stringtie'
    trinity_path: str = 'Trinity'
    pasa_path: str = 'Launch_PASA_pipeline.pl'
    transdecoder_longorfs_path: str = 'TransDecoder.LongOrfs'
    transdecoder_predict_path: str = 'TransDecoder.Predict'
    samtools_path: str = 'samtools'
    
    # 📝 内部属性 | Internal attributes
    base_name: str = 'transcriptome_prediction'
    
    def __post_init__(self):
        """🚀 初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.rna_seq_files = [os.path.normpath(os.path.abspath(f)) for f in self.rna_seq_files]
        
        if self.samples_file:
            self.samples_file = os.path.normpath(os.path.abspath(self.samples_file))
        
        # 同步线程数设置 | Sync thread settings
        if hasattr(self, 'threads') and self.trinity_cpu == 88:
            self.trinity_cpu = self.threads
        if hasattr(self, 'threads') and self.pasa_cpu == 88:
            self.pasa_cpu = self.threads
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genome file
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查RNA-seq文件 | Check RNA-seq files
        for rna_file in self.rna_seq_files:
            if not os.path.exists(rna_file):
                errors.append(f"RNA-seq文件不存在 | RNA-seq file does not exist: {rna_file}")
        
        # 检查样本文件 | Check samples file
        if self.samples_file and not os.path.exists(self.samples_file):
            errors.append(f"样本文件不存在 | Samples file does not exist: {self.samples_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if self.stringtie_min_length <= 0:
            errors.append(f"StringTie最小长度必须为正整数 | StringTie min length must be positive: {self.stringtie_min_length}")
        
        if self.trinity_min_contig_length <= 0:
            errors.append(f"Trinity最小contig长度必须为正整数 | Trinity min contig length must be positive: {self.trinity_min_contig_length}")
        
        if self.transdecoder_min_protein_len <= 0:
            errors.append(f"TransDecoder最小蛋白质长度必须为正整数 | TransDecoder min protein length must be positive: {self.transdecoder_min_protein_len}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True