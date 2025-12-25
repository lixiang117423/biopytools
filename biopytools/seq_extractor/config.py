"""
序列提取工具配置管理模块 🔧 | Sequence Extraction Tool Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal

@dataclass
class ExtractorConfig:
    """序列提取配置类 🎛️ | Sequence Extraction Configuration Class"""
    
    # 输入文件 | Input files
    sequence_file: str
    regions_file: str
    output_file: str
    
    # 序列参数 | Sequence parameters
    sequence_type: Literal["dna", "protein"] = "dna"
    
    # 提取参数 | Extraction parameters
    threads: int = 88
    merge_output: bool = True  # 是否合并输出到一个文件 | Whether to merge output to one file
    include_headers: bool = True  # 是否包含区域信息在序列名中 | Whether to include region info in sequence names
    
    # 序列处理参数 | Sequence processing parameters
    reverse_complement: bool = False  # DNA序列是否需要反向互补 | Whether DNA sequences need reverse complement
    translate_dna: bool = False  # 是否将DNA翻译为蛋白质 | Whether to translate DNA to protein
    
    # 输出参数 | Output parameters
    line_width: int = 80  # FASTA序列每行字符数 | Characters per line in FASTA sequence
    verbose: bool = True  # 详细输出 | Verbose output
    
    # 工具路径 | Tool paths (用于DNA序列) | (for DNA sequences)
    samtools_path: str = 'samtools'
    
    def __post_init__(self):
        """初始化后处理 📋 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.sequence_file = os.path.normpath(os.path.abspath(self.sequence_file))
        self.regions_file = os.path.normpath(os.path.abspath(self.regions_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        # 创建输出目录 | Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 ✅ | Validate configuration parameters"""
        errors = []
        
        # 检查序列文件 | Check sequence file
        if not os.path.exists(self.sequence_file):
            errors.append(f"序列文件不存在 | Sequence file does not exist: {self.sequence_file}")
        
        # 检查区域文件 | Check regions file
        if not os.path.exists(self.regions_file):
            errors.append(f"区域文件不存在 | Regions file does not exist: {self.regions_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if self.line_width <= 0:
            errors.append(f"行宽必须为正整数 | Line width must be positive: {self.line_width}")
        
        # 检查序列类型 | Check sequence type
        if self.sequence_type not in ["dna", "protein"]:
            errors.append(f"序列类型必须是'dna'或'protein' | Sequence type must be 'dna' or 'protein': {self.sequence_type}")
        
        # 检查DNA特有参数 | Check DNA-specific parameters
        if self.sequence_type == "protein" and (self.reverse_complement or self.translate_dna):
            errors.append("蛋白质序列不支持反向互补或翻译 | Protein sequences do not support reverse complement or translation")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
