"""
GenomeThreader 分析配置管理模块 | GenomeThreader Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class GTHConfig:
    """GenomeThreader 分析配置类 | GenomeThreader Analysis Configuration Class"""
    
    # 必需输入文件 | Required input files
    genomic_file: str
    output_dir: str = './gth_output'
    
    # 序列输入文件（至少需要一个） | Sequence input files (at least one required)
    cdna_file: Optional[str] = None
    protein_file: Optional[str] = None
    est_file: Optional[str] = None
    
    # 物种参数 | Species parameters  
    species: Optional[str] = None  # human, mouse, rat, chicken, drosophila, nematode, etc.
    
    # 基本分析参数 | Basic analysis parameters
    threads: int = 88  # 线程数 | Number of threads
    forward_only: bool = False  # 只分析正链 | Analyze forward strand only
    reverse_only: bool = False  # 只分析负链 | Analyze reverse strand only
    cdna_forward: bool = False  # 只比对cDNA正链 | Align cDNA forward strand only
    
    # 比对质量参数 | Alignment quality parameters
    min_alignment_score: float = 0.65  # 最小比对得分 | Minimum alignment score
    gc_min_coverage: int = 50  # 全局链最小覆盖度 | Global chain minimum coverage
    paralogs: bool = False  # 计算旁系同源基因 | Compute paralogous genes
    
    # 性能优化参数 | Performance optimization parameters
    intron_cutout: bool = False  # 启用内含子切除技术 | Enable intron cutout technique
    auto_intron_cutout: int = 0  # 自动内含子切除矩阵大小(MB) | Auto intron cutout matrix size (MB)
    fastdp: bool = False  # 使用快速DP算法 | Use fast DP algorithm
    
    # 输出参数 | Output parameters
    gff3_output: bool = True  # 输出GFF3格式 | Output GFF3 format
    xml_output: bool = False  # 输出XML格式 | Output XML format
    intermediate: bool = False  # 输出中间结果 | Output intermediate results
    skip_alignment_out: bool = False  # 跳过比对输出 | Skip alignment output
    
    # 高级参数 | Advanced parameters
    bssm_file: Optional[str] = None  # BSSM参数文件路径 | BSSM parameter file path
    score_matrix: str = 'BLOSUM62'  # 氨基酸替换评分矩阵 | Amino acid substitution scoring matrix
    translation_table: int = 1  # 密码子翻译表 | Codon translation table
    
    # 序列处理参数 | Sequence processing parameters
    from_pos: int = 1  # 起始位置 | Starting position
    to_pos: Optional[int] = None  # 结束位置 | Ending position
    width: int = 1000000  # 序列宽度 | Sequence width
    
    # 链接参数 | Chaining parameters
    max_gap_width: int = 30000  # 最大间隔宽度 | Maximum gap width
    
    # 输出限制参数 | Output limitation parameters
    first_alignments: int = 0  # 每个基因组序列最大比对数(0=无限制) | Max alignments per genomic sequence (0=unlimited)
    
    # 工具路径 | Tool paths
    gth_path: str = 'gth'
    gthconsensus_path: str = 'gthconsensus'
    
    # 内部属性 | Internal attributes
    base_name: str = 'genomethreader'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genomic_file = os.path.normpath(os.path.abspath(self.genomic_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.cdna_file:
            self.cdna_file = os.path.normpath(os.path.abspath(self.cdna_file))
        if self.protein_file:
            self.protein_file = os.path.normpath(os.path.abspath(self.protein_file))
        if self.est_file:
            self.est_file = os.path.normpath(os.path.abspath(self.est_file))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genomic file
        if not os.path.exists(self.genomic_file):
            errors.append(f"基因组文件不存在 | Genomic file does not exist: {self.genomic_file}")
        
        # 检查是否至少有一个序列文件 | Check if at least one sequence file is provided
        sequence_files = [self.cdna_file, self.protein_file, self.est_file]
        valid_sequence_files = [f for f in sequence_files if f and os.path.exists(f)]
        
        if not valid_sequence_files:
            errors.append("至少需要提供一个cDNA、蛋白质或EST序列文件 | At least one cDNA, protein, or EST sequence file is required")
        
        # 检查序列文件存在性 | Check sequence file existence
        for file_path, file_type in [(self.cdna_file, "cDNA"), (self.protein_file, "蛋白质"), (self.est_file, "EST")]:
            if file_path and not os.path.exists(file_path):
                errors.append(f"{file_type}文件不存在 | {file_type} file does not exist: {file_path}")
        
        # 检查参数范围 | Check parameter ranges
        if not 0 <= self.min_alignment_score <= 1:
            errors.append(f"最小比对得分必须在0-1之间 | Min alignment score must be between 0-1: {self.min_alignment_score}")
        
        if self.gc_min_coverage < 0 or self.gc_min_coverage > 100:
            errors.append(f"全局链最小覆盖度必须在0-100之间 | GC min coverage must be between 0-100: {self.gc_min_coverage}")
        
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.forward_only and self.reverse_only:
            errors.append("不能同时指定只分析正链和只分析负链 | Cannot specify both forward-only and reverse-only")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
