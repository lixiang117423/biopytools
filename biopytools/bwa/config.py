# ===== FILE: bwa_align/config.py =====
"""
BWA比对配置管理模块|BWA Alignment Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

@dataclass
class AlignConfig:
    """BWA比对配置类|BWA Alignment Configuration Class"""
    
    # 必需参数|Required parameters
    genome: str
    input_dir: str
    pattern: str
    
    # 输出参数|Output parameters
    output_dir: str = './bwa_output'
    
    # 性能参数|Performance parameters
    threads: int = 88

    # BWA mem算法参数|BWA mem algorithm parameters
    bwa_k: int = 19  # minimum seed length
    bwa_w: int = 100  # band width
    bwa_d: int = 100  # off-diagonal X-dropoff
    bwa_r: float = 1.5  # internal seed factor
    bwa_y: int = 20  # seed occurrence for 3rd round
    bwa_c: int = 500  # skip seeds with more occurrences
    bwa_drop_ratio: float = 0.50  # 原bwa_D: drop chains shorter than this fraction
    bwa_min_chain: int = 0  # 原bwa_W: discard chain if seeded bases shorter than this
    bwa_m: int = 50  # mate rescue rounds
    bwa_s: bool = False  # 原bwa_S: skip mate rescue
    bwa_p: bool = False  # 原bwa_P: skip pairing
    
    # BWA mem打分参数|BWA mem scoring parameters
    bwa_score_match: int = 1  # 原bwa_A: match score
    bwa_score_mismatch: int = 4  # 原bwa_B: mismatch penalty
    bwa_gap_open: str = "6,6"  # 原bwa_O: gap open penalties
    bwa_gap_ext: str = "1,1"  # 原bwa_E: gap extension penalty
    bwa_clip: str = "5,5"  # 原bwa_L: clipping penalty
    bwa_unpaired: int = 17  # 原bwa_U: unpaired read pair penalty
    
    # BWA mem输入输出参数|BWA mem I/O parameters
    bwa_mark_secondary: bool = False  # 原bwa_M: mark shorter split hits as secondary
    bwa_min_score: int = 30  # 原bwa_T: minimum score to output
    bwa_h: str = "5,200"  # output all alignments if >80% max score
    bwa_all_align: bool = False  # 原bwa_a: output all alignments
    bwa_append_comment: bool = False  # 原bwa_C: append FASTA/FASTQ comment
    bwa_ref_header: bool = False  # 原bwa_V: output reference FASTA header in XR tag
    bwa_soft_clip: bool = False  # 原bwa_Y: use soft clipping for supplementary alignments
    
    # 后处理参数|Post-processing parameters
    markdup: bool = False  # 标记重复序列|Mark duplicates
    remove_dup: bool = False  # 移除重复序列|Remove duplicates
    
    # 覆盖度参数|Coverage parameters
    min_base_quality: int = 0  # samtools depth -q
    min_mapping_quality: int = 0  # samtools depth -Q
    max_depth: int = 0  # samtools depth -d (0=no limit)
    
    # 滑窗参数|Window parameters
    window_size: int = 1000000  # 1M window
    step_size: int = 100000  # 100kb step
    
    # 清理参数|Cleanup parameters
    keep_sam: bool = False  # 保留SAM文件|Keep SAM files
    keep_sorted_bam: bool = True  # 保留排序后的BAM|Keep sorted BAM
    
    # 断点续传|Resume option
    resume: bool = False  # 跳过已完成的样品|Skip completed samples
    
    # 工具路径|Tool paths
    bwa_path: str = 'bwa'
    samtools_path: str = 'samtools'
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directories
        self.output_path = Path(self.output_dir)
        self.bam_dir = self.output_path / "bam"
        self.coverage_dir = self.output_path / "coverage"
        self.window_dir = self.output_path / "windows"
        self.stats_dir = self.output_path / "stats"
        self.log_dir = self.output_path / "logs"
        
        for d in [self.bam_dir, self.coverage_dir, self.window_dir, 
                  self.stats_dir, self.log_dir]:
            d.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.genome = os.path.normpath(os.path.abspath(self.genome))
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome}")

        # 检查输入目录|Check input directory
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Threads must be positive: {self.threads}")

        if self.window_size <= 0 or self.step_size <= 0:
            errors.append(f"窗口大小和步长必须为正数|Window size and step must be positive")

        if self.step_size > self.window_size:
            errors.append(f"步长大于窗口大小，将产生不连续的窗口|Step > window size will create gaps")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def get_bwa_cmd_options(self) -> str:
        """生成BWA mem命令选项|Generate BWA mem command options"""
        options = []
        
        # 算法参数|Algorithm parameters
        options.append(f"-t {self.threads}")
        options.append(f"-k {self.bwa_k}")
        options.append(f"-w {self.bwa_w}")
        options.append(f"-d {self.bwa_d}")
        options.append(f"-r {self.bwa_r}")
        options.append(f"-y {self.bwa_y}")
        options.append(f"-c {self.bwa_c}")
        options.append(f"-D {self.bwa_drop_ratio}")  # 用新参数名，但BWA命令用-D
        options.append(f"-W {self.bwa_min_chain}")   # 用新参数名，但BWA命令用-W
        options.append(f"-m {self.bwa_m}")
        
        if self.bwa_s:  # 小写
            options.append("-S")
        if self.bwa_p:  # 小写
            options.append("-P")
        
        # 打分参数|Scoring parameters
        options.append(f"-A {self.bwa_score_match}")
        options.append(f"-B {self.bwa_score_mismatch}")
        options.append(f"-O {self.bwa_gap_open}")
        options.append(f"-E {self.bwa_gap_ext}")
        options.append(f"-L {self.bwa_clip}")
        options.append(f"-U {self.bwa_unpaired}")
        
        # 输入输出参数|I/O parameters
        if self.bwa_mark_secondary:
            options.append("-M")
        options.append(f"-T {self.bwa_min_score}")
        options.append(f"-h {self.bwa_h}")
        if self.bwa_all_align:
            options.append("-a")
        if self.bwa_append_comment:
            options.append("-C")
        if self.bwa_ref_header:
            options.append("-V")
        if self.bwa_soft_clip:
            options.append("-Y")
        
        return " ".join(options)

# ===== END FILE =====