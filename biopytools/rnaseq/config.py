# """
# RNA-seq分析配置管理模块 | RNA-seq Analysis Configuration Management Module
# """

# import os
# from dataclasses import dataclass
# from pathlib import Path
# from typing import Optional, List

# @dataclass
# class RNASeqConfig:
#     """RNA-seq分析配置类 | RNA-seq Analysis Configuration Class"""
    
#     # 必需文件 | Required files
#     genome_file: str
#     gtf_file: str
#     input_path: str  # FASTQ文件目录或样本信息文件 | FASTQ directory or sample info file
#     output_dir: str
    
#     # 处理参数 | Processing parameters
#     threads: int = 8
#     fastq_pattern: Optional[str] = None  # FASTQ文件命名模式 | FASTQ naming pattern
#     remove_bam: str = "no"  # 是否删除BAM文件 | Whether to remove BAM files
    
#     # 内部属性 | Internal attributes
#     samples: List[dict] = None  # 样本信息列表 | Sample information list
    
#     def __post_init__(self):
#         """初始化后处理 | Post-initialization processing"""
#         self.output_path = Path(self.output_dir)
#         self.output_path.mkdir(parents=True, exist_ok=True)
        
#         # 标准化路径 | Normalize paths
#         self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
#         self.gtf_file = os.path.normpath(os.path.abspath(self.gtf_file))
#         self.input_path = os.path.normpath(os.path.abspath(self.input_path))
#         self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
#     def validate(self):
#         """验证配置参数 | Validate configuration parameters"""
#         errors = []
        
#         # 检查必需文件 | Check required files
#         required_files = [
#             ('基因组文件 | Genome file', self.genome_file),
#             ('GTF文件 | GTF file', self.gtf_file),
#         ]
        
#         for file_desc, file_path in required_files:
#             if not os.path.exists(file_path):
#                 errors.append(f"{file_desc}不存在 | does not exist: {file_path}")
        
#         # 检查输入路径 | Check input path
#         if not os.path.exists(self.input_path):
#             errors.append(f"输入路径不存在 | Input path does not exist: {self.input_path}")
        
#         # 检查线程数 | Check thread count
#         if self.threads <= 0:
#             errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")
        
#         # 检查remove_bam参数 | Check remove_bam parameter
#         if self.remove_bam.lower() not in ["yes", "y", "no", "n"]:
#             errors.append(f"无效的remove_bam参数 | Invalid remove_bam parameter: {self.remove_bam}")
        
#         if errors:
#             raise ValueError("\n".join(errors))
        
#         return True


"""
RNA-seq分析配置管理模块 | RNA-seq Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class RNASeqConfig:
    """RNA-seq分析配置类 | RNA-seq Analysis Configuration Class"""

    # 必需文件 | Required files
    genome_file: str
    gtf_file: str
    input_path: str  # FASTQ文件目录或样本信息文件 | FASTQ directory or sample info file
    output_dir: str

    # 处理参数|Processing parameters
    threads: int = 8
    fastq_pattern: Optional[str] = None  # FASTQ文件命名模式|FASTQ naming pattern
    remove_bam: str = "no"  # 是否删除BAM文件|Whether to remove BAM files
    sample_timeout: int = 21600  # 单个样本处理超时时间（秒），默认6小时|Sample processing timeout in seconds (default: 6 hours)

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    # 高级选项 | Advanced options
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False

    # 内部属性 | Internal attributes
    samples: List[dict] = None  # 样本信息列表 | Sample information list

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.gtf_file = os.path.normpath(os.path.abspath(self.gtf_file))
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查必需文件 | Check required files
        required_files = [
            ('基因组文件 | Genome file', self.genome_file),
            ('GTF文件 | GTF file', self.gtf_file),
        ]

        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在 | does not exist: {file_path}")

        # 检查输入路径 | Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在 | Input path does not exist: {self.input_path}")

        # 检查线程数 | Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")

        # 检查remove_bam参数 | Check remove_bam parameter
        if self.remove_bam.lower() not in ["yes", "y", "no", "n"]:
            errors.append(f"无效的remove_bam参数 | Invalid remove_bam parameter: {self.remove_bam}")

        if errors:
            raise ValueError("\n".join(errors))

        return True