
# """
# 🔧 BLAST分析配置管理模块 | BLAST Analysis Configuration Management Module
# """

# import os
# from dataclasses import dataclass, field
# from pathlib import Path
# from typing import Optional

# @dataclass
# class BLASTConfig:
#     """🔧 BLAST分析配置类 | BLAST Analysis Configuration Class"""
    
#     # 必需参数 | Required parameters
#     target_file: str  # 目标基因序列文件
#     output_dir: str = './blast_output'  # 输出目录
    
#     # 输入数据参数 | Input data parameters
#     input_path: Optional[str] = None  # 输入文件或目录路径
#     sample_map_file: Optional[str] = None  # 样品映射文件
#     auto_generated_map: bool = False  # 是否为自动生成的映射文件
    
#     # BLAST程序参数 | BLAST program parameters
#     blast_type: str = 'blastn'
#     evalue: float = 1e-5
#     max_target_seqs: int = 10
#     word_size: int = 11
#     threads: int = 88
    
#     # 文件格式参数 | File format parameters
#     input_suffix: str = '*.fa'  # 输入文件后缀模式
#     target_db_type: str = 'nucl'
#     outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
    
#     # 过滤参数 | Filtering parameters
#     min_identity: float = 70.0
#     min_coverage: float = 50.0
#     high_quality_evalue: float = 1e-10
    
#     # 样品信息 | Sample information
#     auto_detect_samples: bool = True
#     sample_name_pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$'
    
#     # 工具路径 | Tool paths
#     makeblastdb_path: str = 'makeblastdb'
#     blastn_path: str = 'blastn'
#     blastp_path: str = 'blastp'
#     blastx_path: str = 'blastx'
#     tblastn_path: str = 'tblastn'
#     tblastx_path: str = 'tblastx'
    
#     # 内部属性 | Internal attributes
#     base_name: str = field(default='blast_analysis', init=False)
#     target_base_name: str = field(init=False)  # 新增: 目标文件的基础名称

#     def __post_init__(self):
#         """🔧 初始化后处理"""
#         # 处理参数优先级
#         if self.input_path and self.sample_map_file:
#             print("同时指定了输入路径和样品映射文件，将优先使用样品映射文件")
#             self.input_path = None
        
#         # 创建输出路径
#         self.output_path = Path(self.output_dir)
#         self.output_path.mkdir(parents=True, exist_ok=True)
        
#         # 标准化路径
#         if self.input_path:
#             self.input_path = os.path.normpath(os.path.abspath(self.input_path))
#         self.target_file = os.path.normpath(os.path.abspath(self.target_file))
#         self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
#         if self.sample_map_file:
#             self.sample_map_file = os.path.normpath(os.path.abspath(self.sample_map_file))
        
#         # 设置基础名称
#         self.base_name = f"{self.blast_type}_analysis"
#         self.target_base_name = Path(self.target_file).stem  # 新增: 从目标文件路径获取基础名称
        
#         # 验证BLAST类型
#         if self.blast_type not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
#             raise ValueError(f"不支持的BLAST类型: {self.blast_type}")
    
#     def get_blast_executable(self) -> str:
#         """获取BLAST可执行文件路径"""
#         blast_paths = {
#             'blastn': self.blastn_path,
#             'blastp': self.blastp_path,
#             'blastx': self.blastx_path,
#             'tblastn': self.tblastn_path,
#             'tblastx': self.tblastx_path
#         }
#         return blast_paths[self.blast_type]
    
#     def get_db_type(self) -> str:
#         """获取数据库类型"""
#         db_types = {
#             'blastn': 'nucl',
#             'blastp': 'prot',
#             'blastx': 'prot',
#             'tblastn': 'nucl',
#             'tblastx': 'nucl'
#         }
#         return db_types.get(self.blast_type, self.target_db_type)
    
#     def validate(self):
#         """验证配置参数"""
#         errors = []
        
#         # 检查输入数据参数
#         if not self.input_path and not self.sample_map_file:
#             errors.append("必须指定输入路径(-i)或样品映射文件(-s)中的一个")
        
#         # 检查输入路径
#         if self.input_path and not os.path.exists(self.input_path):
#             errors.append(f"输入路径不存在: {self.input_path}")
        
#         # 检查目标文件
#         if not os.path.exists(self.target_file):
#             errors.append(f"目标基因文件不存在: {self.target_file}")
        
#         # 检查样品映射文件
#         if self.sample_map_file and not self.auto_generated_map and not os.path.exists(self.sample_map_file):
#             errors.append(f"样品映射文件不存在: {self.sample_map_file}")
        
#         # 检查参数范围
#         if self.threads <= 0:
#             errors.append(f"线程数必须为正整数: {self.threads}")
#         if not 0 < self.evalue <= 1:
#             errors.append(f"E-value阈值必须在0-1之间: {self.evalue}")
#         if not 0 <= self.min_identity <= 100:
#             errors.append(f"最小相似度必须在0-100之间: {self.min_identity}")
#         if not 0 <= self.min_coverage <= 100:
#             errors.append(f"最小覆盖度必须在0-100之间: {self.min_coverage}")
        
#         if errors:
#             raise ValueError("\n".join(errors))
#         return True

"""
🔧 BLAST分析配置管理模块 | BLAST Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

@dataclass
class BLASTConfig:
    """🔧 BLAST分析配置类 | BLAST Analysis Configuration Class"""
    
    # 必需参数 | Required parameters
    target_file: str  # 目标基因序列文件 | Target gene sequence file
    output_dir: str = './blast_output'  # 输出目录 | Output directory
    
    # 输入数据参数 | Input data parameters
    input_path: Optional[str] = None  # 输入文件或目录路径 | Input file or directory path
    sample_map_file: Optional[str] = None  # 样品映射文件 | Sample mapping file
    auto_generated_map: bool = False  # 是否为自动生成的映射文件 | Whether it is an auto-generated mapping file
    
    # BLAST程序参数 | BLAST program parameters
    blast_type: str = 'blastn'
    evalue: float = 1e-5
    max_target_seqs: int = 10
    word_size: int = 11
    threads: int = 88
    
    # 文件格式参数 | File format parameters
    input_suffix: str = '*.fa'  # 输入文件后缀模式 | Input file suffix pattern
    target_db_type: str = 'nucl'
    outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
    
    # 过滤参数 | Filtering parameters
    min_identity: float = 70.0
    min_coverage: float = 50.0
    high_quality_evalue: float = 1e-10
    
    # 样品信息 | Sample information
    auto_detect_samples: bool = True
    sample_name_pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$'
    
    # 工具路径 | Tool paths
    makeblastdb_path: str = 'makeblastdb'
    blastn_path: str = 'blastn'
    blastp_path: str = 'blastp'
    blastx_path: str = 'blastx'
    tblastn_path: str = 'tblastn'
    tblastx_path: str = 'tblastx'
    
    # 内部属性 | Internal attributes
    base_name: str = field(default='blast_analysis', init=False)
    target_base_name: str = field(init=False)  # 新增: 目标文件的基础名称 | New: Base name of the target file

    def __post_init__(self):
        """🔧 初始化后处理 | Post-initialization processing"""
        # 处理参数优先级 | Handle parameter priority
        if self.input_path and self.sample_map_file:
            print("同时指定了输入路径和样品映射文件，将优先使用样品映射文件 | Both input path and sample map file are specified, the sample map file will be used with priority")
            self.input_path = None
        
        # 创建输出路径 | Create output path
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Standardize paths
        if self.input_path:
            self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.target_file = os.path.normpath(os.path.abspath(self.target_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        if self.sample_map_file:
            self.sample_map_file = os.path.normpath(os.path.abspath(self.sample_map_file))
        
        # 设置基础名称 | Set base names
        self.base_name = f"{self.blast_type}_analysis"
        self.target_base_name = Path(self.target_file).stem  # 新增: 从目标文件路径获取基础名称 | New: Get base name from target file path
        
        # 验证BLAST类型 | Validate BLAST type
        if self.blast_type not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
            raise ValueError(f"不支持的BLAST类型 | Unsupported BLAST type: {self.blast_type}")
    
    def get_blast_executable(self) -> str:
        """获取BLAST可执行文件路径 | Get BLAST executable path"""
        blast_paths = {
            'blastn': self.blastn_path,
            'blastp': self.blastp_path,
            'blastx': self.blastx_path,
            'tblastn': self.tblastn_path,
            'tblastx': self.tblastx_path
        }
        return blast_paths[self.blast_type]
    
    def get_db_type(self) -> str:
        """获取数据库类型 | Get database type"""
        db_types = {
            'blastn': 'nucl',
            'blastp': 'prot',
            'blastx': 'prot',
            'tblastn': 'nucl',
            'tblastx': 'nucl'
        }
        return db_types.get(self.blast_type, self.target_db_type)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入数据参数 | Check input data parameters
        if not self.input_path and not self.sample_map_file:
            errors.append("必须指定输入路径(-i)或样品映射文件(-s)中的一个 | Must specify one of input path (-i) or sample map file (-s)")
        
        # 检查输入路径 | Check input path
        if self.input_path and not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在 | Input path does not exist: {self.input_path}")
        
        # 检查目标文件 | Check target file
        if not os.path.exists(self.target_file):
            errors.append(f"目标基因文件不存在 | Target gene file does not exist: {self.target_file}")
        
        # 检查样品映射文件 | Check sample mapping file
        if self.sample_map_file and not self.auto_generated_map and not os.path.exists(self.sample_map_file):
            errors.append(f"样品映射文件不存在 | Sample mapping file does not exist: {self.sample_map_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be a positive integer: {self.threads}")
        if not 0 < self.evalue <= 1:
            errors.append(f"E-value阈值必须在0-1之间 | E-value threshold must be between 0 and 1: {self.evalue}")
        if not 0 <= self.min_identity <= 100:
            errors.append(f"最小相似度必须在0-100之间 | Minimum identity must be between 0 and 100: {self.min_identity}")
        if not 0 <= self.min_coverage <= 100:
            errors.append(f"最小覆盖度必须在0-100之间 | Minimum coverage must be between 0 and 100: {self.min_coverage}")
        
        if errors:
            raise ValueError("\n".join(errors))
        return True