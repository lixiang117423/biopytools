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
    
#     def __post_init__(self):
#         """🔧 初始化后处理"""
#         # 处理参数优先级
#         if self.input_path and self.sample_map_file:
#             print("⚠️  同时指定了输入路径和样品映射文件，将优先使用样品映射文件")
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
    target_file: str  # 目标基因序列文件
    output_dir: str = './blast_output'  # 输出目录
    
    # 输入数据参数 | Input data parameters
    input_path: Optional[str] = None  # 输入文件或目录路径
    sample_map_file: Optional[str] = None  # 样品映射文件
    auto_generated_map: bool = False  # 是否为自动生成的映射文件
    
    # BLAST程序参数 | BLAST program parameters
    blast_type: str = 'blastn'
    evalue: float = 1e-5
    max_target_seqs: int = 10
    word_size: int = 11
    threads: int = 88
    
    # 文件格式参数 | File format parameters
    input_suffix: str = '*.fa'  # 输入文件后缀模式
    target_db_type: str = 'nucl'
    outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
    
    # 过滤参数 | Filtering parameters
    min_identity: float = 70.0
    min_coverage: float = 50.0
    high_quality_evalue: float = 1e-10
    
    # 样品信息 | Sample information
    auto_detect_samples: bool = True
    sample_name_pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$'
    
    # 比对可视化参数 | Alignment visualization parameters
    alignment_output: str = 'both'  # 'none', 'text', 'html', 'both'
    alignment_width: int = 80  # 每行显示字符数
    alignment_output_dir: str = 'alignments'  # 输出目录名
    alignment_min_identity: float = 0.0  # 最小相似度（0表示不限制）
    alignment_min_coverage: float = 0.0  # 最小覆盖度
    alignment_max_per_sample: int = 100  # 每个样品最多显示N个比对
    html_theme: str = 'modern'  # 'modern', 'classic', 'dark'
    html_interactive: bool = True  # 启用JavaScript交互功能
    
    # 工具路径 | Tool paths
    makeblastdb_path: str = 'makeblastdb'
    blastn_path: str = 'blastn'
    blastp_path: str = 'blastp'
    blastx_path: str = 'blastx'
    tblastn_path: str = 'tblastn'
    tblastx_path: str = 'tblastx'
    
    # 内部属性 | Internal attributes
    base_name: str = field(default='blast_analysis', init=False)
    
    def __post_init__(self):
        """🔧 初始化后处理"""
        # 处理参数优先级
        if self.input_path and self.sample_map_file:
            print("⚠️  同时指定了输入路径和样品映射文件，将优先使用样品映射文件")
            self.input_path = None
        
        # 创建输出路径
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径
        if self.input_path:
            self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.target_file = os.path.normpath(os.path.abspath(self.target_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        if self.sample_map_file:
            self.sample_map_file = os.path.normpath(os.path.abspath(self.sample_map_file))
        
        # 设置基础名称
        self.base_name = f"{self.blast_type}_analysis"
        
        # 验证BLAST类型
        if self.blast_type not in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']:
            raise ValueError(f"不支持的BLAST类型: {self.blast_type}")
        
        # 验证比对输出类型
        if self.alignment_output not in ['none', 'text', 'html', 'both']:
            raise ValueError(f"不支持的比对输出类型: {self.alignment_output}，必须是 'none', 'text', 'html', 或 'both'")
    
    def get_blast_executable(self) -> str:
        """获取BLAST可执行文件路径"""
        blast_paths = {
            'blastn': self.blastn_path,
            'blastp': self.blastp_path,
            'blastx': self.blastx_path,
            'tblastn': self.tblastn_path,
            'tblastx': self.tblastx_path
        }
        return blast_paths[self.blast_type]
    
    def get_db_type(self) -> str:
        """获取数据库类型"""
        db_types = {
            'blastn': 'nucl',
            'blastp': 'prot',
            'blastx': 'prot',
            'tblastn': 'nucl',
            'tblastx': 'nucl'
        }
        return db_types.get(self.blast_type, self.target_db_type)
    
    def needs_alignment_sequences(self) -> bool:
        """判断是否需要序列数据（用于比对可视化）"""
        return self.alignment_output in ['text', 'html', 'both']
    
    def get_final_outfmt(self) -> str:
        """获取最终的outfmt格式（自动添加qseq和sseq如果需要）"""
        base_outfmt = self.outfmt
        
        if self.needs_alignment_sequences():
            # 将格式字符串分割为字段列表进行精确检查
            fields = base_outfmt.split()
            
            # 检查是否已包含 qseq 和 sseq（精确匹配字段）
            if 'qseq' not in fields:
                base_outfmt += ' qseq'
            if 'sseq' not in fields:
                base_outfmt += ' sseq'
        
        return base_outfmt
    
    def validate(self):
        """验证配置参数"""
        errors = []
        
        # 检查输入数据参数
        if not self.input_path and not self.sample_map_file:
            errors.append("必须指定输入路径(-i)或样品映射文件(-s)中的一个")
        
        # 检查输入路径
        if self.input_path and not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在: {self.input_path}")
        
        # 检查目标文件
        if not os.path.exists(self.target_file):
            errors.append(f"目标基因文件不存在: {self.target_file}")
        
        # 检查样品映射文件
        if self.sample_map_file and not self.auto_generated_map and not os.path.exists(self.sample_map_file):
            errors.append(f"样品映射文件不存在: {self.sample_map_file}")
        
        # 检查参数范围
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数: {self.threads}")
        if not 0 < self.evalue <= 1:
            errors.append(f"E-value阈值必须在0-1之间: {self.evalue}")
        if not 0 <= self.min_identity <= 100:
            errors.append(f"最小相似度必须在0-100之间: {self.min_identity}")
        if not 0 <= self.min_coverage <= 100:
            errors.append(f"最小覆盖度必须在0-100之间: {self.min_coverage}")
        
        # 检查比对可视化参数
        if self.alignment_width < 40:
            errors.append(f"比对宽度必须至少为40: {self.alignment_width}")
        if self.alignment_max_per_sample < 1:
            errors.append(f"每个样品最大比对数必须为正整数: {self.alignment_max_per_sample}")
        
        if errors:
            raise ValueError("\n".join(errors))
        return True