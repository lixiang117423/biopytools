"""
SRA转换配置管理模块 |SRA Conversion Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class ConvertConfig:
    """转换配置类|Conversion Configuration Class"""
    
    # 输入输出|Input/Output
    input_path: str  # 可以是文件或文件夹|Can be file or folder
    output_dir: str = './fastq_output'
    
    # 转换参数|Conversion parameters
    compress: bool = True  # 是否压缩输出|Compress output
    threads: int = 12  # 线程数|Number of threads
    split_files: bool = True  # 是否拆分双端测序|Split paired-end reads
    
    # 工具选择|Tool selection
    use_parallel: bool = True  # 优先使用parallel-fastq-dump|Prefer parallel-fastq-dump
    tool_path: str = 'parallel-fastq-dump'  # 工具路径|Tool path
    
    # 额外参数|Additional parameters
    skip_technical: bool = True  # 跳过技术序列|Skip technical reads
    clip: bool = False  # 剪切adapters|Clip adapters
    min_read_len: int = 0  # 最小读长过滤|Minimum read length filter
    tmpdir: Optional[str] = None  # 临时目录|Temporary directory
    
    # 内部属性|Internal attributes
    input_files: List[str] = None
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 如果指定了临时目录，创建它|Create tmpdir if specified
        if self.tmpdir:
            Path(self.tmpdir).mkdir(parents=True, exist_ok=True)
        
        # 确定输入文件列表|Determine input file list
        self._determine_input_files()
    
    # def _determine_input_files(self):
    #     """确定输入文件列表|Determine input file list"""
    #     input_p = Path(self.input_path)
        
    #     if input_p.is_file():
    #         # 单个文件|Single file
    #         if input_p.suffix == '.sra':
    #             self.input_files = [str(input_p)]
    #         else:
    #             raise ValueError(f"输入文件必须是.sra格式|Input file must be .sra format: {self.input_path}")
    #     elif input_p.is_dir():
    #         # 文件夹，查找所有.sra文件|Folder, find all .sra files
    #         self.input_files = sorted([str(f) for f in input_p.glob('*.sra')])
    #         if not self.input_files:
    #             raise ValueError(f"在目录中未找到.sra文件|No .sra files found in directory: {self.input_path}")
    #     else:
    #         raise ValueError(f"输入路径不存在|Input path does not exist: {self.input_path}")

    def _determine_input_files(self):
        """确定输入文件列表|Determine input file list"""
        input_p = Path(self.input_path)
        
        if input_p.is_file():
            # 单个文件|Single file
            # 接受带或不带.sra后缀的文件|Accept files with or without .sra suffix
            if input_p.suffix == '.sra' or input_p.suffix == '':
                self.input_files = [str(input_p)]
            else:
                raise ValueError(f"输入文件格式不正确|Invalid input file format: {self.input_path}")
        elif input_p.is_dir():
            # 文件夹，查找所有目标文件|Folder, find all target files
            # 先查找.sra文件|First find .sra files
            sra_files = list(input_p.glob('*.sra'))
            
            if sra_files:
                # 如果有.sra文件，使用它们|If .sra files exist, use them
                self.input_files = sorted([str(f) for f in sra_files])
            else:
                # 如果没有.sra文件，查找所有文件（排除隐藏文件和目录）
                # If no .sra files, find all files (exclude hidden files and directories)
                all_files = [f for f in input_p.iterdir() 
                            if f.is_file() and not f.name.startswith('.')]
                
                if not all_files:
                    raise ValueError(f"在目录中未找到文件|No files found in directory: {self.input_path}")
                
                self.input_files = sorted([str(f) for f in all_files])
        else:
            raise ValueError(f"输入路径不存在|Input path does not exist: {self.input_path}")
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查输入路径|Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")
        
        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Number of threads must be positive: {self.threads}")
        
        if self.min_read_len < 0:
            errors.append(f"最小读长不能为负数|Minimum read length cannot be negative: {self.min_read_len}")
        
        if not self.input_files:
            errors.append("没有找到可转换的.sra文件|No .sra files found for conversion")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
