"""
叶绿体基因组组装配置管理模块|Plastome Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class PlastomeConfig:
    """叶绿体基因组组装配置类|Plastome Assembly Configuration Class"""

    # 必需文件/目录|Required files/directories
    input_dir: str

    # 路径配置|Path configuration
    getorganelle_path: str = '~/miniforge3/envs/getorganelle_v.1.7.71/bin/get_organelle_from_reads.py'
    output_dir: str = './plastome_output'

    # GetOrganelle参数|GetOrganelle parameters
    organelle_type: str = 'embplant_pt'  # 叶绿体基因组|Plastome
    max_rounds: int = 15  # 最大扩展轮数|Maximum extension rounds
    kmer_list: str = '21,45,65,85,105'  # kmer列表|Kmer list
    threads: int = 12  # 线程数|Threads

    # Reads文件后缀模式|Reads file suffix patterns (支持通配符|Supports wildcards)
    read1_suffix: str = '_1.clean.fq.gz'  # R1文件后缀|R1 file suffix
    read2_suffix: str = '_2.clean.fq.gz'  # R2文件后缀|R2 file suffix

    # 自动检测的文件|Auto-detected files
    r1_file: Optional[str] = None
    r2_file: Optional[str] = None
    unpaired_files: Optional[str] = None

    # 输出前缀|Output prefix
    output_prefix: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径（保持相对路径以避免中文路径问题）|Normalize paths (keep relative to avoid Chinese character issues)
        self.input_dir = os.path.normpath(self.input_dir)
        # getorganelle_path 需要绝对路径，因为它是外部工具|getorganelle_path needs absolute path as it's an external tool
        self.getorganelle_path = os.path.normpath(os.path.abspath(expand_path(self.getorganelle_path)))
        self.output_dir = os.path.normpath(self.output_dir)

        # 如果未指定输出前缀，使用输入目录名|If output prefix not specified, use input directory name
        if self.output_prefix is None:
            self.output_prefix = os.path.basename(self.input_dir)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        # 检查GetOrganelle路径|Check GetOrganelle path
        if not os.path.exists(self.getorganelle_path):
            errors.append(f"GetOrganelle路径不存在|GetOrganelle path does not exist: {self.getorganelle_path}")

        # 检查reads文件|Check reads files
        if self.r1_file and not os.path.exists(self.r1_file):
            errors.append(f"R1文件不存在|R1 file does not exist: {self.r1_file}")

        if self.r2_file and not os.path.exists(self.r2_file):
            errors.append(f"R2文件不存在|R2 file does not exist: {self.r2_file}")

        if self.unpaired_files:
            for f in self.unpaired_files.split(','):
                if not os.path.exists(f.strip()):
                    errors.append(f"未配对文件不存在|Unpaired file does not exist: {f.strip()}")

        # 验证organelle类型|Validate organelle type
        valid_types = ['embplant_pt', 'embplant_mt', 'embplant_nr', 'other_pt',
                      'animal_mt', 'fungus_mt', 'fungus_nr']
        if self.organelle_type not in valid_types:
            errors.append(f"无效的organelle类型|Invalid organelle type: {self.organelle_type}")

        # 验证线程数|Validate threads
        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Threads must be > 0: {self.threads}")

        # 验证kmer列表|Validate kmer list
        try:
            kmers = [int(k.strip()) for k in self.kmer_list.split(',')]
            if any(k <= 0 for k in kmers):
                errors.append(f"Kmer值必须大于0|Kmer values must be > 0")
        except ValueError:
            errors.append(f"无效的kmer列表格式|Invalid kmer list format: {self.kmer_list}")

        # 检查是否找到了reads文件|Check if reads files were found
        if not self.r1_file and not self.unpaired_files:
            errors.append("未找到任何reads文件|No reads files found in input directory")

        if errors:
            raise ValueError("\n".join(errors))

        return True
