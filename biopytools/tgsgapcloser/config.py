"""
TGS-GapCloser配置管理模块|TGS-GapCloser Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class TGSGapCloserConfig:
    """TGS-GapCloser配置类|TGS-GapCloser Configuration Class"""

    # 必需文件|Required files
    scaff_file: str  # 输入scaffold文件
    reads_file: str  # 输入TGS reads文件
    output_prefix: str  # 输出前缀

    # 工具路径|Tool paths
    tgsgapcloser_path: str = '~/software/TGS-GapCloser2/TGS-GapCloser2-master/tgsgapcloser2'

    # 纠错模式|Error correction mode
    mode: str = 'none'  # 'none', 'racon', 'pilon'

    # TGS类型|TGS type
    tgstype: str = 'ont'  # 'ont', 'pb', 'hifi'

    # 过滤参数|Filter parameters
    min_idy: Optional[float] = None  # 最小同一性，默认根据tgstype自动设置
    min_match: Optional[int] = None  # 最小匹配长度，默认根据tgstype自动设置

    # 线程数|Threads
    threads: int = 12

    # Racon参数|Racon parameters
    racon_path: Optional[str] = None
    racon_round: int = 3

    # Pilon参数|Pilon parameters
    pilon_path: Optional[str] = None
    pilon_mem: str = '300G'
    pilon_round: int = 3
    ngs_file: Optional[str] = None
    java_path: Optional[str] = None
    samtools_path: Optional[str] = None

    # 其他参数|Other parameters
    minmap_arg: Optional[str] = None  # 自定义minimap2参数
    chunk: int = 3  # 分块数量
    g_check: bool = False  # Gap大小差异检查
    min_nread: int = 1  # 最小reads数量
    max_nread: int = -1  # 最大reads数量
    max_candidate: int = 200  # 最大候选数

    # quarTeT gapfiller参数（第2轮填充）|quarTeT gapfiller parameters (2nd round)
    unitig_file: Optional[str] = None  # hifiasm unitig/contig文件
    quartet_path: Optional[str] = None  # quarTeT gapfiller路径
    flanking_len: int = 5000  # flanking序列长度
    min_align_len: int = 1000  # 最小比对长度
    min_identity: int = 40  # 最小比对同一性（百分比）
    max_filling_len: int = 1000000  # 最大填充长度

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.scaff_file = os.path.normpath(os.path.abspath(self.scaff_file))
        self.reads_file = os.path.normpath(os.path.abspath(self.reads_file))
        self.tgsgapcloser_path = os.path.normpath(os.path.abspath(os.path.expanduser(self.tgsgapcloser_path)))

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(self.output_prefix)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 自动设置tgstype相关参数|Auto-set tgstype related parameters
        if self.min_idy is None:
            if self.tgstype == 'ont':
                self.min_idy = 0.3
            else:  # pb or hifi
                self.min_idy = 0.2

        if self.min_match is None:
            if self.tgstype == 'ont':
                self.min_match = 300
            else:  # pb or hifi
                self.min_match = 200

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        required_files = [
            ('Scaffold文件|Scaffold file', self.scaff_file),
            ('TGS reads文件|TGS reads file', self.reads_file),
        ]

        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查TGS-GapCloser路径|Check TGS-GapCloser path
        if not os.path.exists(self.tgsgapcloser_path):
            errors.append(f"TGS-GapCloser不存在|TGS-GapCloser does not exist: {self.tgsgapcloser_path}")

        # 检查mode参数|Check mode parameter
        valid_modes = ['none', 'racon', 'pilon']
        if self.mode not in valid_modes:
            errors.append(f"无效的纠错模式|Invalid error correction mode: {self.mode} (必须是|must be one of {valid_modes})")

        # 检查tgstype参数|Check tgstype parameter
        valid_tgstypes = ['ont', 'pb', 'hifi']
        if self.tgstype not in valid_tgstypes:
            errors.append(f"无效的TGS类型|Invalid TGS type: {self.tgstype} (必须是|must be one of {valid_tgstypes})")

        # 检查Racon相关参数|Check Racon related parameters
        if self.mode == 'racon':
            if self.racon_path and not os.path.exists(self.racon_path):
                errors.append(f"Racon路径不存在|Racon path does not exist: {self.racon_path}")

        # 检查Pilon相关参数|Check Pilon related parameters
        if self.mode == 'pilon':
            if self.pilon_path and not os.path.exists(self.pilon_path):
                errors.append(f"Pilon路径不存在|Pilon path does not exist: {self.pilon_path}")
            if self.ngs_file and not os.path.exists(self.ngs_file):
                errors.append(f"NGS reads文件不存在|NGS reads file does not exist: {self.ngs_file}")
            if self.samtools_path and not os.path.exists(self.samtools_path):
                errors.append(f"Samtools路径不存在|Samtools path does not exist: {self.samtools_path}")
            if self.java_path and not os.path.exists(self.java_path):
                errors.append(f"Java路径不存在|Java path does not exist: {self.java_path}")

        # 检查quarTeT相关参数|Check quarTeT related parameters
        if self.unitig_file and not os.path.exists(self.unitig_file):
            errors.append(f"Unitig文件不存在|Unitig file does not exist: {self.unitig_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
