"""
PanDepth覆盖度计算配置管理模块|PanDepth Coverage Calculation Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class PanDepthConfig:
    """PanDepth覆盖度计算配置类|PanDepth Coverage Calculation Configuration Class"""

    # 必需参数|Required parameters
    input_path: str  # BAM文件或BAM文件目录|BAM file or BAM file directory
    output_dir: str  # 输出目录|Output directory

    # 路径配置|Path configuration
    pandepth_path: str = '~/software/PanDepth-2.26-Linux-x86_64/pandepth'

    # 目标区域选项|Target region options
    gff_file: Optional[str] = None  # GFF/GTF文件用于基因覆盖度|GFF/GTF file for gene coverage
    bed_file: Optional[str] = None  # BED文件用于特定区域|BED file for specific regions
    window_size: Optional[int] = None  # 滑动窗口大小(bp)|Sliding window size in bp

    # GFF/GTF特征类型|GFF/GTF feature type
    feature_type: str = 'CDS'  # CDS 或 exon|CDS or exon

    # 过滤选项|Filter options
    min_mapq: int = 0  # 最小比对质量|Minimum mapping quality
    min_depth: int = 1  # 最小深度用于统计|Minimum depth for statistics
    exclude_flag: int = 1796  # 排除reads的FLAG标志|FLAG bits to exclude reads

    # 其他选项|Other options
    threads: int = 12  # 线程数|Number of threads
    reference: Optional[str] = None  # 参考基因组文件(用于CRAM解码或GC计算)|Reference genome file for CRAM decode or GC calculation
    enable_gc: bool = False  # 启用GC含量计算|Enable GC content calculation
    output_all_sites: bool = False  # 输出所有位点深度|Output all site depths

    # 日志选项|Logging options
    verbose: bool = False
    quiet: bool = False
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.pandepth_path = os.path.normpath(os.path.abspath(expand_path(self.pandepth_path)))

        if self.gff_file:
            self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))
        if self.bed_file:
            self.bed_file = os.path.normpath(os.path.abspath(self.bed_file))
        if self.reference:
            self.reference = os.path.normpath(os.path.abspath(self.reference))

        # 判断输入是文件还是目录|Determine if input is file or directory
        self.is_batch_mode = os.path.isdir(self.input_path)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")

        # 检查PanDepth程序|Check PanDepth program
        if not os.path.exists(self.pandepth_path):
            errors.append(f"PanDepth程序不存在|PanDepth program not found: {self.pandepth_path}")

        # 检查目标区域选项|Check target region options
        target_options = [self.gff_file, self.bed_file, self.window_size]
        if sum(option is not None for option in target_options) > 1:
            errors.append(
                "只能指定一个目标区域选项(-g, -b, -w)|"
                "Only one target region option can be specified (-g, -b, -w)"
            )

        # 检查GFF/GTF文件|Check GFF/GTF file
        if self.gff_file and not os.path.exists(self.gff_file):
            errors.append(f"GFF/GTF文件不存在|GFF/GTF file does not exist: {self.gff_file}")

        # 检查BED文件|Check BED file
        if self.bed_file and not os.path.exists(self.bed_file):
            errors.append(f"BED文件不存在|BED file does not exist: {self.bed_file}")

        # 检查参考基因组|Check reference genome
        if self.reference and not os.path.exists(self.reference):
            errors.append(f"参考基因组文件不存在|Reference genome file does not exist: {self.reference}")

        # 检查特征类型|Check feature type
        if self.feature_type not in ['CDS', 'exon']:
            errors.append(f"无效的特征类型|Invalid feature type: {self.feature_type} (必须为CDS或exon|must be CDS or exon)")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        # 检查窗口大小|Check window size
        if self.window_size is not None and self.window_size <= 0:
            errors.append(f"窗口大小必须为正整数|Window size must be positive integer: {self.window_size}")

        # 如果启用GC计算，必须提供参考基因组|If GC calculation is enabled, reference genome must be provided
        if self.enable_gc and not self.reference:
            errors.append("启用GC计算时必须提供参考基因组文件|Reference genome file required when GC calculation is enabled")

        if errors:
            raise ValueError("\n".join(errors))

        return True
