"""
BAM覆盖度统计配置管理模块|BAM Coverage Statistics Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class BAMCoverageConfig:
    """BAM覆盖度统计配置类|BAM Coverage Statistics Configuration Class"""

    # 必需输入|Required inputs
    chromosome: str  # 染色体名称|Chromosome name
    start: int  # 起始位置|Start position
    input: str  # 输入路径（BAM文件或包含BAM的目录）| Input path (BAM file or directory containing BAM files)

    # 路径配置|Path configuration
    output_dir: str = './bam_coverage_stats_output'
    output_prefix: str = 'coverage'

    # 过滤参数|Filtering parameters
    min_mapq: int = 0  # 最小mapping质量|Minimum mapping quality
    min_baseq: int = 0  # 最小碱基质量|Minimum base quality

    # 输出选项|Output options
    end: Optional[int] = None  # 终止位置（None表示到染色体末端）| End position (None means to chromosome end)
    merge_output: bool = True  # 是否合并所有样本的输出|Whether to merge all samples' output
    generate_summary: bool = True  # 是否生成统计摘要|Whether to generate summary statistics

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    # 高级选项|Advanced options
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    threads: int = 64  # 线程数|Number of threads

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化输入路径|Normalize input path
        self.input = os.path.normpath(os.path.abspath(self.input))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 自动识别输入类型|Auto-detect input type
        self._is_directory = os.path.isdir(self.input)
        self._is_file = os.path.isfile(self.input)

        # 验证输入|Validate input
        if not self._is_directory and not self._is_file:
            raise ValueError(f"输入路径不存在|Input path does not exist: {self.input}")

        # 如果是文件，验证是BAM文件|If it's a file, validate it's a BAM file
        if self._is_file and not self.input.endswith('.bam'):
            raise ValueError(f"输入文件必须是BAM格式 (.bam)|Input file must be BAM format (.bam): {self.input}")

        # 验证位置参数|Validate position parameters
        if self.start < 0:
            raise ValueError(f"起始位置必须非负|Start position must be non-negative: {self.start}")

        if self.end is not None and self.end <= self.start:
            raise ValueError(f"终止位置必须大于起始位置|End position must be greater than start: {self.start} - {self.end}")

        # 验证质量参数|Validate quality parameters
        if self.min_mapq < 0 or self.min_mapq > 255:
            raise ValueError(f"最小mapping质量必须在0-255之间|Min MAPQ must be between 0-255: {self.min_mapq}")

        if self.min_baseq < 0 or self.min_baseq > 93:
            raise ValueError(f"最小碱基质量必须在0-93之间|Min base quality must be between 0-93: {self.min_baseq}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        # 输入路径已在__post_init__中验证|Input path already validated in __post_init__
        return True

    def get_bam_files(self) -> List[str]:
        """获取BAM文件列表|Get list of BAM files"""
        if self._is_file:
            # 单个BAM文件|Single BAM file
            return [self.input]
        else:
            # 从目录获取所有BAM文件|Get all BAM files from directory
            bam_files = []
            for file in os.listdir(self.input):
                if file.endswith('.bam'):
                    bam_path = os.path.join(self.input, file)
                    if os.path.isfile(bam_path):
                        bam_files.append(bam_path)

            if not bam_files:
                raise ValueError(f"目录中没有找到BAM文件|No BAM files found in directory: {self.input}")

            return sorted(bam_files)

    def is_directory(self) -> bool:
        """判断输入是否为目录|Check if input is a directory"""
        return self._is_directory

    def is_file(self) -> bool:
        """判断输入是否为文件|Check if input is a file"""
        return self._is_file
