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

    # 输入配置|Input configuration
    input: str  # 输入路径（BAM文件或包含BAM的目录）| Input path (BAM file or directory containing BAM files)
    bed_file: Optional[str] = None  # BED文件路径，提供批量区间|BED file path for batch intervals

    # 区域配置（BED模式下可选）|Region configuration (optional in BED mode)
    chromosome: Optional[str] = None  # 染色体名称|Chromosome name
    start: Optional[int] = None  # 起始位置（1-based）|Start position (1-based)

    # 路径配置|Path configuration
    output: str = './coverage.txt'  # 输出文件路径|Output file path

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
        # 解析 -o：目录型(已存在目录 或 无扩展名) → 合并文件落 目录/coverage.txt；
        # 文件型(有扩展名) → 直接作为输出文件(向后兼容)
        # |Resolve -o: dir-like (existing dir or no extension) → merged at dir/coverage.txt;
        # file-like (has extension) → use as output file (backward compatible).
        # 否则直接把目录当文件写会触发 [Errno 21] Is a directory。
        # |Otherwise writing to a dir path raises [Errno 21] Is a directory.
        output_abs = os.path.normpath(os.path.abspath(self.output))
        has_extension = '.' in os.path.basename(output_abs)
        if os.path.isdir(output_abs) or not has_extension:
            # 目录型：自动创建，合并文件/临时文件/摘要都落在目录内
            # |dir-like: auto-create; merged/temp/summary all land inside the dir
            os.makedirs(output_abs, exist_ok=True)
            self.output = os.path.join(output_abs, 'coverage.txt')
        else:
            self.output = output_abs
        self.output_path = Path(self.output).parent
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化输入路径|Normalize input path
        self.input = os.path.normpath(os.path.abspath(self.input))

        # 标准化BED文件路径|Normalize BED file path
        if self.bed_file:
            self.bed_file = os.path.normpath(os.path.abspath(self.bed_file))

        # 验证BED模式与手动模式的互斥性|Validate mutual exclusivity of BED mode vs manual mode
        if self.bed_file:
            if not os.path.isfile(self.bed_file):
                raise ValueError(f"BED文件不存在|BED file does not exist: {self.bed_file}")
            if self.chromosome is not None or self.start is not None:
                raise ValueError(
                    "BED模式下不需要指定--chromosome和--start参数|"
                    "--chromosome and --start are not needed in BED mode"
                )
        else:
            if self.chromosome is None or self.start is None:
                raise ValueError(
                    "非BED模式下必须指定--chromosome和--start参数|"
                    "--chromosome and --start are required when --bed is not provided"
                )

        # 自动识别输入类型|Auto-detect input type
        self._is_directory = os.path.isdir(self.input)
        self._is_file = os.path.isfile(self.input)

        # 验证输入|Validate input
        if not self._is_directory and not self._is_file:
            raise ValueError(f"输入路径不存在|Input path does not exist: {self.input}")

        # 如果是文件，验证是BAM文件|If it's a file, validate it's a BAM file
        if self._is_file and not self.input.endswith('.bam'):
            raise ValueError(f"输入文件必须是BAM格式 (.bam)|Input file must be BAM format (.bam): {self.input}")

        # 验证位置参数（仅手动模式）|Validate position parameters (manual mode only)
        if not self.bed_file:
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

    def is_bed_mode(self) -> bool:
        """是否为BED批量模式|Whether in BED batch mode"""
        return self.bed_file is not None

    def is_file(self) -> bool:
        """判断输入是否为文件|Check if input is a file"""
        return self._is_file
