"""
BAM统计分析配置模块|BAM Statistics Analysis Configuration Module
"""

import os
import argparse
from pathlib import Path
from typing import Optional, List

class BAMStatsConfig:
    """BAM统计分析配置类|BAM Statistics Analysis Configuration Class"""

    def __init__(self, **kwargs):
        """初始化配置参数|Initialize configuration parameters"""

        # 输入文件或文件夹|Input file or directory
        self.input_path: str = kwargs.get('input_path', '')
        self.bam_files: List[str] = []
        self.reference_file: Optional[str] = kwargs.get('reference_file', None)
        self.bed_file: Optional[str] = kwargs.get('bed_file', None)

        # 输出设置|Output settings
        self.output_dir: str = kwargs.get('output_dir', 'bam_stats_output')
        self.prefix: str = kwargs.get('prefix', 'sample')

        # 分析参数|Analysis parameters
        self.min_mapq: int = kwargs.get('min_mapq', 20)
        self.min_base_quality: int = kwargs.get('min_base_quality', 20)
        self.coverage_bins: int = kwargs.get('coverage_bins', 100)
        self.window_size: int = kwargs.get('window_size', 1000000)  # 1M窗口大小|1M window size
        self.step_size: int = kwargs.get('step_size', 100000)      # 100kb步长|100kb step size
        self.max_insert_size: int = kwargs.get('max_insert_size', 1000)

        # 分析模块开关|Analysis module switches
        self.skip_alignment_stats: bool = kwargs.get('skip_alignment_stats', False)
        self.skip_quality_stats: bool = kwargs.get('skip_quality_stats', False)
        self.skip_coverage_stats: bool = kwargs.get('skip_coverage_stats', False)
        self.skip_sequence_stats: bool = kwargs.get('skip_sequence_stats', False)
        self.skip_insert_stats: bool = kwargs.get('skip_insert_stats', False)
        self.skip_duplicate_stats: bool = kwargs.get('skip_duplicate_stats', False)
        self.skip_chromosome_stats: bool = kwargs.get('skip_chromosome_stats', False)

        # 可视化选项|Visualization options
        self.generate_plots: bool = kwargs.get('generate_plots', False)
        self.plot_format: str = kwargs.get('plot_format', 'png')

        # 工具路径|Tool paths
        self.samtools_path: str = kwargs.get('samtools_path', 'samtools')
        self.bedtools_path: str = kwargs.get('bedtools_path', 'bedtools')

        # 性能设置|Performance settings
        self.threads: int = kwargs.get('threads', 88)
        self.memory_limit: str = kwargs.get('memory_limit', '4G')
        self.parallel_samples: bool = kwargs.get('parallel_samples', True)
        self.max_workers: int = kwargs.get('max_workers', min(self.threads, 16))  # 限制并发样品数|Limit concurrent samples

        # 处理路径|Process paths
        self._process_paths()

    def _process_paths(self):
        """处理文件路径|Process file paths"""
        # 发现BAM文件|Discover BAM files
        self._discover_bam_files()

        # 规范化参考基因组路径|Normalize reference file path
        if self.reference_file:
            self.reference_file = os.path.normpath(os.path.abspath(self.reference_file))

        # 规范化BED文件路径|Normalize BED file path
        if self.bed_file:
            self.bed_file = os.path.normpath(os.path.abspath(self.bed_file))

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(self.output_path.absolute())

    def _discover_bam_files(self):
        """发现BAM文件|Discover BAM files"""
        if not self.input_path:
            return

        input_path = Path(self.input_path).absolute()

        if input_path.is_file():
            # 单个文件|Single file
            if input_path.suffix.lower() in ['.bam', '.sam']:
                self.bam_files = [str(input_path)]
        elif input_path.is_dir():
            # 文件夹，递归查找BAM文件|Directory, recursively find BAM files
            bam_patterns = ['*.bam', '*.BAM']
            bam_files = []
            for pattern in bam_patterns:
                bam_files.extend(input_path.rglob(pattern))
            self.bam_files = [str(f) for f in sorted(bam_files)]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input path
        if not self.input_path:
            errors.append("必须指定输入路径|Input path must be specified")
        elif not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")

        # 检查发现的BAM文件|Check discovered BAM files
        if not self.bam_files:
            errors.append(f"在输入路径中未找到BAM文件|No BAM files found in input path: {self.input_path}")

        # 检查参考基因组文件|Check reference file
        if self.reference_file and not os.path.exists(self.reference_file):
            errors.append(f"参考基因组文件不存在|Reference file does not exist: {self.reference_file}")

        # 检查BED文件|Check BED file
        if self.bed_file and not os.path.exists(self.bed_file):
            errors.append(f"BED文件不存在|BED file does not exist: {self.bed_file}")

        # 检查参数范围|Check parameter ranges
        if self.min_mapq < 0 or self.min_mapq > 60:
            errors.append(f"MAPQ阈值必须在0-60之间|MAPQ threshold must be between 0-60: {self.min_mapq}")

        if self.min_base_quality < 0 or self.min_base_quality > 60:
            errors.append(f"碱基质量阈值必须在0-60之间|Base quality threshold must be between 0-60: {self.min_base_quality}")

        if self.coverage_bins <= 0:
            errors.append(f"覆盖度分箱数必须为正整数|Coverage bins must be positive: {self.coverage_bins}")

        if self.max_insert_size <= 0:
            errors.append(f"最大插入大小必须为正整数|Max insert size must be positive: {self.max_insert_size}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if self.window_size <= 0:
            errors.append(f"滑窗大小必须为正整数|Window size must be positive: {self.window_size}")

        if self.step_size <= 0:
            errors.append(f"滑窗步长必须为正整数|Step size must be positive: {self.step_size}")

        if self.step_size > self.window_size:
            errors.append(f"滑窗步长不能大于窗口大小|Step size cannot be larger than window size: {self.step_size} > {self.window_size}")

        if self.plot_format not in ['png', 'pdf', 'svg']:
            errors.append(f"不支持的图片格式|Unsupported plot format: {self.plot_format}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
