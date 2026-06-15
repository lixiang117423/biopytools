"""
端粒识别配置管理模块|Telomere Finder Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class TelomereFinderConfig:
    """端粒识别配置类|Telomere Finder Configuration Class"""

    # 必需文件|Required files
    genome_file: str

    # 输出配置|Output configuration
    output_dir: str = './telomere_output'
    output_prefix: str = 'telomere'

    # 软件路径配置|Software path configuration
    tidk_path: str = '~/miniforge3/envs/tidk_v.0.2.65/bin/tidk'

    # 分析模式|Analysis mode
    mode: str = 'pipeline'  # explore, find, search, plot, pipeline (默认 pipeline)

    # Explore 模式参数|Explore mode parameters
    explore_min_length: int = 5
    explore_max_length: int = 12
    explore_threshold: int = 100
    explore_distance: float = 0.01

    # Find 模式参数|Find mode parameters
    clade: Optional[str] = None
    window_size: int = 10000

    # Search 模式参数|Search mode parameters
    search_string: Optional[str] = None
    output_format: str = 'tsv'  # tsv, bedgraph

    # Plot 模式参数|Plot mode parameters
    tsv_file: Optional[str] = None
    plot_height: int = 200
    plot_width: int = 1000
    plot_fontsize: int = 12
    plot_strokewidth: int = 2

    # 通用参数|General parameters
    threads: int = 1
    verbose: bool = False
    log_file: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.tidk_path = os.path.normpath(os.path.abspath(expand_path(self.tidk_path)))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 如果指定了 TSV 文件，标准化路径|If TSV file specified, normalize path
        if self.tsv_file:
            self.tsv_file = os.path.normpath(os.path.abspath(self.tsv_file))

        # 验证模式|Validate mode
        valid_modes = ['explore', 'find', 'search', 'plot', 'pipeline']
        if self.mode not in valid_modes:
            raise ValueError(f"无效的模式|Invalid mode: {self.mode}. 必须是|Must be one of {valid_modes}")

        # 验证输出格式|Validate output format
        valid_formats = ['tsv', 'bedgraph']
        if self.output_format not in valid_formats:
            raise ValueError(f"无效的输出格式|Invalid output format: {self.output_format}. 必须是|Must be one of {valid_formats}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_file}")

        # 检查 tidk 路径|Check tidk path
        if not os.path.exists(self.tidk_path):
            errors.append(f"tidk 路径不存在|tidk path does not exist: {self.tidk_path}")

        # 模式特定验证|Mode-specific validation
        if self.mode == 'find' and not self.clade:
            errors.append("find 模式需要指定 --clade 参数|find mode requires --clade parameter")

        if self.mode == 'search' and not self.search_string:
            errors.append("search 模式需要指定 --search-string 参数|search mode requires --search-string parameter")

        if self.mode == 'plot' and not self.tsv_file:
            errors.append("plot 模式需要指定 --tsv 参数|plot mode requires --tsv parameter")

        if self.mode == 'plot' and self.tsv_file and not os.path.exists(self.tsv_file):
            errors.append(f"TSV 文件不存在|TSV file does not exist: {self.tsv_file}")

        # 验证 explore 距离参数|Validate explore distance parameter
        if not 0 <= self.explore_distance <= 0.5:
            errors.append(f"explore_distance 必须在 0-0.5 之间|explore_distance must be between 0-0.5: {self.explore_distance}")

        # 验证窗口大小|Validate window size
        if self.window_size <= 0:
            errors.append(f"window_size 必须大于 0|window_size must be > 0: {self.window_size}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
