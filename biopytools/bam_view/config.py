"""
BAM比对可视化配置管理模块|BAM Alignment Visualization Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class BamViewConfig:
    """BAM比对可视化配置类|BAM Alignment Visualization Configuration Class"""

    # 必需参数|Required parameters
    bam_file: str
    reference: str
    region: str

    # 软件配置|Software configuration
    alignoth_path: str = '~/miniforge3/envs/alignoth/bin/alignoth'

    # 输出配置|Output configuration
    output_dir: str = './bam_view_output'
    output_format: str = 'html'  # html, json, svg, pdf

    # 可视化参数|Visualization parameters
    max_read_depth: int = 500
    max_width: int = 1024
    mismatch_display_min_percent: float = 1.0

    # 高亮配置|Highlight configuration
    vcf_file: Optional[str] = None
    bed_file: Optional[str] = None
    highlight_intervals: Optional[List[str]] = None

    # 其他选项|Other options
    aux_tags: Optional[List[str]] = None
    no_embed_js: bool = False
    plot_all: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.bam_file = os.path.normpath(os.path.abspath(self.bam_file))
        self.reference = os.path.normpath(os.path.abspath(self.reference))

        if self.vcf_file:
            self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))

        if self.bed_file:
            self.bed_file = os.path.normpath(os.path.abspath(self.bed_file))

        self.alignoth_path = os.path.normpath(os.path.abspath(expand_path(self.alignoth_path)))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 验证输出格式|Validate output format
        valid_formats = ['html', 'json', 'svg', 'pdf']
        if self.output_format not in valid_formats:
            raise ValueError(
                f"无效的输出格式|Invalid output format: {self.output_format}. "
                f"支持的格式|Supported formats: {', '.join(valid_formats)}"
            )

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.bam_file):
            errors.append(f"BAM文件不存在|BAM file not found: {self.bam_file}")

        if not os.path.exists(self.reference):
            errors.append(f"参考序列文件不存在|Reference file not found: {self.reference}")

        # 检查BAM索引|Check BAM index
        bam_index = f"{self.bam_file}.bai"
        if not os.path.exists(bam_index):
            errors.append(f"BAM索引文件不存在|BAM index file not found: {bam_index}")

        # 检查可选文件|Check optional files
        if self.vcf_file and not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file not found: {self.vcf_file}")

            # 检查VCF索引|Check VCF index
            vcf_index = f"{self.vcf_file}.csi"
            vcf_tbi = f"{self.vcf_file}.tbi"
            if not os.path.exists(vcf_index) and not os.path.exists(vcf_tbi):
                errors.append(
                    f"VCF索引文件不存在|VCF index file not found: {vcf_index} 或|or {vcf_tbi}. "
                    f"请使用bcftools index或tabix创建索引|Please create index using bcftools index or tabix"
                )

        if self.bed_file and not os.path.exists(self.bed_file):
            errors.append(f"BED文件不存在|BED file not found: {self.bed_file}")

        # 检查alignoth|Check alignoth
        if not os.path.exists(self.alignoth_path):
            errors.append(
                f"alignoth不存在|alignoth not found: {self.alignoth_path}. "
                f"请先安装alignoth|Please install alignoth first"
            )

        # 检查参数范围|Check parameter ranges
        if self.max_read_depth <= 0:
            errors.append(f"max_read_depth必须为正数|max_read_depth must be positive: {self.max_read_depth}")

        if self.max_width <= 0:
            errors.append(f"max_width必须为正数|max_width must be positive: {self.max_width}")

        if not (0 <= self.mismatch_display_min_percent <= 100):
            errors.append(
                f"mismatch_display_min_percent必须在0-100之间|"
                f"mismatch_display_min_percent must be between 0-100: {self.mismatch_display_min_percent}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
