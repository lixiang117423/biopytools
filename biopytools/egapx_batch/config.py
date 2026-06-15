"""
EGAPx批量配置类|EGAPx Batch Configuration Class
"""

import os
from dataclasses import dataclass, field
from typing import List
from ..common.paths import expand_path


@dataclass
class EGAPxBatchConfig:
    """EGAPx批量运行配置类|EGAPx Batch Processing Configuration Class"""

    # 必需参数|Required parameters
    genome: str
    output_dir: str

    # 可选参数|Optional parameters
    taxid: str = "71234"
    egapx_path: str = "~/software/EGAPX_v.0.4.1-alpha/egapx"
    local_cache: str = "~/software/EGAPX_v.0.4.1-alpha/local_cache"
    sif_image: str = "~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif"
    split_genome: bool = True
    chr_prefix: str = None
    locus_tag_prefix: str = ""
    report_name: str = "EGAPx"
    short_reads: str = ""
    long_reads: str = ""

    # 运行时变量|Runtime variables
    chromosomes: List[str] = field(default_factory=list)
    job_count: int = 0

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.genome = expand_path(self.genome)
        self.egapx_path = expand_path(self.egapx_path)

        # 展开可选路径|Expand optional paths
        if self.local_cache:
            self.local_cache = expand_path(self.local_cache)
        self.sif_image = expand_path(self.sif_image)
        if self.short_reads:
            self.short_reads = expand_path(self.short_reads)
        if self.long_reads:
            self.long_reads = expand_path(self.long_reads)

        # 创建输出目录|Create output directory
        self.output_dir = expand_path(self.output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时|When configuration parameters are invalid
        """
        errors = []

        # 检查必需文件|Check required files
        if not self.genome or not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        if not self.egapx_path or not os.path.exists(self.egapx_path):
            errors.append(f"EGAPx路径不存在|EGAPx path not found: {self.egapx_path}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"EGAPxBatchConfig(\n"
            f"  genome={self.genome!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  taxid={self.taxid!r},\n"
            f"  egapx_path={self.egapx_path!r},\n"
            f"  chr_prefix={self.chr_prefix!r},\n"
            f"  locus_tag_prefix={self.locus_tag_prefix!r},\n"
            f"  report_name={self.report_name!r}\n"
            f")"
        )
