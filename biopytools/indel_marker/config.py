"""INDEL分子标记开发配置模块|INDEL Marker Configuration Module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class IndelMarkerConfig:
    """INDEL分子标记配置类|INDEL Marker Configuration Class"""

    # 必需输入|Required inputs
    vcf_file: str
    samplesheet: str
    genome_fasta: str

    # 输出|Output
    output_dir: str = "./indel_marker_output"

    # 通用|General
    threads: int = 12

    # INDEL过滤|INDEL filtering
    min_indel_size: int = 10
    max_indel_size: int = 100
    min_quality: float = 20.0

    # 候选数限制|Candidate cap (0 = no limit)
    max_candidates: int = 0

    # 群体判定|Population calling
    min_group_consistency: float = 0.9   # 组内纯合一致比例|within-group homozygous consistency
    min_samples_per_group: int = 1
    max_missing_rate: float = 0.2

    # 覆盖度|Coverage
    min_depth: int = 10
    min_mapq: int = 20
    min_baseq: int = 15
    deletion_depth_ratio: float = 0.3    # 骤降阈值|drop threshold
    deletion_size_for_coverage_check: int = 30

    # 引物|Primer
    flank_length: int = 300
    primer_product_min: int = 100
    primer_product_max: int = 600

    # 工具路径|Tool paths
    bcftools_path: str = field(
        default_factory=lambda: get_tool_path('bcftools', 'bcftools', 'BCFTOOLS_PATH')
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path('samtools', 'samtools', 'SAMTOOLS_PATH')
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.vcf_file = expand_path(self.vcf_file)
        self.samplesheet = expand_path(self.samplesheet)
        self.genome_fasta = expand_path(self.genome_fasta)
        self.output_dir = expand_path(self.output_dir)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 创建输出目录|Create output dirs
        self.output_path = Path(self.output_dir)
        self.pipeline_info_dir = self.output_path / "00_pipeline_info"
        self.vcf_extract_dir = self.output_path / "01_vcf_extract"
        self.genotype_dir = self.output_path / "02_genotype_call"
        self.coverage_dir = self.output_path / "03_coverage"
        self.sequence_dir = self.output_path / "04_sequence"
        self.primer_dir = self.output_path / "05_primer"
        self.results_dir = self.output_path / "06_results"
        self.logs_dir = self.output_path / "99_logs"

        for d in [self.pipeline_info_dir, self.vcf_extract_dir, self.genotype_dir,
                  self.coverage_dir, self.sequence_dir, self.primer_dir,
                  self.results_dir, self.logs_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 必需文件|Required files
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file not found: {self.vcf_file}")
        if not os.path.exists(self.samplesheet):
            errors.append(f"samplesheet不存在|samplesheet not found: {self.samplesheet}")
        if not os.path.exists(self.genome_fasta):
            errors.append(f"基因组文件不存在|Genome FASTA not found: {self.genome_fasta}")

        # 参数范围|Parameter ranges
        if self.min_indel_size < 1:
            errors.append(f"最小INDEL长度必须>=1|min_indel_size must be >=1: {self.min_indel_size}")
        if self.max_indel_size < self.min_indel_size:
            errors.append(f"最大INDEL长度必须>=最小|max_indel_size must be >= min_indel_size")
        if not 0 < self.min_group_consistency <= 1:
            errors.append(f"一致性阈值必须在0-1|min_group_consistency must be in (0,1]: {self.min_group_consistency}")
        if not 0 <= self.max_missing_rate <= 1:
            errors.append(f"缺失率必须在0-1|max_missing_rate must be in [0,1]")
        if self.min_samples_per_group < 1:
            errors.append(f"每组最少样品数必须>=1|min_samples_per_group must be >=1")
        if not 0 < self.deletion_depth_ratio <= 1:
            errors.append(f"骤降阈值必须在0-1|deletion_depth_ratio must be in (0,1]")
        if self.threads <= 0:
            errors.append(f"线程数必须为正|threads must be positive")

        if errors:
            raise ValueError("\n".join(errors))
        return True
