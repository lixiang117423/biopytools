"""
甲基化分析配置管理模块 | Methylation Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass
class MethylationConfig:
    """甲基化分析配置类 | Methylation Analysis Configuration Class"""

    # 必需输入文件 | Required input files
    raw_dir: str
    genome_fa: str
    output_dir: str  # <--- 已修改: 移除默认值，使其成为必需参数

    # 可选输入文件 | Optional input files
    target_fa: Optional[str] = None
    target_promoter_fa: Optional[str] = None
    annotation_gff: Optional[str] = None

    # 分析模式 | Analysis mode
    enhanced_mode: bool = False  # 启用增强分析 | Enable enhanced analysis

    # 样品分组 | Sample grouping
    sample_groups: Dict[str, str] = field(
        default_factory=lambda: {
            "FZY4201": "CK",
            "FZY4203": "CK",
            "FZY4205": "CK",
            "FZY4207": "CK",
            "FZY4202": "Treatment",
            "FZY4204": "Treatment",
            "FZY4206": "Treatment",
            "FZY4208": "Treatment",
        }
    )

    # 基础分析参数 | Basic analysis parameters
    threads: int = 8
    min_coverage: int = 5
    min_cytosines: int = 5

    # 差异分析参数 | Differential analysis parameters
    methylation_diff_threshold: float = 0.3
    pvalue_threshold: float = 0.01
    window_size: int = 100
    step_size: int = 50
    merge_distance: int = 200

    # bins分析参数 | Bins analysis parameters
    gene_bins: int = 60
    flank_bins: int = 20
    flank_size: int = 2000

    # 基础工具路径 | Basic tool paths
    fastp_path: str = "fastp"
    bismark_path: str = "bismark"
    bismark_genome_preparation_path: str = "bismark_genome_preparation"
    bowtie2_path: str = "bowtie2"
    deduplicate_bismark_path: str = "deduplicate_bismark"
    bismark_methylation_extractor_path: str = "bismark_methylation_extractor"
    bismark2report_path: str = "bismark2report"
    bismark2summary_path: str = "bismark2summary"

    # 可选工具路径 | Optional tool paths
    makeblastdb_path: str = "makeblastdb"
    blastn_path: str = "blastn"
    multiqc_path: str = "multiqc"
    bedtools_path: str = "bedtools"
    samtools_path: str = "samtools"

    # R环境路径 | R environment path
    r_executable: str = "/share/org/YZWL/yzwl_lixg/miniforge3/envs/methylkit/bin/R"

    # 内部属性 | Internal attributes
    clean_dir: str = field(init=False)
    mapping_dir: str = field(init=False)
    analysis_dir: str = field(init=False)

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.raw_dir = os.path.normpath(os.path.abspath(self.raw_dir))
        self.genome_fa = os.path.normpath(os.path.abspath(self.genome_fa))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 设置子目录 | Set subdirectories
        self.clean_dir = os.path.join(self.output_dir, "clean_data")
        self.mapping_dir = os.path.join(self.output_dir, "mapping")
        self.analysis_dir = os.path.join(self.output_dir, "analysis")  # 通用分析目录

        if self.enhanced_mode:
            self.analysis_dir = os.path.join(self.output_dir, "methylkit")

        # 标准化可选文件路径 | Normalize optional file paths
        if self.target_fa:
            self.target_fa = os.path.normpath(os.path.abspath(self.target_fa))
        if self.target_promoter_fa:
            self.target_promoter_fa = os.path.normpath(
                os.path.abspath(self.target_promoter_fa)
            )
        if self.annotation_gff:
            self.annotation_gff = os.path.normpath(os.path.abspath(self.annotation_gff))

    def validate(self):
        """验证配置 | Validate configuration"""
        # 检查必需文件 | Check required files
        if not os.path.exists(self.raw_dir):
            raise FileNotFoundError(
                f"原始数据目录不存在 | Raw data directory not found: {self.raw_dir}"
            )

        if not os.path.exists(self.genome_fa):
            raise FileNotFoundError(
                f"基因组文件不存在 | Genome file not found: {self.genome_fa}"
            )

        # 检查可选文件 | Check optional files
        if self.target_fa and not os.path.exists(self.target_fa):
            raise FileNotFoundError(
                f"目标序列文件不存在 | Target sequence file not found: {self.target_fa}"
            )

        if self.target_promoter_fa and not os.path.exists(self.target_promoter_fa):
            raise FileNotFoundError(
                f"目标启动子序列文件不存在 | Target promoter file not found: {self.target_promoter_fa}"
            )

        if self.annotation_gff and not os.path.exists(self.annotation_gff):
            raise FileNotFoundError(
                f"基因注释文件不存在 | Annotation file not found: {self.annotation_gff}"
            )

        # 检查增强模式所需条件 | Check enhanced mode requirements
        if self.enhanced_mode and not os.path.exists(self.r_executable):
            raise FileNotFoundError(
                f"R可执行文件不存在 | R executable not found: {self.r_executable}"
            )

        # 验证数值参数合理性 | Validate numeric parameter ranges
        if self.threads <= 0:
            raise ValueError("线程数必须大于0 | Thread count must be greater than 0")

        if self.min_coverage <= 0:
            raise ValueError(
                "最小覆盖度必须大于0 | Minimum coverage must be greater than 0"
            )

        if self.min_cytosines <= 0:
            raise ValueError(
                "最小胞嘧啶数必须大于0 | Minimum cytosines must be greater than 0"
            )

        if not (0 < self.methylation_diff_threshold < 1):
            raise ValueError(
                "甲基化差异阈值必须在0-1之间 | Methylation difference threshold must be between 0-1"
            )

        if not (0 < self.pvalue_threshold < 1):
            raise ValueError(
                "p值阈值必须在0-1之间 | P-value threshold must be between 0-1"
            )

        if self.window_size <= 0:
            raise ValueError("窗口大小必须大于0 | Window size must be greater than 0")

        if self.step_size <= 0:
            raise ValueError("步长必须大于0 | Step size must be greater than 0")

        if self.merge_distance < 0:
            raise ValueError("合并距离不能为负数 | Merge distance cannot be negative")

        if self.gene_bins <= 0:
            raise ValueError(
                "基因体bins数必须大于0 | Gene bins count must be greater than 0"
            )

        if self.flank_bins <= 0:
            raise ValueError(
                "侧翼bins数必须大于0 | Flank bins count must be greater than 0"
            )

        if self.flank_size <= 0:
            raise ValueError(
                "侧翼区域大小必须大于0 | Flank size must be greater than 0"
            )

        # 验证样品分组 | Validate sample grouping
        if not self.sample_groups:
            raise ValueError("样品分组不能为空 | Sample groups cannot be empty")

        # 检查分组中是否有CK和Treatment组 | Check if both CK and Treatment groups exist
        groups = set(self.sample_groups.values())
        if len(groups) < 2:
            raise ValueError(
                "至少需要两个不同的分组 | At least two different groups are required"
            )

        # 警告：如果没有CK或Treatment组 | Warning: if no CK or Treatment group
        if "CK" not in groups:
            print("⚠️  警告: 样品分组中没有CK组，可能影响分析结果")
        if "Treatment" not in groups:
            print("⚠️  警告: 样品分组中没有Treatment组，可能影响分析结果")

    @property
    def ck_samples(self) -> List[str]:
        """获取CK组样品列表 | Get CK group samples"""
        return [sample for sample, group in self.sample_groups.items() if group == "CK"]

    @property
    def treatment_samples(self) -> List[str]:
        """获取Treatment组样品列表 | Get Treatment group samples"""
        return [
            sample
            for sample, group in self.sample_groups.items()
            if group == "Treatment"
        ]

    @property
    def all_groups(self) -> List[str]:
        """获取所有分组列表 | Get all group names"""
        return list(set(self.sample_groups.values()))

    def get_samples_by_group(self, group_name: str) -> List[str]:
        """根据分组名获取样品列表 | Get samples by group name"""
        return [
            sample
            for sample, group in self.sample_groups.items()
            if group == group_name
        ]

    def get_pairwise_comparisons(self) -> List[Tuple[str, str]]:
        """获取所有可能的分组两两比较组合 | Get all possible pairwise group comparisons"""
        groups = self.all_groups
        if len(groups) < 2:
            return []

        # 生成所有两两组合，避免重复（A vs B 和 B vs A是一样的）
        comparisons = list(combinations(groups, 2))
        return comparisons
