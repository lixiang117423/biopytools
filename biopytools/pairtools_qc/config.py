"""
Hi-C数据质量控制评估配置管理模块|Hi-C QC Assessment Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict
from ..common.paths import expand_path


@dataclass
class PairtoolsQCConfig:
    """Hi-C数据质量控制评估配置类|Hi-C QC Assessment Configuration Class"""

    # 必需参数|Required parameters
    pairs_file: str

    # 可选参数|Optional parameters
    pairtools_path: str = '~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools'
    output_dir: str = './pairtools_qc_output'
    chroms_path: str = None  # 染色体大小文件路径（BAM输入时必需）|Chromosome sizes file path (required for BAM input)

    # 质量阈值|Quality thresholds
    max_unmapped_rate: float = 20.0  # 未比对reads比例阈值 (%) | Threshold for unmapped reads rate (%)
    max_single_sided_rate: float = 10.0  # 单端比对比例阈值 (%) | Threshold for single-sided mapping rate (%)
    min_mapped_rate: float = 80.0  # 双端比对率阈值 (%) | Threshold for paired mapping rate (%)
    max_dup_rate: float = 30.0  # PCR重复率阈值 (%) | Threshold for PCR duplication rate (%)
    min_cis_trans_ratio: float = 4.0  # cis/trans比例阈值 | Threshold for cis/trans ratio

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.pairs_file = os.path.normpath(os.path.abspath(self.pairs_file))
        self.pairtools_path = os.path.normpath(os.path.abspath(expand_path(self.pairtools_path)))

        # 标准化chroms_path（如果提供）|Normalize chroms_path if provided
        if self.chroms_path:
            self.chroms_path = os.path.normpath(os.path.abspath(self.chroms_path))

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件是否存在|Check if input file exists
        if not os.path.exists(self.pairs_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.pairs_file}")

        # 检查是否为BAM文件，如果是则必须有chroms_path|Check if BAM file, requires chroms_path
        file_ext = self.pairs_file.lower()
        if (file_ext.endswith('.bam') or file_ext.endswith('.sam')):
            if not self.chroms_path:
                errors.append(
                    f"处理BAM/SAM文件需要提供chrom.sizes文件|"
                    f"Processing BAM/SAM files requires chrom.sizes file\n"
                    f"请使用 --chroms-path 参数指定|"
                    f"Please specify with --chroms-path parameter"
                )
            elif not os.path.exists(self.chroms_path):
                errors.append(f"Chrom.sizes文件不存在|Chrom.sizes file does not exist: {self.chroms_path}")

        # 检查pairtools可执行文件是否存在|Check if pairtools executable exists
        if not os.path.exists(self.pairtools_path):
            errors.append(f"Pairtools可执行文件不存在|Pairtools executable does not exist: {self.pairtools_path}")

        # 检查阈值范围|Check threshold ranges
        if not (0 <= self.max_unmapped_rate <= 100):
            errors.append(f"未比对阈值范围错误|Invalid unmapped threshold: {self.max_unmapped_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.max_single_sided_rate <= 100):
            errors.append(f"单端比对阈值范围错误|Invalid single-sided threshold: {self.max_single_sided_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.min_mapped_rate <= 100):
            errors.append(f"双端比对阈值范围错误|Invalid mapped threshold: {self.min_mapped_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.max_dup_rate <= 100):
            errors.append(f"PCR重复阈值范围错误|Invalid duplication threshold: {self.max_dup_rate} (应为0-100|should be 0-100)")

        if self.min_cis_trans_ratio < 0:
            errors.append(f"cis/trans阈值不能为负数|cis/trans threshold cannot be negative: {self.min_cis_trans_ratio}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_thresholds(self) -> Dict[str, float]:
        """获取所有阈值|Get all thresholds

        Returns:
            Dict[str, float]: 阈值字典|Threshold dictionary
        """
        return {
            'max_unmapped_rate': self.max_unmapped_rate,
            'max_single_sided_rate': self.max_single_sided_rate,
            'min_mapped_rate': self.min_mapped_rate,
            'max_dup_rate': self.max_dup_rate,
            'min_cis_trans_ratio': self.min_cis_trans_ratio
        }
