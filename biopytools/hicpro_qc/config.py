"""
HiC-Pro质量控制评估配置管理模块|HiC-Pro QC Assessment Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional


@dataclass
class HiCProQCConfig:
    """HiC-Pro质量控制评估配置类|HiC-Pro QC Assessment Configuration Class"""

    # 必需参数|Required parameters
    hicpro_dir: str  # HiC-Pro输出目录|HiC-Pro output directory

    # 可选参数|Optional parameters
    output_dir: str = './hicpro_qc_output'
    sample_name: Optional[str] = None  # 样本名称（可选，自动检测）|Sample name (optional, auto-detect)

    # 质量阈值|Quality thresholds
    # Mapping统计阈值|Mapping thresholds
    min_mapping_rate: float = 70.0  # 最低比对率 (%) | Minimum mapping rate (%)
    min_unique_rate: float = 60.0  # 最低唯一比对率 (%) | Minimum unique mapping rate (%)

    # Valid pairs阈值|Valid pairs thresholds
    min_valid_pairs_rate: float = 50.0  # 最低valid pairs比例 (%) | Minimum valid pairs rate (%)
    max_dangling_ends_rate: float = 15.0  # 最高dangling ends比例 (%) | Maximum dangling ends rate (%)
    max_self_ligation_rate: float = 5.0  # 最高self-ligation比例 (%) | Maximum self-ligation rate (%)
    max_religation_rate: float = 10.0  # 最高religation比例 (%) | Maximum religation rate (%)

    # Cis/Trans阈值|Cis/Trans thresholds
    min_cis_trans_ratio: float = 5.0  # 最低cis/trans比例 | Minimum cis/trans ratio

    # PCR重复阈值|PCR duplication threshold
    max_duplication_rate: float = 30.0  # 最高PCR重复率 (%) | Maximum PCR duplication rate (%)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.hicpro_dir = os.path.normpath(os.path.abspath(self.hicpro_dir))

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 自动检测样本名称（如果未提供）|Auto-detect sample name if not provided
        if self.sample_name is None:
            self.sample_name = self._detect_sample_name()

    def _detect_sample_name(self) -> str:
        """自动检测样本名称|Auto-detect sample name from HiC-Pro output

        Returns:
            str: 检测到的样本名称|Detected sample name
        """
        # 方法1：查找.mmapstat文件（合并的mapping统计）
        bowtie_dir = Path(self.hicpro_dir) / 'bowtie_results' / 'bwt2'
        if bowtie_dir.exists():
            mmapstat_files = list(bowtie_dir.glob('*.mmapstat'))
            if mmapstat_files:
                # 提取样本名（去掉.mmapstat后缀）
                sample = mmapstat_files[0].stem
                return sample

        # 方法2：查找.mpairstat文件（合并的pairing统计）
        if bowtie_dir.exists():
            mpairstat_files = list(bowtie_dir.glob('*.mpairstat'))
            if mpairstat_files:
                sample = mpairstat_files[0].stem
                return sample

        # 方法3：查找.mRSstat文件（合并的RS统计）
        hic_results_dir = Path(self.hicpro_dir) / 'hic_results' / 'data'
        if hic_results_dir.exists():
            mrsstat_files = list(hic_results_dir.glob('*.mRSstat'))
            if mrsstat_files:
                sample = mrsstat_files[0].stem
                return sample

        # 如果都找不到，返回默认值
        return "sample"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查HiC-Pro目录是否存在|Check if HiC-Pro directory exists
        if not os.path.exists(self.hicpro_dir):
            errors.append(f"HiC-Pro输出目录不存在|HiC-Pro output directory not found: {self.hicpro_dir}")
        else:
            # 检查必要的子目录|Check required subdirectories
            required_dirs = ['bowtie_results', 'hic_results']
            for req_dir in required_dirs:
                dir_path = Path(self.hicpro_dir) / req_dir
                if not dir_path.exists():
                    errors.append(f"HiC-Pro缺少必要子目录|Missing required HiC-Pro subdirectory: {req_dir}")

        # 检查阈值范围|Check threshold ranges
        if not (0 <= self.min_mapping_rate <= 100):
            errors.append(f"比对率阈值范围错误|Invalid mapping rate threshold: {self.min_mapping_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.min_unique_rate <= 100):
            errors.append(f"唯一比对率阈值范围错误|Invalid unique rate threshold: {self.min_unique_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.min_valid_pairs_rate <= 100):
            errors.append(f"Valid pairs阈值范围错误|Invalid valid pairs threshold: {self.min_valid_pairs_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.max_dangling_ends_rate <= 100):
            errors.append(f"Dangling ends阈值范围错误|Invalid dangling ends threshold: {self.max_dangling_ends_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.max_self_ligation_rate <= 100):
            errors.append(f"Self-ligation阈值范围错误|Invalid self-ligation threshold: {self.max_self_ligation_rate} (应为0-100|should be 0-100)")

        if not (0 <= self.max_religation_rate <= 100):
            errors.append(f"Religation阈值范围错误|Invalid religation threshold: {self.max_religation_rate} (应为0-100|should be 0-100)")

        if self.min_cis_trans_ratio < 0:
            errors.append(f"cis/trans阈值不能为负数|cis/trans threshold cannot be negative: {self.min_cis_trans_ratio}")

        if not (0 <= self.max_duplication_rate <= 100):
            errors.append(f"PCR重复率阈值范围错误|Invalid duplication threshold: {self.max_duplication_rate} (应为0-100|should be 0-100)")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_thresholds(self) -> Dict[str, float]:
        """获取所有阈值|Get all thresholds

        Returns:
            Dict[str, float]: 阈值字典|Threshold dictionary
        """
        return {
            'min_mapping_rate': self.min_mapping_rate,
            'min_unique_rate': self.min_unique_rate,
            'min_valid_pairs_rate': self.min_valid_pairs_rate,
            'max_dangling_ends_rate': self.max_dangling_ends_rate,
            'max_self_ligation_rate': self.max_self_ligation_rate,
            'max_religation_rate': self.max_religation_rate,
            'min_cis_trans_ratio': self.min_cis_trans_ratio,
            'max_duplication_rate': self.max_duplication_rate
        }
