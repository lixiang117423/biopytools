"""
Primer3引物设计配置管理模块|Primer3 Primer Design Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple
from ..common.paths import get_tool_path, expand_path


@dataclass
class Primer3Config:
    """Primer3引物设计配置类|Primer3 Primer Design Configuration Class"""

    # 必需参数|Required parameters
    input_fasta: str
    output_dir: str

    # Primer3路径配置|Primer3 path configuration
    primer3_core_path: str = field(
        default_factory=lambda: get_tool_path(
            'primer3_core',
            '~/miniforge3/envs/primer3_v.2.6.1/bin/primer3_core',
            'PRIMER3_PATH'
        )
    )

    # 引物长度参数|Primer size parameters (默认20-22)
    primer_min_size: int = 20
    primer_opt_size: int = 20
    primer_max_size: int = 22

    # 退火温度参数|Annealing temperature parameters (默认58±5, 即53-63)
    primer_min_tm: float = 53.0
    primer_opt_tm: float = 58.0
    primer_max_tm: float = 63.0

    # 产物大小范围|Product size range
    primer_product_size_range: Tuple[int, int] = (100, 300)

    # 其他设计参数|Other design parameters
    primer_num_return: int = 5
    primer_max_ns_accepted: int = 0
    primer_gc_clamp: int = 1

    # 输出格式|Output format
    output_format: str = 'csv'  # csv, tsv, or xlsx
    output_header_lang: str = 'zh'  # zh (中文) or en (英文)

    # 引物设计策略|Primer design strategy
    method: str = 'all'  # all: 覆盖头尾(Cover entire sequence), random: 随机设计(Random design)
    primer_end_margin: int = 200  # 两端允许的引物位置范围|Allowed margin at ends (only for method=all)

    # 引物位置限制|Primer position constraints (内部使用，由method自动设置|Internal use, auto-set by method)
    primer_force_ends: bool = field(init=False)  # 是否强制引物在序列两端|Force primers at sequence ends

    # 自动产物大小范围|Automatic product size range
    auto_product_size: bool = True  # 自动根据序列长度设置产物大小范围|Auto set product size range based on sequence length
    product_size_min_ratio: float = 0.5  # 产物最小长度占序列长度的比例|Min product size ratio to sequence length
    product_size_max_ratio: float = 1.0  # 产物最大长度占序列长度的比例|Max product size ratio to sequence length

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 根据method自动设置primer_force_ends|Auto-set primer_force_ends based on method
        if self.method == 'all':
            self.primer_force_ends = True
        elif self.method == 'random':
            self.primer_force_ends = False
        else:
            raise ValueError(f"不支持的引物设计策略|Unsupported method: {self.method}, 支持|supported: 'all', 'random'")

        # 展开路径|Expand paths
        self.input_fasta = expand_path(self.input_fasta)
        self.output_dir = expand_path(self.output_dir)
        self.primer3_core_path = expand_path(self.primer3_core_path)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_fasta):
            errors.append(f"输入FASTA文件不存在|Input FASTA file not found: {self.input_fasta}")

        # 检查primer3_core|Check primer3_core
        if not os.path.exists(self.primer3_core_path):
            errors.append(f"primer3_core不存在|primer3_core not found: {self.primer3_core_path}")

        # 验证引物长度参数|Validate primer size parameters
        if not (18 <= self.primer_min_size <= self.primer_opt_size <= self.primer_max_size <= 40):
            errors.append(
                f"引物长度参数无效|Invalid primer size parameters: "
                f"min({self.primer_min_size}) <= opt({self.primer_opt_size}) <= max({self.primer_max_size}), "
                f"范围应为18-40|range should be 18-40"
            )

        # 验证退火温度参数|Validate temperature parameters
        if not (30.0 <= self.primer_min_tm <= self.primer_opt_tm <= self.primer_max_tm <= 90.0):
            errors.append(
                f"退火温度参数无效|Invalid temperature parameters: "
                f"min({self.primer_min_tm}) <= opt({self.primer_opt_tm}) <= max({self.primer_max_tm}), "
                f"范围应为30.0-90.0°C|range should be 30.0-90.0°C"
            )

        # 验证产物大小范围|Validate product size range
        if self.primer_product_size_range[0] >= self.primer_product_size_range[1]:
            errors.append(
                f"产物大小范围无效|Invalid product size range: "
                f"{self.primer_product_size_range[0]} >= {self.primer_product_size_range[1]}"
            )

        # 验证method参数|Validate method parameter
        valid_methods = ['all', 'random']
        if self.method not in valid_methods:
            errors.append(f"引物设计策略无效|Invalid method: {self.method}, 支持|supported: {valid_methods}")

        # 验证输出格式|Validate output format
        valid_formats = ['csv', 'tsv', 'xlsx']
        if self.output_format not in valid_formats:
            errors.append(f"输出格式无效|Invalid output format: {self.output_format}, 支持|supported: {valid_formats}")

        # 验证表头语言|Validate header language
        valid_langs = ['zh', 'en']
        if self.output_header_lang not in valid_langs:
            errors.append(f"表头语言无效|Invalid header language: {self.output_header_lang}, 支持|supported: {valid_langs}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_primer3_settings(self) -> dict:
        """
        获取Primer3设置字典|Get Primer3 settings dictionary

        Returns:
            dict: Primer3参数字典|Primer3 parameters dictionary
        """
        settings = {
            'PRIMER_TASK': 'generic',
            'PRIMER_MIN_SIZE': self.primer_min_size,
            'PRIMER_OPT_SIZE': self.primer_opt_size,
            'PRIMER_MAX_SIZE': self.primer_max_size,
            'PRIMER_MIN_TM': self.primer_min_tm,
            'PRIMER_OPT_TM': self.primer_opt_tm,
            'PRIMER_MAX_TM': self.primer_max_tm,
            'PRIMER_PRODUCT_SIZE_RANGE': f"{self.primer_product_size_range[0]}-{self.primer_product_size_range[1]}",
            'PRIMER_NUM_RETURN': self.primer_num_return,
            'PRIMER_MAX_NS_ACCEPTED': self.primer_max_ns_accepted,
            'PRIMER_GC_CLAMP': self.primer_gc_clamp,
        }

        # 添加引物位置限制|Add primer position constraints
        if self.primer_force_ends:
            # 从序列起始位置开始设计引物|Start from sequence beginning
            settings['PRIMER_FIRST_PRODUCT_INTERVAL'] = 1

        return settings

    def get_sequence_specific_settings(self, seq_length: int) -> dict:
        """
        获取序列特定的Primer3设置|Get sequence-specific Primer3 settings

        Args:
            seq_length: 序列长度|Sequence length

        Returns:
            dict: 序列特定的Primer3参数字典|Sequence-specific Primer3 parameters dictionary
        """
        settings = {}

        # 自动设置产物大小范围|Auto set product size range
        if self.auto_product_size:
            # 计算产物大小范围|Calculate product size range
            min_product = max(
                int(seq_length * self.product_size_min_ratio),
                self.primer_product_size_range[0]  # 不低于全局最小值
            )

            # 当method=all时（primer_force_ends=True），确保产物能覆盖完整序列
            # When method=all (primer_force_ends=True), ensure product can cover full sequence
            if self.primer_force_ends:
                # method=all模式：产物最大值至少为序列长度，不受全局max限制
                # method=all mode: max product at least sequence length, ignore global max limit
                max_product = int(seq_length * self.product_size_max_ratio)
            else:
                # method=random模式：受全局max限制
                # method=random mode: respect global max limit
                max_product = min(
                    int(seq_length * self.product_size_max_ratio),
                    self.primer_product_size_range[1]
                )

            # 确保最小值不大于最大值|Ensure min <= max
            min_product = min(min_product, max_product)

            settings['PRIMER_PRODUCT_SIZE_RANGE'] = f"{min_product}-{max_product}"

        if self.primer_force_ends:
            # 确保引物设计在序列两端|Ensure primers are designed at sequence ends
            # Primer3使用0-based索引|Primer3 uses 0-based indexing
            # 正向引物区域：起始端|Forward primer region: start end
            left_start = 0
            left_length = min(self.primer_end_margin, seq_length // 2)

            # 反向引物区域：终止端|Reverse primer region: end end
            # right_start需要从序列末尾向前计算|right_start calculated from end backwards
            right_start = max(seq_length - self.primer_end_margin, left_start + left_length)
            right_length = seq_length - right_start

            # 设置允许的引物对区域|Set allowed primer pair region
            settings['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = f"{left_start},{left_length},{right_start},{right_length}"

        return settings
