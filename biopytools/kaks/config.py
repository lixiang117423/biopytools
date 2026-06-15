"""
Ka/Ks Calculator配置管理模块|Ka/Ks Calculator Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class KaKsConfig:
    """Ka/Ks分析配置类|Ka/Ks Analysis Configuration Class"""

    # 计算方法参数|Calculation method parameters
    method: str = "GMYN"
    supported_methods: List[str] = field(default_factory=lambda: [
        "GMYN", "MYN", "YN", "NG", "LWL", "LPB",
        "MLWL", "MLPB", "GY", "MS", "MA", "GNG",
        "GLWL", "GLPB", "GMLWL", "GMLPB", "GYN"
    ])

    # 文件扩展名|File extensions
    fasta_extensions: List[str] = field(default_factory=lambda: ['.fasta', '.fa', '.fas', '.fna', '.cds'])
    pair_extensions: List[str] = field(default_factory=lambda: ['.txt', '.tsv', '.csv'])
    axt_extension: str = '.axt'

    # 输出文件名|Output file names
    output_files: Dict[str, str] = field(default_factory=lambda: {
        'summary': 'kaks_summary.xlsx',
        'detailed': 'kaks_detailed.tsv',
        'detailed_csv': 'kaks_detailed.csv',
        'statistics': 'summary_stats.json',
        'report': 'analysis_report.html',
        'log': 'pipeline.log',
        'temp_axt': 'kaks_input.axt'
    })

    # 质量阈值|Quality thresholds
    min_sequence_length: int = 50
    max_sequence_length: int = 50000
    min_similarity_threshold: float = 0.3
    max_missing_ratio: float = 0.1

    # 选择压力分类|Selection pressure classification
    selection_thresholds: Dict[str, float] = field(default_factory=lambda: {
        'strong_negative': 0.5,
        'moderate_negative': 1.0,
        'neutral_lower': 0.95,
        'neutral_upper': 1.05,
        'weak_positive': 2.0
    })

    # 计算参数|Calculation parameters
    calculation_timeout: int = 86400  # 24小时超时|24 hours timeout
    max_batch_size: int = 1000  # 最大批处理大小|Maximum batch size

    # 工具路径|Tool paths
    kaks_path: str = field(
        default_factory=lambda: get_tool_path(
            'KaKs_Calculator',
            '~/miniforge3/envs/kakscalculator/bin/KaKs_Calculator',
            'KAKS_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.kaks_path = expand_path(self.kaks_path)

    def get_method_description(self, method: str) -> str:
        """获取计算方法描述|Get calculation method description"""
        descriptions = {
            "GMYN": "推荐：Gamma分布修正的Yang-Nielsen方法|Recommended: Gamma-corrected Yang-Nielsen",
            "MYN": "修正的Yang-Nielsen方法|Modified Yang-Nielsen method",
            "YN": "经典Yang-Nielsen方法|Classic Yang-Nielsen method",
            "NG": "Nei-Gojobori方法|Nei-Gojobori method",
            "LWL": "Li-Wu-Luo方法|Li-Wu-Luo method",
            "LPB": "Li-Pamilo-Bianchi方法|Li-Pamilo-Bianchi method",
            "MLWL": "修正的LWL方法|Modified LWL method",
            "MLPB": "修正的LPB方法|Modified LPB method",
            "GY": "Goldman-Yang最大似然方法|Goldman-Yang maximum likelihood",
            "MS": "基于AICc的模型选择|Model Selection based on AICc",
            "MA": "候选模型的模型平均|Model Averaging on candidate models",
            "GNG": "Gamma-Nei-Gojobori方法|Gamma-Nei-Gojobori method",
            "GLWL": "Gamma-LWL方法|Gamma-LWL method",
            "GLPB": "Gamma-LPB方法|Gamma-LPB method",
            "GMLWL": "Gamma修正的MLWL方法|Gamma-modified MLWL method",
            "GMLPB": "Gamma修正的MLPB方法|Gamma-modified MLPB method",
            "GYN": "Gamma-Yang-Nielsen方法|Gamma-Yang-Nielsen method"
        }
        return descriptions.get(method, f"{method}计算方法|{method} calculation method")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if self.method not in self.supported_methods:
            errors.append(
                f"不支持的计算方法|Unsupported calculation method: {self.method}. "
                f"支持的方法|Supported methods: {', '.join(self.supported_methods)}"
            )

        if self.min_sequence_length < 1:
            errors.append(f"最小序列长度必须>=1|Minimum sequence length must be >=1")

        if self.max_sequence_length <= self.min_sequence_length:
            errors.append(f"最大序列长度必须>=最小序列长度|Maximum sequence length must be >= minimum")

        if not 0 < self.min_similarity_threshold <= 1:
            errors.append(f"最小相似度阈值必须在0-1之间|Min similarity threshold must be between 0-1")

        if not 0 <= self.max_missing_ratio <= 1:
            errors.append(f"最大缺失率必须在0-1之间|Max missing ratio must be between 0-1")

        if self.calculation_timeout <= 0:
            errors.append(f"计算超时必须为正数|Calculation timeout must be positive")

        if self.max_batch_size <= 0:
            errors.append(f"最大批处理大小必须为正数|Max batch size must be positive")

        if errors:
            raise ValueError("\n".join(errors))

        return True
