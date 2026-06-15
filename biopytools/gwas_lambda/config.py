"""
GWAS Lambda计算配置模块|GWAS Lambda Calculator Configuration Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class GWASLambdaConfig:
    """GWAS Lambda计算配置类|GWAS Lambda Calculator Configuration Class"""

    # 必需参数|Required parameters
    search_pattern: str = "feture_*/GWAS_Result.mlm.manht_input"

    # 可选参数|Optional parameters
    output_file: str = "Batch_Lambda_Assessment.txt"
    significance_threshold: float = 1e-5
    p_value_column: int = 3  # 第4列，索引为3|Column 4, index 3
    min_p_value: float = 1e-300
    max_p_value: float = 1.0
    expected_median: float = 0.4549364  # df=1时的期望中位数|Expected median when df=1
    output_dir: str = "./gwas_lambda_output"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 确保输出目录存在|Ensure output directory exists
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def get_output_path(self) -> str:
        """获取输出文件的完整路径|Get full path of output file"""
        return os.path.join(self.output_dir, self.output_file)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if self.significance_threshold <= 0 or self.significance_threshold >= 1:
            errors.append("显著性阈值必须在0和1之间|Significance threshold must be between 0 and 1")

        if self.p_value_column < 0:
            errors.append("P值列索引必须非负|P-value column index must be non-negative")

        if self.min_p_value <= 0 or self.max_p_value > 1:
            errors.append("P值范围无效|Invalid P-value range")

        if self.expected_median <= 0:
            errors.append("期望中位数必须为正数|Expected median must be positive")

        if errors:
            raise ValueError("\n".join(errors))

    def __str__(self) -> str:
        """配置的字符串表示|String representation of configuration"""
        return f"""GWASLambdaConfig:
    搜索模式|Search Pattern: {self.search_pattern}
    输出文件|Output File: {self.get_output_path()}
    显著性阈值|Significance Threshold: {self.significance_threshold}
    P值列|P-value Column: {self.p_value_column}
    P值范围|P-value Range: [{self.min_p_value}, {self.max_p_value}]
    期望中位数|Expected Median: {self.expected_median}
"""
