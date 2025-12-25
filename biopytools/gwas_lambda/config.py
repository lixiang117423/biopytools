"""
GWAS Lambda计算配置模块 | GWAS Lambda Calculator Configuration Module
"""

import os
from typing import Optional


class GWASLambdaConfig:
    """GWAS Lambda计算配置类 | GWAS Lambda Calculator Configuration Class"""

    def __init__(
        self,
        search_pattern: str = "feture_*/GWAS_Result.mlm.manht_input",
        output_file: str = "Batch_Lambda_Assessment.txt",
        significance_threshold: float = 1e-5,
        p_value_column: int = 3,  # 第4列，索引为3
        min_p_value: float = 1e-300,
        max_p_value: float = 1.0,
        expected_median: float = 0.4549364,  # df=1时的期望中位数
        output_dir: str = "./gwas_lambda_output"
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            search_pattern: 文件搜索模式 | File search pattern
            output_file: 输出文件名 | Output filename
            significance_threshold: 显著性阈值 | Significance threshold
            p_value_column: P值所在列（从0开始）| P-value column index (0-based)
            min_p_value: 最小有效P值 | Minimum valid P-value
            max_p_value: 最大有效P值 | Maximum valid P-value
            expected_median: 期望卡方中位数 | Expected chi-square median
            output_dir: 输出目录 | Output directory
        """
        self.search_pattern = search_pattern
        self.output_file = output_file
        self.significance_threshold = significance_threshold
        self.p_value_column = p_value_column
        self.min_p_value = min_p_value
        self.max_p_value = max_p_value
        self.expected_median = expected_median
        self.output_dir = output_dir

        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def get_output_path(self) -> str:
        """获取输出文件的完整路径 | Get full path of output file"""
        return os.path.join(self.output_dir, self.output_file)

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        if self.significance_threshold <= 0 or self.significance_threshold >= 1:
            raise ValueError("显著性阈值必须在0和1之间 | Significance threshold must be between 0 and 1")

        if self.p_value_column < 0:
            raise ValueError("P值列索引必须非负 | P-value column index must be non-negative")

        if self.min_p_value <= 0 or self.max_p_value > 1:
            raise ValueError("P值范围无效 | Invalid P-value range")

        if self.expected_median <= 0:
            raise ValueError("期望中位数必须为正数 | Expected median must be positive")

    def __str__(self) -> str:
        """配置的字符串表示 | String representation of configuration"""
        return f"""GWASLambdaConfig:
    搜索模式 | Search Pattern: {self.search_pattern}
    输出文件 | Output File: {self.get_output_path()}
    显著性阈值 | Significance Threshold: {self.significance_threshold}
    P值列 | P-value Column: {self.p_value_column}
    P值范围 | P-value Range: [{self.min_p_value}, {self.max_p_value}]
    期望中位数 | Expected Median: {self.expected_median}
"""