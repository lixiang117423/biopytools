"""
GWAS结果分析器模块 | GWAS Result Analyzer Module
"""

import os
import numpy as np
from scipy.stats import chi2
from typing import Tuple, Optional, List, Dict


class GWASResultAnalyzer:
    """GWAS结果分析器类 | GWAS Result Analyzer Class"""

    def __init__(self, significance_threshold: float = 1e-5,
                 p_value_column: int = 3,
                 min_p_value: float = 1e-300,
                 max_p_value: float = 1.0,
                 expected_median: float = 0.4549364):
        """
        初始化分析器 | Initialize analyzer

        Args:
            significance_threshold: 显著性阈值 | Significance threshold
            p_value_column: P值所在列（从0开始）| P-value column index (0-based)
            min_p_value: 最小有效P值 | Minimum valid P-value
            max_p_value: 最大有效P值 | Maximum valid P-value
            expected_median: 期望卡方中位数 | Expected chi-square median
        """
        self.significance_threshold = significance_threshold
        self.p_value_column = p_value_column
        self.min_p_value = min_p_value
        self.max_p_value = max_p_value
        self.expected_median = expected_median

    def analyze_file(self, file_path: str) -> Tuple[Optional[float], int, int, Optional[str]]:
        """
        分析单个GWAS结果文件 | Analyze a single GWAS result file

        Args:
            file_path: 文件路径 | File path

        Returns:
            Tuple: (lambda_gc, total_snps, sig_count, error_msg)
        """
        try:
            p_values = []
            sig_count = 0

            if not os.path.exists(file_path):
                return None, 0, 0, f"文件不存在 | File not found: {file_path}"

            with open(file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split()
                    if len(parts) <= self.p_value_column:
                        continue

                    # 获取P值 | Get P-value
                    val_str = parts[self.p_value_column]

                    try:
                        p = float(val_str)
                        if self.min_p_value < p <= self.max_p_value:
                            p_values.append(p)

                            # 统计显著位点 | Count significant variants
                            if p < self.significance_threshold:
                                sig_count += 1
                    except ValueError:
                        # 跳过表头或无效数据 | Skip headers or invalid data
                        continue

            if not p_values:
                return None, 0, 0, "无有效P值 | No valid P-values found"

            total_snps = len(p_values)

            # 如果没有显著位点，不计算Lambda | If no significant variants, don't calculate Lambda
            if sig_count == 0:
                return None, total_snps, 0, "无显著位点 | No significant variants found"

            # 计算Lambda GC | Calculate Lambda GC
            p_values_arr = np.array(p_values)

            # 转换为卡方值 (1自由度) | Convert to chi-square values (df=1)
            chisq_stats = chi2.ppf(1 - p_values_arr, df=1)

            # 计算中位数 | Calculate median
            obs_median = np.median(chisq_stats)
            lambda_gc = obs_median / self.expected_median

            return lambda_gc, total_snps, sig_count, None

        except Exception as e:
            return None, 0, 0, str(e)

    def get_status(self, lambda_val: Optional[float], sig_count: int) -> str:
        """
        根据Lambda值和显著位点数量判断结果状态 | Determine status based on Lambda and significant count

        Args:
            lambda_val: Lambda GC值 | Lambda GC value
            sig_count: 显著位点数量 | Number of significant variants

        Returns:
            str: 状态描述 | Status description
        """
        # 首先检查是否有显著位点 | First check if there are significant variants
        if sig_count == 0:
            return "No Signals (无显著关联)"

        if lambda_val is None:
            return "Error (计算失败 | Calculation failed)"

        # 判断Lambda质量 | Assess Lambda quality
        if 0.95 <= lambda_val <= 1.05:
            return "Ideal (完美)"
        elif 0.90 <= lambda_val <= 1.10:
            return "Acceptable (可接受)"
        elif lambda_val > 1.10:
            return "Inflated (膨胀/假阳性高)"
        elif 0.80 <= lambda_val < 0.90:
            return "Deflated (轻微过度校正)"
        else:
            return "Warning: Extreme Deflated (异常/模型失效)"

    def should_highlight(self, status: str) -> bool:
        """
        判断是否需要高亮显示 | Determine if highlighting is needed

        Args:
            status: 状态描述 | Status description

        Returns:
            bool: 是否需要高亮 | Whether to highlight
        """
        highlight_keywords = ["Warning", "Inflated", "Extreme", "Error"]
        return any(keyword in status for keyword in highlight_keywords)