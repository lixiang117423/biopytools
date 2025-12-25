"""
滑动窗口分析模块 | Sliding Window Analysis Module
用于生成类似Figure 3B的折线图 | Used to generate line plots similar to Figure 3B
"""

import numpy as np
import pandas as pd
import logging
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from .config import SNPIndexConfig


class SlidingWindowAnalyzer:
    """滑动窗口分析器 | Sliding Window Analyzer"""

    def __init__(self, data: List[Dict], config: SNPIndexConfig):
        """
        初始化滑动窗口分析器 | Initialize sliding window analyzer

        Args:
            data: SNP index数据 | SNP index data
            config: 配置对象 | Configuration object
        """
        self.data = data
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self.windows = []
        self.chromosome_lengths = {}

        # 默认滑动窗口参数 | Default sliding window parameters
        self.window_size = config.window_size if hasattr(config, 'window_size') else 1000000  # 1Mb
        self.step_size = config.step_size if hasattr(config, 'step_size') else 100000        # 100kb
        self.min_snps = config.min_window_snps if hasattr(config, 'min_window_snps') else 5

    def calculate_chromosome_lengths(self) -> None:
        """计算染色体长度 | Calculate chromosome lengths"""
        chrom_max_positions = defaultdict(int)
        for record in self.data:
            chrom_max_positions[record['chromosome']] = max(
                chrom_max_positions[record['chromosome']],
                record['position']
            )

        self.chromosome_lengths = dict(chrom_max_positions)
        self.logger.info(f"计算得到染色体长度 | Calculated chromosome lengths: {self.chromosome_lengths}")

    def create_sliding_windows(self) -> List[Dict]:
        """
        创建滑动窗口 | Create sliding windows

        Returns:
            list: 窗口数据列表 | Window data list
        """
        self.logger.info("创建滑动窗口 | Creating sliding windows")
        self.logger.info(f"窗口大小 | Window size: {self.window_size/1000:.0f} kb")
        self.logger.info(f"步长 | Step size: {self.step_size/1000:.0f} kb")
        self.logger.info(f"最少SNP数 | Minimum SNPs: {self.min_snps}")

        if not self.chromosome_lengths:
            self.calculate_chromosome_lengths()

        windows = []

        for chrom, length in self.chromosome_lengths.items():
            chrom_windows = self._create_chromosome_windows(chrom, length)
            windows.extend(chrom_windows)
            self.logger.info(f"染色体 {chrom} 创建了 {len(chrom_windows)} 个窗口 | "
                           f"Chromosome {chrom}: {len(chrom_windows)} windows created")

        self.windows = windows
        return windows

    def _create_chromosome_windows(self, chrom: str, length: int) -> List[Dict]:
        """
        为单个染色体创建窗口 | Create windows for single chromosome

        Args:
            chrom: 染色体名 | Chromosome name
            length: 染色体长度 | Chromosome length

        Returns:
            list: 窗口列表 | Window list
        """
        windows = []

        # 获取该染色体的SNP数据 | Get SNP data for this chromosome
        chrom_snps = [snp for snp in self.data if snp['chromosome'] == chrom]

        # 按位置排序 | Sort by position
        chrom_snps.sort(key=lambda x: x['position'])

        # 创建窗口 | Create windows
        start = 0
        window_id = 0

        while start < length:
            end = min(start + self.window_size, length)
            center = (start + end) / 2

            # 查找窗口内的SNP | Find SNPs in window
            window_snps = [
                snp for snp in chrom_snps
                if start <= snp['position'] < end
            ]

            if len(window_snps) >= self.min_snps:
                # 计算窗口统计量 | Calculate window statistics
                delta_indices = [snp['delta_snp_index'] for snp in window_snps]

                window_data = {
                    'window_id': window_id,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'center': center,
                    'snp_count': len(window_snps),
                    'mean_delta': np.mean(delta_indices),
                    'median_delta': np.median(delta_indices),
                    'std_delta': np.std(delta_indices),
                    'min_delta': np.min(delta_indices),
                    'max_delta': np.max(delta_indices),
                    'positions': [snp['position'] for snp in window_snps],
                    'delta_values': delta_indices
                }
                windows.append(window_data)
                window_id += 1

            start += self.step_size

        return windows

    def calculate_confidence_intervals(self, confidence_level: float = 0.95) -> Dict:
        """
        计算置信区间 | Calculate confidence intervals

        Args:
            confidence_level: 置信水平 | Confidence level

        Returns:
            dict: 置信区间数据 | Confidence interval data
        """
        self.logger.info(f"计算 {confidence_level*100:.0f}% 置信区间 | "
                        f"Calculating {confidence_level*100:.0f}% confidence intervals")

        if not self.windows:
            self.create_sliding_windows()

        # 合并所有窗口的数据 | Combine all window data
        all_deltas = [w['mean_delta'] for w in self.windows]

        if len(all_deltas) < 2:
            self.logger.warning("数据点太少，无法计算置信区间 | Too few data points for confidence intervals")
            return {}

        # 计算置信区间 | Calculate confidence intervals
        alpha = 1 - confidence_level
        from scipy import stats

        mean_delta = np.mean(all_deltas)
        std_delta = np.std(all_deltas)
        n = len(all_deltas)

        # 使用t分布计算置信区间 | Use t-distribution for confidence intervals
        t_value = stats.t.ppf(1 - alpha/2, n - 1)
        margin_error = t_value * std_delta / np.sqrt(n)

        ci_lower = mean_delta - margin_error
        ci_upper = mean_delta + margin_error

        # 计算基于百分位的置信区间 | Calculate percentile-based confidence intervals
        ci_lower_percentile = np.percentile(all_deltas, (1 - confidence_level) * 100 / 2)
        ci_upper_percentile = np.percentile(all_deltas, 100 - (1 - confidence_level) * 100 / 2)

        confidence_intervals = {
            'mean': mean_delta,
            'std': std_delta,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'ci_lower_percentile': ci_lower_percentile,
            'ci_upper_percentile': ci_upper_percentile,
            'confidence_level': confidence_level,
            'margin_error': margin_error
        }

        self.logger.info(f"置信区间计算完成 | Confidence intervals calculated: "
                        f"[{ci_lower:.4f}, {ci_upper:.4f}]")

        return confidence_intervals

    def identify_candidate_regions(self, threshold: float = None) -> List[Dict]:
        """
        识别候选区域 | Identify candidate regions

        Args:
            threshold: 显著性阈值 | Significance threshold

        Returns:
            list: 候选区域列表 | Candidate region list
        """
        if threshold is None:
            threshold = self.config.region_threshold

        self.logger.info(f"识别候选区域 (阈值 | threshold: {threshold}) | "
                        f"Identifying candidate regions (threshold: {threshold})")

        if not self.windows:
            self.create_sliding_windows()

        # 找到超过阈值的窗口 | Find windows exceeding threshold
        significant_windows = [
            w for w in self.windows
            if abs(w['mean_delta']) > threshold
        ]

        # 合并连续的窗口 | Merge consecutive windows
        candidate_regions = self._merge_consecutive_windows(significant_windows)

        self.logger.info(f"识别到 {len(candidate_regions)} 个候选区域 | "
                        f"Identified {len(candidate_regions)} candidate regions")

        return candidate_regions

    def _merge_consecutive_windows(self, windows: List[Dict]) -> List[Dict]:
        """
        合并连续的窗口 | Merge consecutive windows

        Args:
            windows: 窗口列表 | Window list

        Returns:
            list: 合并后的区域列表 | Merged region list
        """
        if not windows:
            return []

        # 按染色体和位置排序 | Sort by chromosome and position
        windows.sort(key=lambda x: (x['chromosome'], x['center']))

        regions = []
        current_region = None

        for window in windows:
            if (current_region is None or
                window['chromosome'] != current_region['chromosome'] or
                window['start'] > current_region['end']):

                # 开始新区域 | Start new region
                if current_region:
                    regions.append(current_region)

                current_region = {
                    'chromosome': window['chromosome'],
                    'start': window['start'],
                    'end': window['end'],
                    'window_count': 1,
                    'mean_delta': window['mean_delta'],
                    'max_delta': abs(window['mean_delta']),
                    'windows': [window]
                }
            else:
                # 扩展当前区域 | Extend current region
                current_region['end'] = window['end']
                current_region['window_count'] += 1
                current_region['windows'].append(window)

                # 更新平均ΔSNP index | Update mean ΔSNP index
                all_deltas = [w['mean_delta'] for w in current_region['windows']]
                current_region['mean_delta'] = np.mean(all_deltas)
                current_region['max_delta'] = max(abs(w['mean_delta']) for w in current_region['windows'])

        # 添加最后一个区域 | Add last region
        if current_region:
            regions.append(current_region)

        return regions

    def export_windows(self, output_file: str) -> bool:
        """
        导出窗口数据 | Export window data

        Args:
            output_file: 输出文件路径 | Output file path

        Returns:
            bool: 导出是否成功 | Whether export succeeded
        """
        try:
            if not self.windows:
                self.create_sliding_windows()

            # 准备数据 | Prepare data
            export_data = []
            for window in self.windows:
                export_data.append({
                    'Chromosome': window['chromosome'],
                    'Start': window['start'],
                    'End': window['end'],
                    'Center': window['center'],
                    'SNP_Count': window['snp_count'],
                    'Mean_Delta_SNP_Index': window['mean_delta'],
                    'Median_Delta_SNP_Index': window['median_delta'],
                    'Std_Delta_SNP_Index': window['std_delta'],
                    'Min_Delta_SNP_Index': window['min_delta'],
                    'Max_Delta_SNP_Index': window['max_delta']
                })

            # 写入文件 | Write to file
            df = pd.DataFrame(export_data)
            df.to_csv(output_file, sep='\t', index=False)

            self.logger.info(f"窗口数据已导出 | Window data exported: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"导出窗口数据时出错 | Error exporting window data: {str(e)}")
            return False

    def get_summary_statistics(self) -> Dict:
        """
        获取摘要统计 | Get summary statistics

        Returns:
            dict: 摘要统计信息 | Summary statistics
        """
        if not self.windows:
            self.create_sliding_windows()

        all_means = [w['mean_delta'] for w in self.windows]
        all_snp_counts = [w['snp_count'] for w in self.windows]

        return {
            'total_windows': len(self.windows),
            'mean_delta_mean': np.mean(all_means),
            'mean_delta_std': np.std(all_means),
            'mean_delta_min': np.min(all_means),
            'mean_delta_max': np.max(all_means),
            'mean_snp_count': np.mean(all_snp_counts),
            'min_snp_count': np.min(all_snp_counts),
            'max_snp_count': np.max(all_snp_counts),
            'chromosome_count': len(self.chromosome_lengths)
        }