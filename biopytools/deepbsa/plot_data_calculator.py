"""
DeepBSA绘图数据计算模块|DeepBSA Plot Data Calculator Module

复现DeepBSA的平滑和阈值计算逻辑，在Python中完成所有计算
Reproduce DeepBSA's smoothing and threshold calculation logic in Python
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple
from pathlib import Path
import logging


class PlotDataCalculator:
    """绘图数据计算器|Plot Data Calculator"""

    def __init__(self, logger=None):
        """初始化计算器|Initialize calculator

        Args:
            logger: 日志器|Logger
        """
        self.logger = logger or logging.getLogger(__name__)

    def calculate_threshold(self, data: np.ndarray) -> float:
        """计算阈值（复现DeepBSA逻辑）|Calculate threshold (reproduce DeepBSA logic)

        源码位置|Source location: Statistic_Methods.py:282-291

        公式|Formula: threshold = median + 3 * std
        备选|Alternative: 90th percentile (if all data < threshold)

        Args:
            data: 平滑后的数据数组|Smoothed data array

        Returns:
            float: 阈值|Threshold value
        """
        from decimal import Decimal

        # 计算标准差和中位数|Calculate std and median
        std = np.std(data)
        med = np.median(data)
        threshold = med + 3 * std

        # 如果所有数据都低于阈值，使用90%分位数
        # If all data below threshold, use 90th percentile
        if np.all(data < round(threshold, 4)):
            threshold = np.percentile(data, 90)

        # 四舍五入到4位小数|Round to 4 decimal places
        threshold = float(Decimal(str(threshold)).quantize(Decimal("0.0001")))

        self.logger.debug(f"阈值计算|Threshold calculation: median={med:.4f}, std={std:.4f}, threshold={threshold:.4f}")

        return threshold

    def smooth_lowess(self, data: np.ndarray, frac: float = 0.1) -> np.ndarray:
        """LOWESS平滑（复现DeepBSA逻辑）|LOWESS smoothing (reproduce DeepBSA logic)

        源码位置|Source location: Statistic_Methods.py:135-138

        Args:
            data: 原始数据数组|Raw data array
            frac: 窗口大小比例|Window size fraction (default: 0.1)

        Returns:
            np.ndarray: 平滑后的数据|Smoothed data
        """
        try:
            from statsmodels.nonparametric.lowess import lowess
            # LOWESS平滑|LOWESS smoothing
            z = lowess(data, np.arange(len(data)), frac=frac, it=0, delta=0, return_sorted=False)
            return np.array(z)
        except ImportError:
            self.logger.warning("statsmodels未安装，使用scipy实现的LOWESS|statsmodels not installed, using scipy implementation")
            return self._smooth_lowess_scipy(data, frac)

    def _smooth_lowess_scipy(self, data: np.ndarray, frac: float) -> np.ndarray:
        """使用scipy实现LOWESS平滑|LOWESS smoothing using scipy

        Args:
            data: 原始数据数组|Raw data array
            frac: 窗口大小比例|Window size fraction

        Returns:
            np.ndarray: 平滑后的数据|Smoothed data
        """
        from scipy.interpolate import interp1d
        from scipy.signal import savgol_filter

        # 计算窗口大小|Calculate window size
        window_size = max(3, int(len(data) * frac))
        if window_size % 2 == 0:
            window_size += 1  # 确保是奇数|Ensure odd number

        # 使用Savitzky-Golay滤波器（近似LOWESS）|Use Savitzky-Golay filter (approximate LOWESS)
        try:
            smoothed = savgol_filter(data, window_size, polyorder=2, mode='interp')
            return smoothed
        except:
            # 如果窗口太大，使用简单的移动平均|If window too large, use simple moving average
            from numpy.convolve import convolve
            kernel = np.ones(window_size) / window_size
            smoothed = np.convolve(data, kernel, mode='same')
            return smoothed

    def smooth_tri_cube_kernel(self, data: np.ndarray, window_size: float = 0.1, step: int = 1) -> np.ndarray:
        """Tri-cube核平滑（复现DeepBSA逻辑）|Tri-cube kernel smoothing (reproduce DeepBSA logic)

        源码位置|Source location: Statistic_Methods.py:91-125

        Args:
            data: 原始数据数组|Raw data array
            window_size: 窗口大小比例|Window size fraction (default: 0.1)
            step: 步长|Step size (default: 1)

        Returns:
            np.ndarray: 平滑后的数据|Smoothed data
        """
        # Tri-cube核函数|Tri-cube kernel function
        weight_cal = lambda x: (1 - x ** 3) ** 3

        # 计算窗口大小|Calculate window size
        n_window = int(len(data) * window_size)

        # 确保窗口大小是奇数|Ensure window size is odd
        if n_window % 2 == 0:
            n_window += 1

        center = int((n_window + 1) / 2)
        left_offset = center - 1
        right_offset = n_window - center

        # 数据填充|Data padding
        padding_value = 0
        head_padding = np.array([padding_value] * left_offset)
        tail_padding = np.array([padding_value] * (right_offset + step - (len(data) - 1) % step))
        padded_data = np.concatenate([head_padding, data, tail_padding])

        # 构造距离向量|Construct distance vector
        distance_vector = np.arange(len(padded_data))
        result = []

        # 滑动窗口处理|Sliding window processing
        for c_index in range(center - 1, len(padded_data) - right_offset - (len(data) - 1) % step, step):
            # 提取窗口数据|Extract window data
            window_data = padded_data[c_index - left_offset: c_index + right_offset + 1]
            window_distance = np.abs(
                distance_vector[c_index - left_offset: c_index + right_offset + 1] -
                distance_vector[c_index]
            )

            # 计算距离权重|Calculate distance weights
            if np.max(window_distance) > 0:
                distance_weight = weight_cal(window_distance / np.max(window_distance))
            else:
                distance_weight = np.ones_like(window_distance)

            # 归一化权重|Normalize weights
            k_weight = distance_weight / np.sum(distance_weight)

            # 加权平均|Weighted average
            smoothed_value = np.sum(np.dot(k_weight, window_data))
            result.append(smoothed_value)

        result = np.array(result)
        result = np.nan_to_num(result)

        return result

    def calculate_all_for_method(self, df: pd.DataFrame, method: str,
                                  smooth_func: str = "LOWESS",
                                  smooth_frac: float = 0.1) -> pd.DataFrame:
        """为一个方法计算所有绘图数据|Calculate all plot data for one method

        Args:
            df: 原始数据DataFrame|Raw data DataFrame
            method: 方法名|Method name
            smooth_func: 平滑函数|Smoothing function ("LOWESS" or "Tri-kernel")
            smooth_frac: 平滑窗口比例|Smoothing window fraction

        Returns:
            pd.DataFrame: 包含原始值、平滑值、阈值的DataFrame|DataFrame with raw, smooth, threshold
        """
        self.logger.info(f"计算|Calculating {method}: 平滑|smoothing={smooth_func}, 窗口|window={smooth_frac}")

        result_data = []

        for chrom in df['Chromosome'].unique():
            df_chrom = df[df['Chromosome'] == chrom].sort_values('Position')
            positions = df_chrom['Position'].values
            raw_values = df_chrom['Value'].values

            # 平滑处理|Smoothing
            if smooth_func == "LOWESS":
                smooth_values = self.smooth_lowess(raw_values, frac=smooth_frac)
            elif smooth_func == "Tri-kernel":
                smooth_values = self.smooth_tri_cube_kernel(raw_values, window_size=smooth_frac)
            else:
                self.logger.warning(f"未知的平滑函数|Unknown smoothing function: {smooth_func}，使用LOWESS")
                smooth_values = self.smooth_lowess(raw_values, frac=smooth_frac)

            # 确保长度一致|Ensure length consistency
            if len(smooth_values) != len(raw_values):
                self.logger.warning(f"平滑后长度不匹配|Smoothed length mismatch: {len(smooth_values)} vs {len(raw_values)}")
                # 调整长度|Adjust length
                if len(smooth_values) < len(raw_values):
                    smooth_values = np.pad(smooth_values, (0, len(raw_values) - len(smooth_values)), 'edge')
                else:
                    smooth_values = smooth_values[:len(raw_values)]

            # 添加到结果|Add to results
            for pos, raw_val, smooth_val in zip(positions, raw_values, smooth_values):
                result_data.append({
                    'Method': method,
                    'Chromosome': chrom,
                    'Position': pos,
                    'Position_Mb': pos / 1e6,
                    'Raw_Value': raw_val,
                    'Smooth_Value': smooth_val
                })

        result_df = pd.DataFrame(result_data)

        # 计算阈值|Calculate threshold
        all_smooth_values = result_df['Smooth_Value'].values
        threshold = self.calculate_threshold(all_smooth_values)

        # 添加阈值列|Add threshold column
        result_df['Threshold'] = threshold

        # 添加阈值指示|Add threshold indicator
        result_df['Above_Threshold'] = result_df['Smooth_Value'] >= threshold

        self.logger.info(f"  阈值|Threshold: {threshold:.4f}")
        self.logger.info(f"  超阈值的点|Points above threshold: {result_df['Above_Threshold'].sum()}")

        return result_df

    def process_all_methods(self, df: pd.DataFrame,
                            smooth_func: str = "LOWESS",
                            smooth_frac: float = 0.1) -> pd.DataFrame:
        """处理所有方法|Process all methods

        Args:
            df: 原始数据DataFrame|Raw data DataFrame
            smooth_func: 平滑函数|Smoothing function
            smooth_frac: 平滑窗口比例|Smoothing window fraction

        Returns:
            pd.DataFrame: 包含所有方法计算结果的DataFrame|DataFrame with all methods' results
        """
        self.logger.info("=" * 60)
        self.logger.info("计算平滑曲线和阈值|Calculating smooth curves and thresholds")
        self.logger.info("=" * 60)
        self.logger.info(f"平滑方法|Smoothing method: {smooth_func}")
        self.logger.info(f"窗口比例|Window fraction: {smooth_frac}")

        all_results = []

        for method in df['Method'].unique():
            df_method = df[df['Method'] == method]
            result_df = self.calculate_all_for_method(df_method, method, smooth_func, smooth_frac)
            all_results.append(result_df)

        # 合并所有方法|Merge all methods
        final_df = pd.concat(all_results, ignore_index=True)

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("计算完成|Calculation completed")
        self.logger.info("=" * 60)
        self.logger.info(f"总计|Total: {len(final_df)} 个数据点|data points")
        self.logger.info(f"方法数量|Methods: {len(final_df['Method'].unique())}")

        return final_df


def extract_and_calculate_plot_data(each_dir: Path, methods: List[str],
                                    logger=None,
                                    smooth_func: str = "LOWESS",
                                    smooth_frac: float = 0.1) -> pd.DataFrame:
    """提取并计算绘图数据|Extract and calculate plot data

    Args:
        each_dir: each目录路径|Path to each directory
        methods: 方法列表|Method list
        logger: 日志器|Logger
        smooth_func: 平滑函数|Smoothing function
        smooth_frac: 平滑窗口比例|Smoothing window fraction

    Returns:
        pd.DataFrame: 包含原始值、平滑值、阈值的数据表|Data table with raw, smooth, threshold
    """
    calculator = PlotDataCalculator(logger)

    # 1. 提取原始数据|Extract raw data
    all_raw_data = []

    for method in methods:
        method_dir = each_dir / method
        if not method_dir.exists():
            if logger:
                logger.warning(f"方法目录不存在|Method directory not found: {method_dir}")
            continue

        # 查找values.txt文件|Find values.txt files
        txt_files = list(method_dir.glob("Results/variant/* values.txt"))

        if not txt_files:
            if logger:
                logger.warning(f"未找到绘图数据文件|No plot data file found for {method}")
            continue

        for txt_file in txt_files:
            if logger:
                logger.info(f"读取|Reading: {txt_file}")

            try:
                # 读取values.txt文件|Read values.txt file
                data_list = []
                with open(txt_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split('\t')
                        if len(parts) == 3:
                            chromosome, position, value = parts
                            data_list.append({
                                'Method': method,
                                'Chromosome': chromosome,
                                'Position': int(position),
                                'Value': float(value)
                            })

                if data_list:
                    all_raw_data.extend(data_list)
                    if logger:
                        logger.info(f"  提取|Extracted {len(data_list)} 个数据点|data points")

            except Exception as e:
                if logger:
                    logger.error(f"  读取文件失败|Failed to read file: {e}")

    if not all_raw_data:
        if logger:
            logger.error("没有找到任何绘图数据|No plot data found")
        return None

    # 创建DataFrame|Create DataFrame
    df = pd.DataFrame(all_raw_data)

    # 修改SNP方法名为δSNP|Change SNP method name to δSNP
    df['Method'] = df['Method'].apply(lambda x: 'δSNP' if x == 'SNP' else x)

    # 2. 计算平滑曲线和阈值|Calculate smooth curves and thresholds
    result_df = calculator.process_all_methods(df, smooth_func=smooth_func, smooth_frac=smooth_frac)

    return result_df
