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
        """计算阈值（精确复现DeepBSA逻辑）|Calculate threshold (exact reproduction of DeepBSA logic)

        源码位置|Source location: Statistic_Methods.py:284-293

        公式|Formula: threshold = median + 3 * std
        备选|Alternative: 90th percentile (if all data < threshold)
        舍入|Rounding: ROUND_HALF_UP to 4 decimal places

        Args:
            data: 平滑后的数据数组|Smoothed data array

        Returns:
            float: 阈值|Threshold value
        """
        from decimal import Decimal

        std = np.std(data)
        med = np.median(data)
        threshold = med + 3 * std

        if np.all(data < round(threshold, 4)):
            threshold = np.percentile(data, 90)

        threshold = float(Decimal(str(threshold)).quantize(Decimal("0.0001"), rounding="ROUND_HALF_UP"))

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


def _read_chrom_positions_from_values(values_path: Path, logger=None) -> tuple:
    """从 values.txt 读取染色体名和位置（按出现顺序分组）|Read chromosome names and positions from values.txt

    values.txt 格式|Format: chr\\tpos\\tvalue，每行一个数据点，按染色体连续排列

    Args:
        values_path: values.txt 文件路径|Path to values.txt
        logger: 日志器|Logger

    Returns:
        tuple: (chrom_names, chrom_positions) - 染色体名列表和对应的位置列表
    """
    chrom_names = []
    chrom_positions = []
    current_chrom = None
    current_positions = []

    with open(values_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 3:
                continue
            chrom, pos, value = parts
            if chrom != current_chrom:
                if current_chrom is not None:
                    chrom_names.append(current_chrom)
                    chrom_positions.append(current_positions)
                current_chrom = chrom
                current_positions = []
            current_positions.append(int(pos))

    if current_chrom is not None:
        chrom_names.append(current_chrom)
        chrom_positions.append(current_positions)

    return chrom_names, chrom_positions


def extract_plot_data_from_npy(input_dir: Path, methods: List[str],
                                logger=None) -> pd.DataFrame:
    """从 npy 文件读取绘图数据（与原始软件完全一致）|Read plot data from npy files (exact reproduction)

    从 DeepBSA 原始输出的 npy 文件读取原始值和平滑值，从 values.txt 读取染色体名和位置，
    计算阈值后组装 DataFrame。结果与原始软件绘图数据完全一致。

    文件结构|File structure:
        input_dir/all_data_for_plot_{method}.npy       # 原始值（按染色体分组）
        input_dir/smooth_data_for_plot_{method}.npy     # 平滑值（按染色体分组）
        input_dir/{method} values.txt                   # 位置和染色体名

    Args:
        input_dir: 包含 npy 和 values.txt 的目录|Directory containing npy and values.txt
        methods: 方法列表|Method list (e.g. ["DL", "ED4", "SNP"])
        logger: 日志器|Logger

    Returns:
        pd.DataFrame: 包含绘图数据的数据表|Data table with plot data
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("=" * 60)
    logger.info("从npy文件读取绘图数据|Reading plot data from npy files")
    logger.info("=" * 60)

    all_rows = []

    for method in methods:
        # SNP 方法在 npy 文件名中用 ΔSNP
        npy_method_name = "\u0394SNP" if method == "SNP" else method
        display_name = "\u0394SNP" if method == "SNP" else method

        # 查找 npy 文件|Find npy files
        raw_path = input_dir / f"all_data_for_plot_{npy_method_name}.npy"
        smooth_path = input_dir / f"smooth_data_for_plot_{npy_method_name}.npy"

        if not raw_path.exists() or not smooth_path.exists():
            logger.warning(f"npy文件不存在|npy files not found for {method}: {raw_path.name}, {smooth_path.name}")
            continue

        # 查找 values.txt|Find values.txt
        values_pattern = f"{npy_method_name} values.txt"
        values_files = list(input_dir.glob(values_pattern))
        if not values_files:
            logger.warning(f"values.txt不存在|values.txt not found for {method}: {values_pattern}")
            continue

        values_path = values_files[0]

        logger.info(f"处理|Processing: {method}")
        logger.info(f"  原始值|Raw: {raw_path.name}")
        logger.info(f"  平滑值|Smooth: {smooth_path.name}")
        logger.info(f"  位置|Position: {values_path.name}")

        try:
            # 加载 npy|Load npy files
            raw = np.load(raw_path, allow_pickle=True)
            smooth = np.load(smooth_path, allow_pickle=True)

            n_chroms_raw = len(raw)
            n_chroms_smooth = len(smooth)

            # 从 values.txt 读取染色体名和位置|Read chrom names and positions from values.txt
            chrom_names, chrom_positions = _read_chrom_positions_from_values(values_path, logger)
            n_chroms_txt = len(chrom_names)

            if n_chroms_raw != n_chroms_smooth or n_chroms_raw != n_chroms_txt:
                logger.error(f"  染色体数量不匹配|Chromosome count mismatch: "
                             f"raw={n_chroms_raw}, smooth={n_chroms_smooth}, txt={n_chroms_txt}")
                continue

            # 验证每条染色体长度一致|Verify length consistency per chromosome
            for i in range(n_chroms_raw):
                len_raw = len(raw[i])
                len_smooth = len(smooth[i])
                len_pos = len(chrom_positions[i])
                if len_raw != len_smooth or len_raw != len_pos:
                    logger.error(f"  Chr{chrom_names[i]} 长度不匹配|length mismatch: "
                                 f"raw={len_raw}, smooth={len_smooth}, pos={len_pos}")
                    continue

            # 计算阈值|Calculate threshold
            smooth_concat = np.concatenate([np.array(smooth[i], dtype=np.float64) for i in range(n_chroms_raw)])
            calculator = PlotDataCalculator(logger)
            threshold = calculator.calculate_threshold(smooth_concat)

            logger.info(f"  阈值|Threshold: {threshold:.4f}")

            # 组装数据|Assemble data
            above_count = 0
            for chrom_idx in range(n_chroms_raw):
                chrom_name = chrom_names[chrom_idx]
                positions = chrom_positions[chrom_idx]
                raw_vals = raw[chrom_idx]
                smooth_vals = smooth[chrom_idx]

                for pos, raw_val, smooth_val in zip(positions, raw_vals, smooth_vals):
                    is_above = float(smooth_val) >= threshold
                    if is_above:
                        above_count += 1
                    all_rows.append({
                        'Method': display_name,
                        'Chromosome': chrom_name,
                        'Position': int(pos),
                        'Position_Mb': pos / 1e6,
                        'Raw_Value': float(raw_val),
                        'Smooth_Value': float(smooth_val),
                        'Threshold': threshold,
                        'Above_Threshold': is_above
                    })

            logger.info(f"  数据点|Data points: {n_chroms_raw} 条染色体, "
                        f"{sum(len(chrom_positions[i]) for i in range(n_chroms_raw))} 个数据点")
            logger.info(f"  超阈值点|Points above threshold: {above_count}")

        except Exception as e:
            logger.error(f"  处理失败|Processing failed for {method}: {e}")

    if not all_rows:
        logger.error("没有从npy文件中读取到任何数据|No data read from npy files")
        return None

    result_df = pd.DataFrame(all_rows)

    logger.info("")
    logger.info("=" * 60)
    logger.info("npy数据读取完成|npy data reading completed")
    logger.info("=" * 60)
    logger.info(f"总计|Total: {len(result_df)} 个数据点|data points")
    logger.info(f"方法|Methods: {', '.join(sorted(result_df['Method'].unique()))}")
    logger.info("")
    logger.info("各方法阈值|Thresholds by method:")
    for method in sorted(result_df['Method'].unique()):
        threshold = result_df[result_df['Method'] == method]['Threshold'].iloc[0]
        count_above = result_df[result_df['Method'] == method]['Above_Threshold'].sum()
        logger.info(f"  {method}: 阈值|threshold={threshold:.4f}, 超阈值点|above threshold={count_above}")

    return result_df


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

        # 查找values.txt文件|Find values.txt files (递归查找，因为子目录名基于输入文件名|Recursive search, subdir name based on input filename)
        txt_files = list(method_dir.glob("Results/**/* values.txt"))

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

    # 修改SNP方法名为ΔSNP|Change SNP method name to ΔSNP
    df['Method'] = df['Method'].apply(lambda x: 'ΔSNP' if x == 'SNP' else x)

    # 2. 计算平滑曲线和阈值|Calculate smooth curves and thresholds
    result_df = calculator.process_all_methods(df, smooth_func=smooth_func, smooth_frac=smooth_frac)

    return result_df
