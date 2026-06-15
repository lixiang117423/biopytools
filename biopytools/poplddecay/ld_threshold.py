"""
LD阈值推荐模块|LD Threshold Recommendation Module

基于Hill & Weir (1988)模型拟合LD衰减曲线，自动推荐合适的r²阈值
Based on Hill & Weir (1988) model to fit LD decay curve and recommend appropriate r² threshold

References:
    Hill, W. G., & Weir, B. S. (1988). Variances and covariances of squared
    linkage disequilibria in finite populations. Theoretical Population
    Biology, 33(1), 54-78.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
import warnings


class LDThresholdRecommender:
    """LD阈值推荐类|LD Threshold Recommender Class"""

    def __init__(self,
                 n: int,
                 dist_col: str = "dist",
                 value_col: str = "r2",
                 dist_unit: str = "bp",
                 max_dist: Optional[float] = None,
                 bg_quantile: float = 0.9,
                 start_rho: List[float] = None,
                 candidate_thresholds: List[float] = None,
                 verbose: bool = True):
        """初始化LD阈值推荐器|Initialize LD Threshold Recommender

        Args:
            n: 个体数量|Number of individuals
            dist_col: 距离列名|Distance column name
            value_col: r²值列名|R² value column name
            dist_unit: 距离单位|Distance unit ("bp" or "kb")
            max_dist: 最大距离（kb）|Maximum distance in kb
            bg_quantile: 背景LD分位数|Background LD quantile
            start_rho: rho参数起始值列表|Starting values for rho parameter
            candidate_thresholds: 候选阈值列表|Candidate threshold list
            verbose: 是否输出详细信息|Whether to print verbose info
        """
        self.n = n
        self.dist_col = dist_col
        self.value_col = value_col
        self.dist_unit = dist_unit
        self.max_dist = max_dist
        self.bg_quantile = bg_quantile
        self.start_rho = start_rho if start_rho else [0.1, 0.01, 1, 0.001, 10]
        self.candidate_thresholds = candidate_thresholds if candidate_thresholds else [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02]
        self.verbose = verbose

        # 结果存储|Result storage
        self.fit_result = None
        self.rho_estimate = None
        self.background_r2 = None
        self.decay_table = None
        self.recommended = None
        self.gwas_window_kb = None

    def hill_weir_formula(self, dist_kb: np.ndarray, rho: float) -> np.ndarray:
        """Hill & Weir期望公式|Hill & Weir expectation formula

        Args:
            dist_kb: 距离（kb）|Distance in kb
            rho: 重组率参数|Recombination rate parameter

        Returns:
            期望r²值|Expected r² values
        """
        numerator = (10 + rho * dist_kb)
        denominator = (22 + 13 * rho * dist_kb + rho**2 * dist_kb**2)
        correction = (1 + ((3 + rho * dist_kb) / (self.n * denominator)))
        return (numerator / denominator) * correction

    def prepare_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """准备数据|Prepare data

        Args:
            data: 输入数据|Input data frame

        Returns:
            清理后的数据|Cleaned data frame
        """
        df = data[[self.dist_col, self.value_col]].copy()
        df.columns = ['dist', 'r2']

        # 过滤无效值|Filter invalid values
        df = df[
            (df['dist'].notna()) &
            (df['r2'].notna()) &
            (df['r2'] > 0) &
            (df['dist'] >= 0)
        ].copy()

        # 转换单位到kb|Convert unit to kb
        if self.dist_unit == "bp":
            df['dist_kb'] = df['dist'] / 1000
        else:
            df['dist_kb'] = df['dist']

        # 过滤最大距离|Filter maximum distance
        if self.max_dist is not None:
            df = df[df['dist_kb'] <= self.max_dist]

        if len(df) < 10:
            raise ValueError(
                f"过滤后有效观测值少于10个|Fewer than 10 valid observations after filtering. "
                f"请检查数据或放宽max_dist限制|Check your data or relax max_dist constraint."
            )

        if self.verbose:
            print(f"[1/4] 数据准备完成: {len(df)}个观测值用于拟合|Data prepared: {len(df)} observations used for fitting.")

        return df

    def fit_model(self, df: pd.DataFrame) -> Tuple[object, float]:
        """拟合Hill & Weir模型|Fit Hill & Weir model

        Args:
            df: 数据|Data frame

        Returns:
            (拟合结果, 使用的rho起始值)|(fit result, rho starting value used)
        """
        fit_result = None
        rho_used = None

        for rho0 in self.start_rho:
            try:
                # 使用curve_fit进行非线性最小二乘拟合|Use curve_fit for nonlinear least squares fitting
                popt, pcov = curve_fit(
                    f=self.hill_weir_formula,
                    xdata=df['dist_kb'].values,
                    ydata=df['r2'].values,
                    p0=[rho0],
                    maxfev=1000,
                    ftol=1e-6,
                    bounds=(0, np.inf)
                )
                fit_result = {'rho': popt[0], 'pcov': pcov}
                rho_used = rho0
                break
            except (RuntimeError, ValueError):
                continue

        if fit_result is None:
            raise ValueError(
                f"Hill & Weir模型对所有起始值都无法收敛|Hill & Weir model failed to converge for all starting values: "
                f"{self.start_rho}. 请提供不同的start_rho值或检查数据质量|Consider providing different start_rho values or checking data quality."
            )

        rho_est = fit_result['rho']

        if self.verbose:
            print(f"[2/4] Hill & Weir模型收敛: rho = {rho_est:.6f} (起始值|start = {rho_used:.4f})")

        return fit_result, rho_est

    def estimate_background_ld(self, df: pd.DataFrame) -> float:
        """估计背景LD|Estimate background LD

        Args:
            df: 数据|Data frame

        Returns:
            背景r²均值|Background r² mean
        """
        dist_bg_cutoff = np.quantile(df['dist_kb'], self.bg_quantile)

        background_r2 = df[df['dist_kb'] >= dist_bg_cutoff]['r2'].mean()

        if self.verbose:
            print(f"[3/4] 背景LD: mean r² = {background_r2:.4f} "
                  f"(距离|distances >= {dist_bg_cutoff:.1f} kb; "
                  f"top {self.bg_quantile*100:.0f}% 分位数|quantile)")

        return background_r2

    def calculate_decay_table(self, df: pd.DataFrame, background_r2: float) -> pd.DataFrame:
        """计算衰减表|Calculate decay table

        Args:
            df: 数据|Data frame
            background_r2: 背景r²|Background r²

        Returns:
            衰减表|Decay table
        """
        # 生成拟合曲线|Generate fit curve
        dist_seq = np.linspace(0, df['dist_kb'].max(), 10000)
        r2_fit = self.hill_weir_formula(dist_seq, self.rho_estimate)
        fit_curve = pd.DataFrame({'dist_kb': dist_seq, 'r2_fit': r2_fit})

        # 计算每个候选阈值的衰减距离|Calculate decay distance for each threshold
        decay_rows = []
        for thr in sorted(self.candidate_thresholds, reverse=True):
            hit = fit_curve[fit_curve['r2_fit'] <= thr]
            if len(hit) == 0:
                continue
            decay_kb = round(hit['dist_kb'].iloc[0], 1)
            decay_rows.append({
                'threshold': thr,
                'decay_kb': decay_kb,
                'bg_ratio': round(thr / background_r2, 2)
            })

        decay_table = pd.DataFrame(decay_rows)

        # 添加推荐建议|Add recommendation
        decay_table['recommendation'] = decay_table['bg_ratio'].apply(
            lambda x: (
                "not recommended: below background LD" if x < 1
                else "caution: close to background LD" if x < 1.5
                else "recommended" if x < 4
                else "caution: far above background LD"
            )
        )

        return decay_table

    def select_best_threshold(self, decay_table: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
        """选择最佳阈值|Select best threshold

        Args:
            decay_table: 衰减表|Decay table

        Returns:
            (最佳阈值行, GWAS窗口kb)|(best threshold row, GWAS window kb)
        """
        recommended_rows = decay_table[decay_table['recommendation'] == "recommended"]

        if len(recommended_rows) == 0:
            warnings.warn(
                "没有阈值符合推荐标准，返回bg_ratio最接近2的阈值|"
                "No threshold met recommendation criteria, returning threshold with bg_ratio closest to 2."
            )
            recommended_rows = decay_table

        # 选择bg_ratio最接近2的行|Select row with bg_ratio closest to 2
        best_idx = np.abs(recommended_rows['bg_ratio'] - 2).idxmin()
        best_row = recommended_rows.loc[[best_idx]].reset_index(drop=True)

        # 计算GWAS窗口（向上取整到50kb）|Calculate GWAS window (round up to 50kb)
        gwas_window_kb = int(np.ceil(best_row['decay_kb'].iloc[0] / 50) * 50)

        return best_row, gwas_window_kb

    def recommend(self, data: Union[pd.DataFrame, str, Path]) -> Dict:
        """推荐LD阈值|Recommend LD threshold

        Args:
            data: 输入数据（DataFrame或文件路径）|Input data (DataFrame or file path)

        Returns:
            结果字典|Result dictionary containing:
                - fit: 拟合结果|Fit result
                - rho: rho估计值|Estimated rho
                - background_r2: 背景r²|Background r²
                - decay_table: 衰减表|Decay table
                - recommended: 推荐阈值|Recommended threshold
                - gwas_window_kb: GWAS窗口kb|GWAS window in kb
        """
        # 读取数据|Read data
        if isinstance(data, (str, Path)):
            # 从文件读取|Read from file
            data_path = Path(data)
            if data_path.suffix == '.gz':
                df = pd.read_csv(data_path, sep='\t', compression='gzip')
            else:
                df = pd.read_csv(data_path, sep='\t')
        else:
            df = data.copy()

        # 准备数据|Prepare data
        df = self.prepare_data(df)

        # 拟合模型|Fit model
        self.fit_result, self.rho_estimate = self.fit_model(df)

        # 估计背景LD|Estimate background LD
        self.background_r2 = self.estimate_background_ld(df)

        # 计算衰减表|Calculate decay table
        self.decay_table = self.calculate_decay_table(df, self.background_r2)

        # 选择最佳阈值|Select best threshold
        self.recommended, self.gwas_window_kb = self.select_best_threshold(self.decay_table)

        if self.verbose:
            print("[4/4] 阈值推荐摘要|Threshold recommendation summary:\n")
            print(self.decay_table.to_string(index=False))
            print(f"\n  最佳阈值|Best threshold: r² = {self.recommended['threshold'].iloc[0]:.2f}  |  "
                  f"衰减距离|decay distance = {self.recommended['decay_kb'].iloc[0]:.1f} kb  |  "
                  f"bg_ratio = {self.recommended['bg_ratio'].iloc[0]:.2f}")
            print(f"  建议的GWAS窗口|Suggested GWAS window: lead SNP +/- {self.gwas_window_kb} kb")

        return {
            'fit': self.fit_result,
            'rho': self.rho_estimate,
            'background_r2': self.background_r2,
            'decay_table': self.decay_table,
            'recommended': self.recommended,
            'gwas_window_kb': self.gwas_window_kb
        }


def recommend_ld_threshold(data: Union[pd.DataFrame, str, Path],
                           n: int,
                           dist_col: str = "dist",
                           value_col: str = "r2",
                           dist_unit: str = "bp",
                           max_dist: Optional[float] = None,
                           bg_quantile: float = 0.9,
                           start_rho: Optional[List[float]] = None,
                           candidate_thresholds: Optional[List[float]] = None,
                           verbose: bool = True) -> Dict:
    """推荐LD阈值（便捷函数）|Recommend LD threshold (convenience function)

    基于Hill & Weir (1988)模型拟合LD衰减曲线，自动推荐合适的r²阈值
    Fit LD decay curve using Hill & Weir (1988) model and recommend appropriate r² threshold

    Args:
        data: 输入数据（DataFrame或文件路径）|Input data (DataFrame or file path)
        n: 个体数量|Number of individuals
        dist_col: 距离列名|Distance column name (default: "dist")
        value_col: r²值列名|R² value column name (default: "r2")
        dist_unit: 距离单位|Distance unit ("bp" or "kb", default: "bp")
        max_dist: 最大距离（kb）|Maximum distance in kb (default: None)
        bg_quantile: 背景LD分位数|Background LD quantile (default: 0.9)
        start_rho: rho参数起始值列表|Starting values for rho parameter
        candidate_thresholds: 候选阈值列表|Candidate threshold list
        verbose: 是否输出详细信息|Whether to print verbose info (default: True)

    Returns:
        结果字典|Result dictionary

    Examples:
        >>> # 从DataFrame推荐|Recommend from DataFrame
        >>> result = recommend_ld_threshold(ld_data, n=805)
        >>> print(result['recommended'])
        >>> print(result['gwas_window_kb'])

        >>> # 从文件推荐|Recommend from file
        >>> result = recommend_ld_threshold("ld.stat.gz", n=805, dist_unit="bp")
    """
    recommender = LDThresholdRecommender(
        n=n,
        dist_col=dist_col,
        value_col=value_col,
        dist_unit=dist_unit,
        max_dist=max_dist,
        bg_quantile=bg_quantile,
        start_rho=start_rho,
        candidate_thresholds=candidate_thresholds,
        verbose=verbose
    )

    return recommender.recommend(data)
