"""LDThresholdRecommender单元测试|Unit tests for LDThresholdRecommender

覆盖背景r²偏高时不可达阈值(如0.5)被误推荐的bug修复|Covers the fix where
structurally unattainable thresholds (e.g. r2=0.5 above Hill-Weir curve max)
were recommended with a degenerate decay distance of 0.
"""

import numpy as np
import pandas as pd

from biopytools.poplddecay.ld_threshold import LDThresholdRecommender


def _make_decay_data(rho: float = 0.006, n: int = 60, max_dist_kb: float = 300,
                     n_points: int = 120, noise: float = 0.02, seed: int = 42) -> pd.DataFrame:
    """合成H&W形状的LD衰减数据|Synthesize LD decay data following the Hill-Weir shape

    返回PopLDdecay stat格式列(#Dist为bp, Mean_r^2)|Returns PopLDdecay stat-format columns.
    """
    rng = np.random.default_rng(seed)
    dist_kb = np.linspace(0, max_dist_kb, n_points)
    x = rho * dist_kb
    r2 = (10 + x) / (22 + 13 * x + x ** 2) * (1 + (3 + x) / (n * (22 + 13 * x + x ** 2)))
    r2 = np.clip(r2 + rng.normal(0, noise, n_points), 0.01, 0.9)
    return pd.DataFrame({'#Dist': dist_kb * 1000, 'Mean_r^2': r2})


def _make_recommender(n: int = 60) -> LDThresholdRecommender:
    return LDThresholdRecommender(
        n=n, dist_col='#Dist', value_col='Mean_r^2', dist_unit='bp', verbose=False)


class TestHillWeirCurveMax:
    """H&W公式结构性质|Structural properties of the Hill-Weir formula"""

    def test_curve_max_always_below_half(self):
        """d=0处最大值恒<0.5,故r²=0.5候选结构性不可达|Max at d=0 is always <0.5, so r2=0.5 is unattainable"""
        rec = _make_recommender()
        for n in [2, 10, 60, 1000]:
            rec.n = n
            r2_at_zero = rec.hill_weir_formula(np.array([0.0]), 0.001)[0]
            assert r2_at_zero < 0.5, f"n={n}时曲线起点应<0.5|curve start should be <0.5 at n={n}"


class TestUnattainableThreshold:
    """高背景r²场景的bug修复|Bug fix for high-background-r2 scenario"""

    def test_unattainable_threshold_not_recommended(self):
        """背景r²>0.22时不可选衰减距离为0的0.5|Must not recommend 0.5 with decay 0 when bg r2>0.22"""
        # rho=0.006 -> 背景r²约0.24,复现真实案例|rho=0.006 -> bg r2 ~0.24, reproduces the real case
        data = _make_decay_data(rho=0.006)
        result = _make_recommender().recommend(data)

        best = result['recommended'].iloc[0]
        assert best['threshold'] != 0.5, "不可达阈值不应被推荐|Unattainable threshold must not be recommended"
        assert np.isfinite(best['decay_kb']) and best['decay_kb'] > 0, \
            f"推荐行衰减距离应>0|Recommended decay must be >0, got {best['decay_kb']}"
        assert result['gwas_window_kb'] > 0, \
            f"GWAS窗口应>0|GWAS window must be >0, got {result['gwas_window_kb']}"

    def test_unattainable_threshold_marked_in_decay_table(self):
        """不可达行应在衰减表中显式标记|Unattainable rows must be explicitly marked in the decay table"""
        data = _make_decay_data(rho=0.006)
        result = _make_recommender().recommend(data)

        row = result['decay_table'][result['decay_table']['threshold'] == 0.5]
        assert len(row) == 1, "默认候选含0.5|Default candidates include 0.5"
        assert pd.isna(row['decay_kb'].iloc[0]), \
            "不可达行衰减距离应为NaN|Unattainable row decay_kb should be NaN"
        assert 'above fitted curve max' in row['recommendation'].iloc[0], \
            f"不可达行应标记原因|Unattainable row should say why, got {row['recommendation'].iloc[0]}"

    def test_attainable_row_keeps_decay_distance(self):
        """可达行(0.4)仍给出真实衰减距离|Attainable row (0.4) still reports a real decay distance"""
        data = _make_decay_data(rho=0.006)
        result = _make_recommender().recommend(data)

        row = result['decay_table'][result['decay_table']['threshold'] == 0.4]
        assert len(row) == 1
        assert 10 < row['decay_kb'].iloc[0] < 100, \
            f"rho~0.006时0.4阈值衰减距离应约30-50kb|Expected ~30-50kb at thr 0.4, got {row['decay_kb'].iloc[0]}"


class TestSelectBestThreshold:
    """select_best_threshold回退路径|Fallback path of select_best_threshold"""

    def test_fallback_prefers_attainable_rows(self):
        """无recommended行时回退应避开不可达行|Fallback must avoid unattainable rows"""
        rec = _make_recommender()
        decay_table = pd.DataFrame({
            'threshold': [0.5, 0.4],
            'decay_kb': [np.nan, 30.0],
            'bg_ratio': [2.1, 1.7],
            'recommendation': [
                'not recommended: above fitted curve max',
                'caution: close to background LD',
            ],
        })
        best, window = rec.select_best_threshold(decay_table)
        assert best['threshold'].iloc[0] == 0.4, "回退应选可达行|Fallback should pick the attainable row"
        assert window == 50, f"30kb向上取整到50kb|30kb rounds up to 50kb, got {window}"


class TestNormalScenario:
    """低背景r²正常场景回归守卫|Regression guard for the normal low-background scenario"""

    def test_low_background_recommendation_unchanged(self):
        """背景r²~0.06时仍推荐0.1并给出真实衰减距离|bg r2~0.06 still recommends 0.1 with real decay"""
        # rho=0.05快速衰减 -> 背景r²约0.06|rho=0.05 fast decay -> bg r2 ~0.06
        data = _make_decay_data(rho=0.05)
        result = _make_recommender().recommend(data)

        best = result['recommended'].iloc[0]
        assert best['threshold'] == 0.1, \
            f"bg_ratio~1.7应选0.1|bg_ratio~1.7 should pick 0.1, got {best['threshold']}"
        assert best['decay_kb'] > 50, f"衰减距离应为正|Decay should be positive, got {best['decay_kb']}"
        assert result['gwas_window_kb'] > 0
