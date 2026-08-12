"""mixrace 判读引擎(纯逻辑)|mixrace verdict engine (pure logic).

四指标透明投票 → single_genotype / mixed_genotype / uncertain + 置信 + 人话依据。
指标基于读计数,与倍性无关;阈值全部可配(默认参考值,需用已知纯样品校准)。
|Transparent 4-metric voting -> verdict + confidence + rationale. Read-count based,
ploidy-independent; thresholds configurable (defaults need calibration with known-pure samples).
"""
import statistics
from typing import Optional


def _band(value: float, low: float, high: float) -> str:
    """三段判定|three-band classification: <low=pure, >high=mixed, else gray."""
    if value < low:
        return "pure"
    if value > high:
        return "mixed"
    return "gray"


def judge(metrics: dict, thr: dict) -> dict:
    """判读|judge mixed vs single genotype.

    Args:
        metrics: {vaf_mid_ratio, multiallelic_ratio, fws(float|None), hw, mean_depth, heterozygosity?}
        thr: {vaf_mid_low, vaf_mid_high, multiallelic_low, multiallelic_high, fws_cutoff, min_depth}

    Returns:
        {verdict, confidence, votes, rationale};多等位 mixed 权重最高(任一 mixed 即判混合)。
        |multiallelic carries highest weight (any mixed vote -> mixed_genotype).
    """
    vaf = metrics["vaf_mid_ratio"]
    ma = metrics["multiallelic_ratio"]
    fws = metrics.get("fws")
    depth = metrics["mean_depth"]

    votes = {
        "vaf_mid_ratio": _band(vaf, thr["vaf_mid_low"], thr["vaf_mid_high"]),
        "multiallelic_ratio": _band(ma, thr["multiallelic_low"], thr["multiallelic_high"]),
    }
    # Fws:>=cutoff 纯,<cutoff 杂;None(单样本)不投票|Fws vote (None=n/a, single sample)
    if fws is None:
        votes["fws"] = "n/a"
    else:
        votes["fws"] = "pure" if fws >= thr["fws_cutoff"] else "mixed"

    low_depth = depth < thr["min_depth"]
    active = [v for v in votes.values() if v != "n/a"]
    if any(v == "mixed" for v in active):
        verdict = "mixed_genotype"
    elif all(v == "pure" for v in active):
        verdict = "single_genotype"
    else:
        verdict = "uncertain"
    confidence = "low" if low_depth else "high"
    return {
        "verdict": verdict,
        "confidence": confidence,
        "votes": votes,
        "rationale": _rationale(verdict, metrics, thr, low_depth),
    }


def _rationale(verdict: str, m: dict, thr: dict, low_depth: bool) -> str:
    """生成人话依据串|build human-readable rationale."""
    label = {"single_genotype": "疑似纯", "mixed_genotype": "疑似混合",
             "uncertain": "不确定"}[verdict]
    parts = [f"{label}（{verdict}）："]
    parts.append(
        f"VAF中间频率占比 {m['vaf_mid_ratio']*100:.1f}%"
        f"（阈值<{thr['vaf_mid_low']*100:.0f}%纯 / >{thr['vaf_mid_high']*100:.0f}%杂）")
    parts.append(
        f"多等位位点占比 {m['multiallelic_ratio']*100:.2f}%"
        f"（<{thr['multiallelic_low']*100:.1f}%纯 / >{thr['multiallelic_high']*100:.0f}%杂，硬指标）")
    if m.get("fws") is not None:
        parts.append(f"Fws={m['fws']:.3f}（cutoff {thr['fws_cutoff']}）")
    else:
        parts.append(f"Fws=N/A（单样本无群体基线，仅参考Hw={m.get('hw', 0):.3f}）")
    if low_depth:
        parts.append(f"平均深度 {m['mean_depth']:.1f}x（深度不足，VAF指标可能失真）")
    else:
        parts.append(f"平均深度 {m['mean_depth']:.1f}x（充足）")
    if m.get("heterozygosity") is not None:
        parts.append(f"GenomeScope杂合度 {m['heterozygosity']*100:.2f}%")
    return "；".join(parts) + "。"


def calibrate_thresholds(pure_metrics: list) -> dict:
    """用已知纯样品按 mean+2SD 算建议阈值|calibrate suggested thresholds (mean+2SD).

    Args:
        pure_metrics: [{vaf_mid_ratio, multiallelic_ratio}, ...] 已知纯样品指标
    Returns:
        {vaf_mid_low, multiallelic_low}(样本<2 或缺键返回 None)|thresholds (None if <2 samples)
    """
    def _mean_plus_2sd(key: str) -> Optional[float]:
        xs = [p[key] for p in pure_metrics if key in p]
        if len(xs) < 2:
            return None
        return statistics.mean(xs) + 2 * statistics.pstdev(xs)

    return {
        "vaf_mid_low": _mean_plus_2sd("vaf_mid_ratio"),
        "multiallelic_low": _mean_plus_2sd("multiallelic_ratio"),
    }
