"""mixrace 判读引擎(单倍体)|mixrace verdict engine (haploid).

指标:杂合率 het_rate(主指标,<0.01%/0.1%/1% 阈值)+ AFS 形态 + 优势株占比 + 平均深度。
输出 single_genotype/mixed_genotype/uncertain + 置信 + 人话依据。Fws 已弃用(疟原虫二倍体指标)。
|Metrics: het_rate (primary) + AFS shape + dominant proportion + depth.
Fws dropped (malaria diploid metric, invalid here).
"""
from typing import Optional


def judge(metrics: dict, thr: dict) -> dict:
    """判读|judge pure vs mixed genotype (haploid).

    Args:
        metrics: {het_rate, afs_shape, dominant_proportion?, mean_depth, heterozygosity?, no_data?}
        thr: {het_pure, het_suspicious, het_impure, min_depth}
    """
    if metrics.get("no_data"):
        # 无变异数据 ≠ 纯:数据缺失与真纯必须可区分,方向性判错不可接受|
        # No variant data != pure: missing data must be distinguishable from truly pure.
        return {"verdict": "uncertain", "confidence": "low",
                "rationale": "未获取到变异数据(过滤后VCF缺失或查询失败),无法判读;"
                             "请先完成 --step 5|No variant data (filtered VCF missing or "
                             "query failed); run --step 5 first."}
    het = metrics["het_rate"]
    shape = metrics.get("afs_shape", "monoclonal")
    dom = metrics.get("dominant_proportion")
    depth = metrics["mean_depth"]

    # 杂合率为主要判据;AFS 形态辅助细化灰色区。
    # |het_rate is the primary criterion; AFS shape refines the gray zone.
    if het < thr["het_pure"]:
        verdict = "single_genotype"        # 杂合率极低 → 基本纯(少量中频位点算噪声)|near-zero het -> pure
    elif het >= thr["het_impure"]:
        verdict = "mixed_genotype"         # >1% 不纯|impure
    elif het >= thr["het_suspicious"]:     # 0.1%–1% 可疑|suspicious
        verdict = "mixed_genotype" if shape != "monoclonal" else "uncertain"
    else:                                  # 0.01%–0.1% 排查区|investigate
        verdict = "uncertain"

    low_depth = depth < thr["min_depth"]
    confidence = "low" if low_depth else "high"
    return {"verdict": verdict, "confidence": confidence,
            "rationale": _rationale(verdict, metrics, thr, low_depth, dom)}


def _rationale(verdict: str, m: dict, thr: dict, low_depth: bool, dom: Optional[float]) -> str:
    label = {"single_genotype": "疑似纯(单克隆)", "mixed_genotype": "疑似混合",
             "uncertain": "不确定(灰色区,需排查 repeat/旁系同源/污染)"}[verdict]
    parts = [f"{label}|{verdict}"]   # 勿带尾冒号,join("；")已分隔|no trailing colon (join adds ；)
    parts.append(
        f"全基因组杂合率 {m['het_rate']*100:.4f}%（判定标准 <{thr['het_pure']*100:.3f}%纯 / "
        f">{thr['het_impure']*100:.0f}%不纯），共 {m.get('het_sites', 0)} 个杂合位点")
    parts.append(f"AFS 形态|AFS shape: {m.get('afs_shape','monoclonal')}")
    if dom is not None and verdict != "single_genotype":
        # 仅混合/可疑时报优势株占比;纯样品杂合位点太少(多为噪声),占比无意义
        # |report dominant only for mixed/uncertain; pure samples have too few (noisy) het sites
        parts.append(f"优势株占比 ≈ {dom*100:.1f}%（混合位点主等位频率中位数,方法B）")
    if low_depth:
        parts.append(f"平均深度 {m['mean_depth']:.1f}x（< {thr['min_depth']:.0f}x,建议≥50x 以保证结果可靠,置信降级）")
    else:
        parts.append(f"平均深度 {m['mean_depth']:.1f}x（充足）")
    if m.get("heterozygosity") is not None:
        parts.append(f"GenomeScope 杂合度 {m['heterozygosity']*100:.2f}%")
    return "；".join(parts) + "。"


def calibrate_thresholds(pure_metrics: list) -> dict:
    """用已知纯样品按 mean+2SD 校准 het 阈值|calibrate het thresholds (mean+2SD of pure samples)."""
    import statistics as st
    def _ms(key):
        xs = [p[key] for p in pure_metrics if key in p]
        if len(xs) < 2:
            return None
        return st.mean(xs) + 2 * st.pstdev(xs)
    cal = _ms("het_rate")
    return {"het_pure": cal} if cal is not None else {}
