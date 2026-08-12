"""mixrace VAF / 多等位 / Fws 分析(纯逻辑)|mixrace VAF / multiallelic / Fws analysis (pure).

step05: 解析 AD → 计算 VAF → 中间频率占比 + 样本内杂合度 Hw。
step06: 多等位位点统计 + Fws(队列池化频率,单样本退 Hw)。
核心指标基于读计数,与倍性无关。|step05: AD→VAF→intermediate ratio + within-sample
heterozygosity Hw. step06: multiallelic stats + Fws (cohort-pooled; Hw fallback if single).
Read-count based, ploidy-independent.
"""
from typing import Dict, List, Optional

# VAF 中间频率区间(用户定义,与判读阈值 vaf_mid_low/high 是两回事)|
# VAF intermediate-frequency band (user-defined; distinct from verdict thresholds)
VAF_MID_LO = 0.15
VAF_MID_HI = 0.85


def parse_vcf_ad(text: str) -> List[dict]:
    """解析 bcftools query 输出(CHROM/POS/REF/ALT/AD)|parse bcftools query output.

    期望格式: -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n'
    每行: CHROM<TAB>POS<TAB>REF<TAB>ALT(逗号)<TAB>AD(ref,alt1,alt2,...)
    |Each line: CHROM<TAB>POS<TAB>REF<TAB>ALT(comma)<TAB>AD(ref,alt1,...)
    """
    recs = []
    for line in text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, alt_str, ad_str = f[0], int(f[1]), f[2], f[3], f[4]
        alts = alt_str.split(",")
        try:
            ads = [int(x) for x in ad_str.split(",") if x and x != "."]
        except ValueError:
            ads = []
        recs.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alts": alts,
            "ad_ref": ads[0] if ads else 0,
            "ad_alts": ads[1:],
            "is_multiallelic": len(alts) >= 2,
        })
    return recs


def compute_vaf(records: List[dict]) -> List[dict]:
    """计算每位点每 ALT 的 VAF = AD_alt / (AD_ref + ΣAD_alt)|per-ALT VAF."""
    out = []
    for r in records:
        tot = (r["ad_ref"] + sum(r["ad_alts"])) or 1
        vafs = [a / tot for a in r["ad_alts"]]
        out.append({**r, "dp": tot, "vafs": vafs})
    return out


def vaf_metrics(all_vafs: List[float], lo: float = VAF_MID_LO, hi: float = VAF_MID_HI) -> dict:
    """中间频率占比 + 样本内杂合度 Hw|intermediate-frequency ratio + within-sample het Hw.

    intermediate_ratio = |{v in [lo,hi]}| / total;Hw = mean(2·v·(1−v))。
    |intermediate_ratio = fraction of VAFs in [lo,hi]; Hw = mean(2·v·(1-v)).
    """
    vs = [v for v in all_vafs if v is not None]
    n = len(vs) or 1
    inter = sum(1 for v in vs if lo <= v <= hi)
    hw = sum(2 * v * (1 - v) for v in vs) / n
    return {
        "total_sites": len(vs),
        "intermediate": inter,
        "intermediate_ratio": inter / n,
        "hw": hw,
    }


def multiallelic_stats(records: List[dict]) -> dict:
    """统计 ≥3 等位(ALT 数≥2)位点|count sites with >=2 ALT alleles (>=3 total)."""
    total = len(records)
    c = sum(1 for r in records if r.get("is_multiallelic"))
    return {"count": c, "ratio": c / (total or 1), "total": total}


def compute_fws(per_sample_vafs: Dict[str, Dict[tuple, float]]) -> Dict[str, dict]:
    """计算 Fws(队列池化)|compute Fws (cohort-pooled).

    Args:
        per_sample_vafs: {sample: {(chrom,pos): p_alt}} 各样本主 ALT 频率|per-sample dominant ALT freq

    Returns:
        {sample: {"hw": float, "fws": float|None}};队列<2 样本时 fws=None(单样本无群体基线)
        |fws=None when cohort <2 (single sample has no population baseline)
    """
    samples = list(per_sample_vafs)
    sites_union = set()
    for pv in per_sample_vafs.values():
        sites_union |= set(pv)
    P = {}
    for site in sites_union:
        ps = [per_sample_vafs[s][site] for s in samples if site in per_sample_vafs[s]]
        P[site] = sum(ps) / len(ps)
    Hb = sum(2 * p * (1 - p) for p in P.values()) / (len(P) or 1)
    out = {}
    for s in samples:
        pv = per_sample_vafs[s]
        Hw = sum(2 * p * (1 - p) for p in pv.values()) / (len(pv) or 1)
        fws = (1 - Hw / Hb) if (len(samples) >= 2 and Hb > 0) else None
        out[s] = {"hw": Hw, "fws": fws}
    return out
