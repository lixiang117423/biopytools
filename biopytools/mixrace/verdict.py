"""mixrace 三分支判读|mixrace three-branch verdict.

分支(导师 v4 报告口径):
  pure 纯菌(杂合率<阈值,可保存) / divergent 优势菌株·参考差异型(可保存) /
  contaminated 混杂菌株(强混合伴侣信号,需再分离纯化) / uncertain(位点不足)。
  Pb9-pattern 判据:存在样本B,A 的杂合位点上 B 携带ALT 比例≥阈值 且其中 B 纯合
  1/1 占比≥阈值(区分「B型+参考型混合」与「群内互为杂合」)。
  |pure (keep) / divergent (keep) / contaminated (re-isolate) / uncertain.
"""
from typing import List, Optional

import numpy as np


def judge(row: dict, partner_alt: np.ndarray, partner_hom: np.ndarray,
          samples: List[str], params: dict) -> dict:
    """三分支判读|three-branch verdict.

    Args:
        row: L1 统计行(sample/n_sites/n_het/het_rate/robust_rate/median_altfrac...)
        partner_alt/partner_hom: N×N 伴侣矩阵(het_eval.compute_shared_partner)
        samples: 与矩阵索引对应的样本名
        params: {pure_het, partner_alt_min, partner_hom_min, min_sites}
    """
    pure_het = params.get("pure_het", 0.001)
    alt_min = params.get("partner_alt_min", 0.8)
    hom_min = params.get("partner_hom_min", 0.5)
    min_sites = params.get("min_sites", 1000)
    i = samples.index(row["sample"])
    base = {"partner": None, "mix_proportion": None, "subtag": ""}

    if row.get("n_sites", 0) < min_sites:
        return {**base, "verdict": "uncertain",
                "advice": "位点不足,建议加深测序或检查数据|insufficient sites",
                "rationale": f"有GT位点仅 {row.get('n_sites', 0)} < {min_sites},"
                             f"统计不稳定|called sites below minimum"}
    if row["het_rate"] < pure_het:
        return {**base, "verdict": "pure",
                "advice": "建议: 可保存|Advice: safe to keep",
                "rationale": f"总杂合率 {row['het_rate']*100:.4f}% < {pure_het*100:.1f}%,"
                             f"单倍体背景下无混杂信号|het rate below purity threshold"}

    # 混合伴侣:取 ALT 携带率最高的样本 B|best partner by ALT-carrier rate
    alt_row = np.asarray(partner_alt[i], dtype=float)
    hom_row = np.asarray(partner_hom[i], dtype=float)
    if alt_row.size:
        j = int(np.argmax(alt_row))
        if alt_row[j] >= alt_min and hom_row[j] >= hom_min:
            partner = samples[j]
            mix = row.get("median_altfrac")
            mix = float(mix) if mix is not None and np.isfinite(mix) else None
            comp = (f"≈{mix*100:.0f}% {partner}型 + {(1-mix)*100:.0f}% 参考型"
                    if mix is not None else "成分比例未知(无杂合 altfrac)")
            return {**base, "verdict": "contaminated", "partner": partner,
                    "mix_proportion": mix,
                    "advice": "建议: 需再分离纯化|Advice: re-isolate to purify",
                    "rationale": (f"杂合位点 {alt_row[j]*100:.1f}% 被 {partner} 携带ALT"
                                  f"(其中 {hom_row[j]*100:.1f}% 为纯合1/1),"
                                  f"成分推断 {comp};符合 Pb9-Pb22 型混合模式"
                                  f"|partner-carrier pattern")}
    subtag = "轻度" if row.get("robust_rate", 1.0) < pure_het else ""
    note = "杂合多为低深度错误特征" if subtag else "存在真实杂合信号(见稳健杂合率)"
    return {**base, "verdict": "divergent", "subtag": subtag,
            "advice": "建议: 可保存;高精度下游可强制纯合化(取altfrac>0.5优势等位)"
                      "|Advice: keep; force-homozygose if needed",
            "rationale": (f"总杂合率 {row['het_rate']*100:.2f}% ≥ {pure_het*100:.1f}% "
                          f"但无强混合伴侣(最高携带率 "
                          f"{float(np.nanmax(alt_row))*100:.1f}% 或纯合占比不足),{note},"
                          f"判为优势菌株/参考差异型|divergent-lineage pattern")}
