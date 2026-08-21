"""mixrace 报告层(v0.3 三分支)|mixrace reporting (three-branch verdict).

判读汇总表(TSV/HTML)、单样本 md 报告(证据链+实验建议)、自包含 HTML(关键图内嵌)。
|Summary tables, per-sample evidence-chain reports, self-contained HTML.
"""
import base64
import html as _html
from typing import Tuple

_SUMMARY_COLS = ["sample", "verdict", "advice", "het_rate", "robust_rate",
                 "shared_only_rate", "top_partner", "mix_proportion",
                 "het_rate_after_hotspot", "dp_ratio", "host_rate",
                 "pathogen_map_rate", "contamination_rate", "mean_depth", "breadth_1x"]

_VERDICT_CN = {"pure": "纯菌", "divergent": "优势菌株/参考差异型",
               "contaminated": "混杂菌株", "uncertain": "不确定"}

_SUBTAG_CN = {"mild": "轻度"}   # subtag 枚举值→显示名|subtag enum -> display label

_METRIC_EXPLAIN = {
    "het_rate": "总杂合率:单倍体每个位点本该纯合,出现杂合=混合或错误。这是判读的主指标,像混合比例的指纹。",
    "robust_rate": "稳健杂合率:只统计 altAD>=5 且 altfrac>=0.2 的杂合位点(排除低深度测序错误)后的杂合率。",
    "shared_only_rate": "shared-only 杂合率:只统计 ALT 在其他样品也出现的杂合位点(排除样品特异噪声)。",
    "top_partner": "混合伴侣:在我的杂合位点上携带 ALT 最多的样品——相当于在群体里找与你共享变异的『另一半』。",
    "mix_proportion": "成分推断:若判混杂,样品约等于 该比例x伴侣型 + 余量x参考型 的混合。",
    "het_rate_after_hotspot": "排除热点后杂合率:剔除所有混杂样品共享的高差异窗口后剩下的杂合率,仍高则混杂证据坚实。",
    "dp_ratio": "DP检验:杂合位点深度/纯合位点深度。>1.5 说明杂合位点深度反而更高,是真实混合信号而非低覆盖错误。",
    "host_rate": "寄主reads占比:比对上寄主基因组的reads比例,代表样品里植物DNA污染程度。",
    "pathogen_map_rate": "病原mapping率:剔除寄主后比对上病原基因组的reads比例,代表病原数据质量。",
    "contamination_rate": "污染reads占比:既没比对上寄主、也没比对上病原基因组的reads比例,可能来自其他微生物或低质量reads。",
    "mean_depth": "平均测序深度:每个位置平均被测了多少遍。",
    "breadth_1x": "覆盖广度:基因组至少被测到1次的碱基占比;广度低说明有大片区域没测到。",
}

_EVIDENCE_KEYS = ["het_rate", "robust_rate", "shared_only_rate", "top_partner",
                  "mix_proportion", "het_rate_after_hotspot", "dp_ratio",
                  "host_rate", "pathogen_map_rate", "contamination_rate",
                  "mean_depth", "breadth_1x"]

_KEY_FIGURES = ["het_heatmap_100kb", "pca_3view", "nj_tree", "eval_3panel", "manhattan_grid"]


def _html_escape(text: str) -> str:
    return _html.escape(str(text))


def _fmt(key, val) -> str:
    """指标值格式化(None→—;比率→百分比)|format metric value."""
    if val is None or val == "":
        return "—"
    if isinstance(val, str):
        return val
    if key in ("het_rate", "robust_rate", "shared_only_rate", "het_rate_after_hotspot"):
        return f"{val*100:.4f}%" if val < 0.01 else f"{val*100:.2f}%"
    if key in ("host_rate", "pathogen_map_rate", "contamination_rate"):
        return f"{val*100:.2f}%"
    if key == "mix_proportion":
        return f"{val*100:.0f}%"
    if key == "dp_ratio":
        return f"{val:.2f}"
    if key == "mean_depth":
        return f"{val:.1f}x"
    if key == "breadth_1x":
        return f"{val:.2f}%"
    return str(val)


def build_summary_table(rows: list) -> Tuple[str, str]:
    """判读汇总表(TSV+HTML)|verdict summary table (tsv + html)."""
    disp = []
    for r in rows:
        d = dict(r)
        for k in ("het_rate", "robust_rate", "shared_only_rate",
                  "het_rate_after_hotspot", "host_rate", "pathogen_map_rate",
                  "contamination_rate", "mix_proportion", "dp_ratio",
                  "mean_depth", "breadth_1x"):
            if k in d:
                d[k] = _fmt(k, d[k])
        d["verdict"] = f"{_VERDICT_CN.get(r.get('verdict'), r.get('verdict', ''))}" \
                       f"{_SUBTAG_CN.get(r.get('subtag', ''), r.get('subtag', ''))}"
        disp.append(d)
    tsv = "\t".join(_SUMMARY_COLS) + "\n"
    for d in disp:
        tsv += "\t".join(str(d.get(c, "")) for c in _SUMMARY_COLS) + "\n"
    v_idx = _SUMMARY_COLS.index("verdict")
    _COLS_CN = {"sample": "样品", "verdict": "判读", "advice": "建议",
                "het_rate": "总杂合率", "robust_rate": "稳健杂合率",
                "shared_only_rate": "shared-only杂合率", "top_partner": "混合伴侣",
                "mix_proportion": "成分推断", "het_rate_after_hotspot": "排除热点后杂合率",
                "dp_ratio": "DP检验", "host_rate": "寄主占比",
                "pathogen_map_rate": "病原mapping率", "contamination_rate": "污染reads",
                "mean_depth": "平均深度", "breadth_1x": "覆盖广度"}
    html = '<table border=1><tr>' + "".join(
        f"<th>{_COLS_CN.get(c, c)}<br>{c}</th>" for c in _SUMMARY_COLS) + "</tr>"
    for d, r in zip(disp, rows):
        color = {"pure": "#2e7d32", "divergent": "#ef6c00",
                 "contaminated": "#c62828", "uncertain": "#9e9e9e"}.get(
                     str(r.get("verdict", "")), "#333")
        cells = []
        for ci, c in enumerate(_SUMMARY_COLS):
            style = f' style="color:{color};font-weight:bold"' if ci == v_idx else ""
            cells.append(f"<td{style}>{_html_escape(str(d.get(c, '')))}</td>")
        html += "<tr>" + "".join(cells) + "</tr>"
    html += "</table>"
    return tsv, html


def build_sample_report(sample: str, row: dict, figures: dict) -> str:
    """单样本 md 报告(证据链+建议+图)|per-sample markdown (evidence + advice + figures)."""
    verdict = str(row.get("verdict", "uncertain"))
    md = [f"# {sample} 混杂评估报告|{sample} contamination report", ""]
    md.append(f"**判读|Verdict: {_VERDICT_CN.get(verdict, verdict)}"
              f"{_SUBTAG_CN.get(row.get('subtag', ''), row.get('subtag', '') or '')}**")
    md.append("")
    md.append(f"**{row.get('advice', '')}**")
    md.append("")
    md.append("## 依据|Rationale")
    md.append("")
    md.append(str(row.get("rationale", "")))
    md.append("")
    md.append("## 证据链|Evidence")
    md.append("")
    md.append("| 指标<br>Metric | 值<br>Value | 通俗解释<br>Plain words |")
    md.append("|---|---|---|")
    for k in _EVIDENCE_KEYS:
        if k in row:
            md.append(f"| {k} | {_fmt(k, row.get(k))} | {_METRIC_EXPLAIN.get(k, '')} |")
    md.append("")
    if figures:
        md.append("## 图|Figures")
        md.append("")
        for stem, p in figures.items():
            md.append(f"![{stem}]({p})")
            md.append("")
    return "\n".join(md)


def _embed_image(path: str):
    """PNG → base64(单文件HTML内嵌)|PNG to base64 for self-contained HTML."""
    try:
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode("ascii")
    except OSError:
        return None


def build_html_report(title: str, rows: list, figures: dict,
                         verdict_note: str = "") -> str:
    """自包含 HTML(判读表+关键图内嵌+逐样品证据链折叠)|self-contained HTML report.

    verdict_note: 判读口径文案(阈值可配,由调用方按 config 生成防漂移)
    |verdict_note: threshold note built from config by the caller (no drift).
    """
    _, table_html = build_summary_table(rows)
    parts = [f"<html><head><meta charset='utf-8'><title>{_html_escape(title)}</title>",
             "<style>body{font-family:sans-serif;margin:24px}"
             "table{border-collapse:collapse;width:100%;margin:12px 0;font-size:0.9em}"
             "td,th{border:1px solid #999;padding:4px 8px;text-align:left}"
             "details{margin:8px 0}img{max-width:100%}"
             "</style></head><body>",
             f"<h1>{_html_escape(title)}</h1>"]
    if verdict_note:
        parts.append(f"<p>{_html_escape(verdict_note)}</p>")
    parts.append(table_html)
    for stem in _KEY_FIGURES:
        b64 = _embed_image(figures.get(stem, "")) if figures.get(stem) else None
        if b64:
            parts.append(f"<h2>{_html_escape(stem)}</h2>"
                         f"<img src='data:image/png;base64,{b64}'>")
    parts.append("<h2>逐样品证据链|Per-sample evidence</h2>")
    for r in rows:
        verdict = str(r.get("verdict", "uncertain"))
        parts.append(f"<details><summary><b>{_html_escape(str(r.get('sample', '')))}</b> "
                     f"— {_VERDICT_CN.get(verdict, verdict)} "
                     f"{_html_escape(str(r.get('subtag', '') or ''))}</summary>"
                     f"<p>{_html_escape(str(r.get('advice', '')))}</p>"
                     f"<p>{_html_escape(str(r.get('rationale', '')))}</p>")
        parts.append("<table border=1><tr><th>指标</th><th>值</th><th>解释</th></tr>")
        for k in _EVIDENCE_KEYS:
            if k in r:
                parts.append(f"<tr><td>{k}</td><td>{_html_escape(_fmt(k, r.get(k)))}</td>"
                             f"<td>{_html_escape(_METRIC_EXPLAIN.get(k, ''))}</td></tr>")
        parts.append("</table></details>")
    parts.append("</body></html>")
    return "".join(parts)
