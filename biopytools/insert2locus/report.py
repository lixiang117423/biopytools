"""面向非专业人员的自包含HTML报告|Self-contained HTML report for non-specialists"""

import html as _html
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from .locus_builder import CompleteLocus, JunctionInfo
from .verifier import SegmentStat

_GRADE_LABEL = {
    "PASS": "重构成功,覆盖验证通过|Reconstructed and verified",
    "WARN": "部分完成,请看注意事项|Partially complete, see caveats",
    "FAIL": "未能重构完整插入位点|Complete locus not reconstructed",
}

_STYLE = """
*{box-sizing:border-box}
body{margin:0;background:#FBFCFB;color:#1B211C;font-family:"Noto Sans SC",
"PingFang SC","Microsoft YaHei",system-ui,sans-serif;font-size:16px;line-height:1.7}
.wrap{max-width:860px;margin:0 auto;padding:0 22px}
header{padding:44px 0 26px;border-bottom:1px solid #DDE3D9}
h1{font-size:26px;margin:0 0 10px}
.meta{color:#5C645B;font-size:13.5px;display:flex;flex-wrap:wrap;gap:6px 22px}
section{padding:30px 0;border-bottom:1px solid #DDE3D9}
h2{font-size:20px;margin:4px 0 12px}
table{width:100%;border-collapse:collapse;margin:12px 0;font-size:14px}
th,td{padding:8px 10px;border-bottom:1px solid #DDE3D9;text-align:left}
th{color:#5C645B;font-weight:500;font-size:12px;text-transform:uppercase}
.badge{display:inline-block;padding:3px 14px;border-radius:16px;color:#fff;
font-weight:700;font-size:14px}
.badge.PASS{background:#2E7D32}.badge.WARN{background:#B8860B}
.badge.FAIL{background:#B3402A}
.callout{border:1px solid #CFD9CB;border-left:4px solid #2E4A32;background:#F3F6F1;
padding:12px 16px;border-radius:6px;margin:14px 0}
.callout.warn{border-left-color:#B8860B;background:#F8F2DC;border-color:#E6D9B0}
.seq-card{border:1px solid #CFD9CB;border-radius:6px;margin:10px 0;background:#fff}
.seq-head{display:flex;align-items:center;gap:10px;padding:8px 12px;flex-wrap:wrap}
.seq-title{font-family:monospace;font-size:13px}
.seq-info{color:#5C645B;font-size:12.5px;flex:1}
.copy-btn{border:none;background:#2E4A32;color:#fff;padding:5px 12px;border-radius:5px;
cursor:pointer;font-size:13px}
.copy-btn.copied{background:#B8860B}
.seq-body{border-top:1px solid #CFD9CB;max-height:240px;overflow:auto}
.seq-body code{display:block;white-space:pre-wrap;word-break:break-all;font-size:11.5px;
line-height:1.45;padding:10px 12px;font-family:monospace}
details{margin:12px 0}
summary{cursor:pointer;font-weight:700;color:#2E4A32}
footer{padding:24px 0 50px;color:#5C645B;font-size:12.5px}
svg{width:100%;height:auto;display:block}
"""

_SCRIPT = """
function copySeq(btn){
var el=document.getElementById(btn.dataset.target);
if(!el){return;}
var text=el.textContent;
function done(){btn.textContent='已复制';btn.classList.add('copied');
setTimeout(function(){btn.textContent='复制序列';btn.classList.remove('copied')},1500);}
if(navigator.clipboard&&navigator.clipboard.writeText){
navigator.clipboard.writeText(text).then(done);
}else{
var ta=document.createElement('textarea');ta.value=text;
ta.style.position='fixed';ta.style.top='-1000px';document.body.appendChild(ta);
ta.select();try{document.execCommand('copy')}catch(e){}
document.body.removeChild(ta);done();}
}
"""


def _esc(s) -> str:
    """HTML转义|HTML escape"""
    return _html.escape(str(s), quote=True)


def backbone_warnings(record_depth: dict, anchor_record: Optional[str]) -> List[str]:
    """非锚定record有覆盖→backbone整合警示|Covered non-anchor records imply backbone"""
    warns = []
    for record, stats in (record_depth or {}).items():
        if record == anchor_record:
            continue
        if stats.get("covered_frac", 0) > 0:
            warns.append(
                f"载体序列 {record} (长度{stats.get('length', '?')}bp, "
                f"平均深度{stats.get('mean', 0)}x) 也有reads覆盖, "
                f"检测到载体骨架序列可能也整合进基因组, 建议另做骨架特异PCR确认|"
                f"Vector sequence {record} also has read coverage; backbone may be "
                f"integrated, confirm with backbone-specific PCR")
    return warns


def _diagram_svg(locus: Optional[CompleteLocus],
                 segments: List[SegmentStat]) -> str:
    """结论区结构图|Structure diagram for the summary card"""
    if locus is None:
        return ('<div class="callout warn">未获得完整插入位点序列, 无法绘制结构图|'
                'No complete locus reconstructed, diagram unavailable</div>')
    depths = {s.name: s.mean_depth for s in segments}
    parts = [("植物左侧翼 LB", locus.lead, "#E9EEE6", "#2E4A32",
              depths.get("LB", 0)),
             ("插入序列 insert", locus.insert_len, "#2E4A32", "#FBFCFB",
              depths.get("insert", 0)),
             ("植物右侧翼 RB", locus.trail, "#E9EEE6", "#2E4A32",
              depths.get("RB", 0))]
    total = sum(p[1] for p in parts) or 1
    view_w, x = 720.0, 10.0
    rects = []
    for label, length, fill, text_fill, depth in parts:
        w = max(view_w * length / total, 60.0)
        rects.append(
            f'<rect x="{x:.0f}" y="30" width="{w:.0f}" height="34" '
            f'fill="{fill}" stroke="#3F6444"/>'
            f'<text x="{x + w / 2:.0f}" y="51" text-anchor="middle" '
            f'font-size="11" fill="{text_fill}">'
            f'{_esc(label)} {length}bp</text>'
            f'<text x="{x + w / 2:.0f}" y="82" text-anchor="middle" '
            f'font-size="10.5" fill="#5C645B">~{depth}x</text>')
        x += w
    return (f'<svg viewBox="0 0 740 95" role="img" aria-label="插入位点结构">'
            f'{"".join(rects)}</svg>')


def _seq_card(title: str, info: str, seq: str, target_id: str) -> str:
    """序列卡片(复制按钮)|Sequence card with copy button"""
    return (
        f'<div class="seq-card"><div class="seq-head">'
        f'<span class="seq-title">{_esc(title)}</span>'
        f'<span class="seq-info">{_esc(info)}</span>'
        f'<button class="copy-btn" onclick="copySeq(this)" '
        f'data-target="{_esc(target_id)}">复制序列</button></div>'
        f'<div class="seq-body"><code id="{_esc(target_id)}">'
        f'{_esc(seq)}</code></div></div>')


def render_report(sample: str, locus: Optional[CompleteLocus],
                  segments: List[SegmentStat], grade: str,
                  lb_junction_reads: int, rb_junction_reads: int,
                  junctions: List[JunctionInfo], record_depth: dict,
                  warnings: List[str], out_dir, tool_versions: dict,
                  junction_seqs: Optional[dict] = None,
                  anchor_label: Optional[str] = None) -> Path:
    """渲染单样本报告,返回文件路径|Render per-sample report, return path"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{sample}.insert2locus.report.html"

    # backbone判定锚:优先region模式标签,退回locus锚定record名|
    # Anchor for backbone check: region label first, else locus anchor record
    anchor = anchor_label or (locus.insert_record if locus else None)
    bb_warns = backbone_warnings(record_depth, anchor)
    all_warnings = list(warnings) + bb_warns

    # 结论区|Summary
    if locus is None:
        summary_html = (
            f'<div class="callout warn">未能重构出贯穿"植物-插入-植物"的完整序列。'
            f'可能原因:插入拷贝数异常/测序深度不足/边界重复。可检查最佳junction '
            f'contigs(见下方)或加大测序数据量后重跑。|Failed to reconstruct a '
            f'plant-insert-plant spanning contig. Check best junction contigs below '
            f'or rerun with more data.</div>')
    else:
        summary_html = _diagram_svg(locus, segments)

    # 结果表|Results table
    if segments:
        rows = []
        source = {"LB": "植物DNA|Plant DNA", "insert": "插入序列|Inserted DNA",
                  "RB": "植物DNA|Plant DNA"}
        for s in segments:
            cont = ("✓ 连续|Continuous" if s.zero_windows == 0
                    else f"⚠ 缺口{s.zero_windows}处|{s.zero_windows} gap(s)")
            rows.append(
                f"<tr><td>{_esc(s.name)}</td><td>{s.start}–{s.end}</td>"
                f"<td>{s.length} bp</td><td>{_esc(source.get(s.name, ''))}</td>"
                f"<td>~{s.mean_depth}x</td><td>{cont}</td></tr>")
        table_html = (
            "<table><tr><th>区段</th><th>位置</th><th>长度</th><th>来源</th>"
            f"<th>平均深度</th><th>覆盖连续性</th></tr>{''.join(rows)}</table>")
    else:
        table_html = "<p>无覆盖统计数据|No coverage statistics available</p>"

    # 序列卡片|Sequence cards
    cards = []
    if locus:
        cards.append(_seq_card(
            f"{sample}_LB_left_flank", f"{locus.lead} bp 左边界植物侧翼|Left flank",
            locus.seq[:locus.lead], "seq_lb"))
        cards.append(_seq_card(
            f"{sample}_RB_right_flank", f"{locus.trail} bp 右边界植物侧翼|Right flank",
            locus.seq[len(locus.seq) - locus.trail:], "seq_rb"))
        cards.append(_seq_card(
            f"{sample}_complete_locus",
            f"{len(locus.seq)} bp 左翼+插入+右翼|Full locus",
            locus.seq, "seq_full"))
    else:
        junction_seqs = junction_seqs or {}
        for i, j in enumerate(junctions[:10]):
            seq = junction_seqs.get(j.contig, "")
            cards.append(_seq_card(
                j.contig,
                f"border={j.border} flank={j.flank_bp}bp "
                f"{'junction contig(序列见junction_contigs.fasta)' if not seq else 'junction contig'}",
                seq, f"junc_{i}"))
    cards_html = ("".join(cards) if cards
                  else "<p>无可复制序列|No sequences to copy</p>")

    # 引物指引|Primer guidance
    if locus:
        primer_html = (
            f"<ul><li>左边界LB:植物端引物落在第 <b>1–{locus.lead}</b> "
            f"(左侧翼),插入端引物落在插入段(第{locus.insert_start}bp起)</li>"
            f"<li>右边界RB:植物端引物落在第 "
            f"<b>{locus.insert_end + 1}–{len(locus.seq)}</b> (右侧翼),"
            f"插入端引物落在插入段(第{locus.insert_end}bp止)</li>"
            f"<li>引物常规要求:18–25bp,Tm 55–65℃成对接近,"
            f"GC 40–60%,3′端避免连续单一碱基;建议Primer3设计+"
            f"BLAST查特异性</li>"
            f"<li>设计前请对照插入序列已知方向,确认引物取向正确</li></ul>")
    else:
        primer_html = ("<p>无完整位点,暂无法给出引物区间;可基于上方junction "
                       "contigs的侧翼段设计|No locus available for primer ranges; "
                       "use junction contig flanks instead</p>")

    # 注意事项|Caveats
    notes = []
    if locus and (locus.lead < 50 or locus.trail < 50):
        notes.append("侧翼较短(<50bp),放引物可能偏短,建议再做一轮步移拉长|"
                     "Flanks are short (<50bp); another walking run may extend them")
    for s in segments:
        if s.zero_windows > 0:
            notes.append(f"{s.name} 段存在{s.zero_windows}处零覆盖缺口"
                         f"(位置{s.start}–{s.end}内),拼接可能在缺口处有误|"
                         f"{s.name} has {s.zero_windows} zero-coverage gap(s); "
                         f"assembly may be wrong at gaps")
    if locus and (lb_junction_reads == 0 or rb_junction_reads == 0):
        notes.append("某一边界缺少跨边界reads支持,该junction证据不足|"
                     "One boundary lacks spanning-read support")
    notes.extend(all_warnings)
    notes_html = ("<ul>" + "".join(f"<li>{_esc(n)}</li>" for n in notes)
                  + "</ul>") if notes else "<p>无特别注意事项|No special caveats</p>"

    # 工具版本表|Tool versions
    ver_rows = "".join(
        f"<tr><td>{_esc(k)}</td><td>{_esc(v)}</td></tr>"
        for k, v in (tool_versions or {}).items())

    html_doc = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_esc(sample)} 转基因插入位点分析报告</title>
<style>{_STYLE}</style></head><body>
<header class="wrap"><h1>{_esc(sample)} 转基因插入位点分析报告</h1>
<div class="meta"><span>样本|Sample: <b>{_esc(sample)}</b></span>
<span>评级|Grade: <span class="badge {_esc(grade)}">{_esc(grade)}</span></span>
<span>{_esc(_GRADE_LABEL.get(grade, grade))}</span>
<span>日期|Date: {_esc(datetime.now().strftime("%Y-%m-%d"))}</span></div>
</header><main class="wrap">
<section><h2>结论|Summary</h2>
{summary_html}
<p>跨边界reads支持|Junction-spanning reads: LB边界 {lb_junction_reads} 条,
RB边界 {rb_junction_reads} 条。</p>
</section>
<section><h2>结果明细|Results</h2>{table_html}</section>
<section><h2>序列(可一键复制)|Sequences (copy with one click)</h2>{cards_html}</section>
<section><h2>引物设计指引|Primer design guidance</h2>{primer_html}</section>
<section><h2>注意事项|Caveats</h2>
<div class="callout{' warn' if notes else ''}">{notes_html}</div></section>
<details><summary>方法与工具版本|Methods &amp; tool versions</summary>
<p>方法:以插入序列(载体+片段)为参考比对WGS reads,提取跨边界soft-clip reads
与侧翼候选,迭代步移局部组装(SPAdes)延伸侧翼,招募全部覆盖插入序列的reads
重构完整位点,最后将WGS reads比回重构序列检查覆盖连续性。|
Method: map WGS reads to the insert reference, fish junction soft-clip reads,
extend flanks by iterative local assembly (SPAdes), recruit all insert-covering
reads to rebuild the complete locus, then verify end-to-end coverage.</p>
<table><tr><th>工具|Tool</th><th>版本|Version</th></tr>{ver_rows}</table>
</details>
</main><footer class="wrap">由 biopytools insert2locus 生成|Generated by
biopytools insert2locus</footer>
<script>{_SCRIPT}</script></body></html>"""
    out_path.write_text(html_doc)
    return out_path


def render_index(samples: List[dict], out_path) -> Path:
    """多样本导航页|Multi-sample navigation page"""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in samples:
        grade = s.get("grade", "?")
        err = (f'<br><span style="color:#B3402A">{_esc(s["error"])}</span>'
               if s.get("error") else "")
        link = (f'<a href="{_esc(s["report"])}">{_esc(s["sample"])} 报告</a>'
                if s.get("report") else _esc(s["sample"]))
        rows.append(
            f'<tr><td>{link}</td>'
            f'<td><span class="badge {_esc(grade)}">{_esc(grade)}</span></td>'
            f"<td>{_esc(s.get('stop_reason', ''))}{err}</td></tr>")
    html_doc = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<title>insert2locus 样本导航|Sample index</title>
<style>{_STYLE}</style></head><body>
<header class="wrap"><h1>insert2locus 分析结果导航|Analysis index</h1></header>
<main class="wrap"><section>
<table><tr><th>样本|Sample</th><th>评级|Grade</th><th>说明|Notes</th></tr>
{''.join(rows)}</table></section></main>
<footer class="wrap">由 biopytools insert2locus 生成|Generated by
biopytools insert2locus</footer></body></html>"""
    out_path.write_text(html_doc)
    return out_path
