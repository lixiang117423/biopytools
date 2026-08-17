"""mixrace 报告生成|mixrace report generation.

step07 判读配套:VAF 直方图 R 脚本生成、per-sample md 报告、跨样品对比表(tsv+html)、
合并自包含 HTML 报告(图片 base64 内嵌)、genomescope summary 解析。
|Verdict-stage reporting: VAF histogram R script, per-sample markdown report,
cross-sample comparison table (tsv+html), merged self-contained HTML, genomescope parsing.
"""
from html import escape as _html_escape
from pathlib import Path
from typing import Tuple

_SUMMARY_COLS = ["sample", "verdict", "confidence", "het_rate", "afs_shape",
                 "dominant_proportion", "mean_depth"]


def parse_genomescope_summary(text: str) -> dict:
    """解析 genomescope2 summary.txt(基因组大小/杂合度)|parse summary.txt.

    容错 tab/:/= 与空格分隔;归一化 genome_size 键。|Tolerates tab/:/= and whitespace;
    normalizes a genome_size key.
    """
    s = {}
    if not text:
        return s
    for line in text.splitlines():
        k = v = None
        for sep in ("\t", ":", "="):
            if sep in line:
                k, v = line.split(sep, 1)
                break
        if k is None:
            parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            k, v = parts[0], parts[1]
        v = v.replace("%%", "").replace("%", "").replace(",", "").strip()
        try:
            s[k.strip()] = float(v)
        except ValueError:
            pass
    for k, val in list(s.items()):
        kl = k.lower()
        if ("haploid length" in kl or "genome size" in kl) and "genome_size" not in s:
            s["genome_size"] = val
    return s


def generate_vaf_histogram_r(vaf_tsv: str, png_out: str, rscript_path: str) -> str:
    """生成 VAF 直方图 R 脚本(返回脚本路径,由 main 用 build_conda_command 执行)
    |generate the VAF-histogram R script; main runs it via build_conda_command.

    vaf_tsv 须含 vafs 列(逗号分隔的 VAF)|vaf_tsv must have a vafs column (comma-joined VAFs).
    """
    script_path = str(Path(png_out).with_suffix(".R"))
    rcode = f'''#!/usr/bin/env Rscript
library(ggplot2)
d <- read.delim("{vaf_tsv}")
vafs <- as.numeric(unlist(strsplit(as.character(d$vafs), ",")))
vafs <- vafs[!is.na(vafs)]
ggplot(data.frame(vafs = vafs), aes(vafs)) +
  geom_histogram(binwidth = 0.02, fill = "#4C78A8", color = "white") +
  xlim(0, 1) +
  labs(x = "VAF", y = "Site count", title = "VAF distribution") +
  theme_minimal()
ggsave("{png_out}", width = 6, height = 4)
'''
    Path(script_path).write_text(rcode)
    return script_path


def build_sample_report(sample: str, metrics: dict, verdict: dict, paths: dict) -> str:
    """生成单样品 markdown 报告|build per-sample markdown report.

    Args:
        sample: 样本名|sample name
        metrics: 指标 dict(het_rate/het_sites/afs_shape/dominant_proportion/maf/mean_depth/heterozygosity?)
        verdict: judge() 返回(verdict/confidence/rationale)
        paths: {label: 相对图路径} 嵌入报告|{label: relative image path}
    """
    md = [f"# {sample} 混合小种检测报告|{sample} mixed-race report", ""]
    md.append(f"**判读|Verdict: {verdict['verdict']}**"
              f"（置信|confidence: {verdict.get('confidence', 'n/a')}）")
    md.append("")
    md.append("## 依据|Rationale")
    md.append("")
    md.append(verdict.get("rationale", ""))
    md.append("")
    md.append("## 指标|Metrics")
    md.append("")
    md.append("| 指标|Metric | 值|Value |")
    md.append("|---|---|")
    for k, v in metrics.items():
        md.append(f"| {k} | {_fmt(k, v)} |")
    md.append("")
    if paths:
        md.append("## 图|Figures")
        md.append("")
        for label, p in paths.items():
            md.append(f"![{label}]({p})")
            md.append("")
    return "\n".join(md)


def build_summary_table(rows: list) -> Tuple[str, str]:
    """跨样品对比表|build cross-sample comparison table.

    Returns:
        (tsv_text, html_text)
    """
    tsv = "\t".join(_SUMMARY_COLS) + "\n"
    for r in rows:
        tsv += "\t".join(str(r.get(c, "")) for c in _SUMMARY_COLS) + "\n"
    html = '<table border=1><tr>' + "".join(f"<th>{c}</th>" for c in _SUMMARY_COLS) + "</tr>"
    for r in rows:
        # 单元格转义防 HTML 注入(样名/文本含 <>&)|escape cells (names/text may contain <>&)
        html += "<tr>" + "".join(
            f"<td>{_html_escape(str(r.get(c, '')))}</td>" for c in _SUMMARY_COLS) + "</tr>"
    html += "</table>"
    return tsv, html


# ---------- 合并 HTML 报告(单文件,图片内嵌)|merged self-contained HTML report ----------

_VERDICT_STYLE = {
    "single_genotype": ("疑似纯 · 单一基因型", "#2e7d32"),
    "mixed_genotype": ("疑似混合 · 多基因型", "#c62828"),
    "uncertain": ("不确定 · 需进一步排查", "#f57f17"),
}

# 指标通俗解释(给无生信基础的研究员)|plain-language metric explanations
_METRIC_EXPLAIN = {
    "het_rate": "全基因组杂合率:同时测到两种碱基的位点占整个基因组的比例。根肿菌是单倍体(每个位点本该只有一种碱基),所以越低越纯。<0.01% 基本纯,0.1–1% 可疑,>1% 明显不纯。",
    "het_sites": "出现两种碱基共存的位点个数(即上面杂合率的绝对计数)。",
    "afs_shape": "等位频率谱形态:把所有变异位点的次要碱基频率画成分布图的形状。单峰贴在两头=纯;中间有峰=存在混合。",
    "dominant_proportion": "优势株占比:若为混合,最主要菌株大约占群体的百分比(基于各位点主碱基频率)。仅判为混合/可疑时显示;纯样品杂合位点太少(多为噪声),该占比无统计意义,显示为—。",
    "maf": "次等位频率(MAF):次要碱基的平均出现频率。纯样品接近 0;混合样品明显升高。",
    "mean_depth": "平均测序深度:每个位置平均被测了多少遍。≥50x 结论可靠;<50x 需谨慎(可能漏检低频信号)。",
    "heterozygosity": "GenomeScope k-mer 杂合度:独立于比对的另一视角,从 k-mer 频谱估计的基因组杂合度。",
    "genome_size": "GenomeScope 估计的基因组大小(bp)。",
    "ploidy": "Smudgeplot 推断的倍性(根肿菌预期单倍体=1)。",
}


def _embed_image(path: str):
    """读 PNG → base64(嵌入 HTML,单文件可分享)|read PNG as base64 for inline embed."""
    import base64
    try:
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode("ascii")
    except OSError:
        return None


# 判读中文短标签(用于树叶名)|verdict short label for tree tips
_TREE_VERDICT_CN = {"single_genotype": "纯", "mixed_genotype": "混合",
                    "uncertain": "不确定"}


def build_annotations_tsv(annotations: dict) -> str:
    """生成树注解 TSV 文本(sample/verdict/het;het 为空写 '-',避免叶标签尾空格)。
    |build annotation TSV text (het empty -> '-', avoids trailing space in tip label)."""
    lines = ["sample\tverdict\thet"]
    for sample in sorted(annotations):
        a = annotations[sample]
        cn = _TREE_VERDICT_CN.get(a.get("verdict", ""), "")
        het = a.get("het_rate")
        het_s = f"{het*100:.4f}%" if het is not None else "-"
        lines.append(f"{sample}\t{cn}\t{het_s}")
    return "\n".join(lines) + "\n"


def generate_tree_png_r(nwk_path: str, png_path: str, annotations: dict) -> str:
    """生成 ggtree 静态树 PNG 的 R 脚本;叶标签 = 样品编号 + 判读 + 杂合率。
    |generate ggtree static-tree R script; tip label = sample + verdict + het rate.

    annotations: {sample: {"verdict": str, "het_rate": float}}
    注解写到 <png>.ann.tsv,R 里按样名 join(ape read.tree 已解析引号样名)。
    |annotations written to <png>.ann.tsv; joined by tip label in R.
    """
    script_path = str(Path(png_path).with_suffix(".R"))
    ann_path = str(Path(png_path).with_suffix(".ann.tsv"))
    Path(ann_path).write_text(build_annotations_tsv(annotations), encoding="utf-8")
    rcode = f'''#!/usr/bin/env Rscript
try(Sys.setlocale("LC_ALL", "C.UTF-8"), silent = TRUE)   # C locale 下中文防乱码|avoid mojibake
library(ape)
library(ggtree)
library(ggplot2)
tree <- read.tree("{nwk_path}")
# 叶标签 = 编号[判读]杂合率(由 ann.tsv join,见下)。严禁对 x 轴设上限:
# 树深是 root→tip 各段枝长之和,天然就是 max_edge 的数倍,任何 max_edge*N 的上限都
# 可能小于"树深+offset+文字宽度";ggplot 的 xlim 是硬性丢弃超范围数据(实测整树
# Removed 25 rows,标签全灭)。不设 xlim 时 geom_tiplab 的 x 参与坐标训练,面板
# 自动扩到标签处。offset 用相对量(绝对 0.003 会比 ~1e-4 的枝长宽几十倍)。
# |tips carry sample[verdict]het joined from ann.tsv. NEVER cap the x-axis: root-to-tip depth is a path
# SUM, several times max_edge, so any max_edge*N cap can fall below depth+offset+
# text width; ggplot xlim hard-drops out-of-range rows (observed: Removed 25 rows,
# all labels gone). Without xlim the tiplab x trains the scale and the panel
# auto-extends. offset stays relative (absolute 0.003 dwarfed ~1e-4 edges).
max_edge <- max(tree$edge.length, na.rm = TRUE)
ann <- read.delim("{ann_path}", sep = "\t", header = TRUE, colClasses = "character",
                  encoding = "UTF-8", check.names = FALSE)
# 叶标签 = 编号[判读]杂合率(het 为 '-' 只留判读,避免尾空格)|tip = id[verdict]het
# |tip label = id[verdict]het ('-' het keeps no trailing space)
ann$tiplab <- ifelse(ann$het == "-", paste0(ann$sample, " [", ann$verdict, "]"),
                     paste0(ann$sample, " [", ann$verdict, "] ", ann$het))
p <- ggtree(tree) %<+% ann
p <- p + geom_tiplab(aes(label = tiplab), size = 3.2,
                     offset = max_edge * 0.25, align = TRUE)
has_support <- !is.null(tree$node.label) && any(nzchar(tree$node.label))
if (has_support) {{
  p <- p + geom_text(aes(label = label),
                     data = function(d) d[!d$isTip & d$label != "", , drop = FALSE],
                     hjust = -0.3, size = 2.5, color = "#777777")
}}
p <- p + labs(title = "样品聚类(系统发育树)|Sample clustering") +
  theme_tree2() + theme(plot.title = element_text(hjust = 0.5))
ggsave("{png_path}", p, width = 14, height = max(5, 0.3 * ape::Ntip(tree)),
       dpi = 150, limitsize = FALSE, device = png, type = "cairo")
'''
    Path(script_path).write_text(rcode, encoding="utf-8")
    return script_path


def _fmt(key, val):
    """指标值格式化(百分比/小数)|format metric value."""
    if val is None:
        return "—"
    if key in ("het_rate", "heterozygosity"):
        return f"{val*100:.4f}%" if val < 0.01 else f"{val*100:.2f}%"
    if key in ("dominant_proportion", "maf"):
        return f"{val*100:.1f}%"
    if key == "mean_depth":
        return f"{val:.1f}x"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


def _shape_cn(shape: str) -> str:
    return {
        "monoclonal": "单克隆(纯):频率集中在 0 或 1,几乎无中间值",
        "two_clone_50_50": "两个基因型约各占一半:中间频率集中在 0.5 附近",
        "dominant_minor": "一个优势株 + 少量其他:主峰偏高,伴次要峰",
        "smeared": "多个基因型混杂:频率从 0 到 1 连续分布,无明显主峰",
    }.get(shape, shape)


def build_html_report(title: str, samples_data: list, tree_png: str = None) -> str:
    """生成合并的自包含 HTML 报告(所有样品、图片 base64 内嵌、中文通俗解释;可选进化树)。
    |build merged self-contained HTML (all samples, base64 images, Chinese explanations;
    optional tree)."""
    rows_html = []
    for s in samples_data:
        v = s.get("verdict", "uncertain")
        label, color = _VERDICT_STYLE.get(v, _VERDICT_STYLE["uncertain"])
        rows_html.append(
            f'<tr><td>{_html_escape(str(s["sample"]))}</td>'
            f'<td style="color:{color};font-weight:bold">{label}</td>'
            f'<td>{s.get("confidence","")}</td><td>{_fmt("het_rate", s.get("het_rate"))}</td>'
            f'<td>{_shape_cn(s.get("afs_shape",""))}</td>'
            f'<td>{_fmt("dominant_proportion", s.get("dominant_proportion"))}</td>'
            f'<td>{_fmt("mean_depth", s.get("mean_depth"))}</td></tr>')

    parts = []
    parts.append("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
    parts.append(f"<title>{title}</title><style>{_CSS()}</style></head><body>")
    parts.append(f"<h1>{title}</h1>")
    # 汇总表|summary table
    parts.append("<h2>一、各样品判读汇总</h2>")
    parts.append("<table class='sum'><tr><th>样品</th><th>判读结论</th><th>置信</th>"
                 "<th>杂合率</th><th>频率谱形态</th><th>优势株占比</th><th>平均深度</th></tr>")
    parts.extend(rows_html)
    parts.append("</table>")

    # 指标说明|metric legend
    parts.append("<h2>二、指标含义(通俗版)</h2><dl class='legend'>")
    for k, expl in _METRIC_EXPLAIN.items():
        parts.append(f"<dt>{k}</dt><dd>{expl}</dd>")
    parts.append("</dl>")

    # 样品聚类树(可选,静态 PNG 内嵌;有树才占"三")|clustering tree (optional static PNG)
    has_tree = False
    if tree_png:
        b64 = _embed_image(tree_png)
        if b64:
            parts.append("<h2>三、样品聚类(系统发育树)</h2>")
            parts.append("<p class='muted'>用全部样品的 SNP 变异构建的聚类树:分支越近的编号,"
                         "基因组越相似。叶名格式为 样品编号[判读]杂合率(如 Pb5 [纯] 0.0005%);"
                         "节点数值为统计支持值(越高越可靠)。注意:树基于各位置的主碱基型,"
                         "混合样品可能表现为较长的分支。</p>")
            parts.append(f"<div class='tree'><img src='data:image/png;base64,{b64}'></div>")
            has_tree = True

    # 各样品详情(折叠,点击展开;动态编号防跳号)|per-sample details (dynamic numbering)
    details_num = "四" if has_tree else "三"
    parts.append(f"<h2>{details_num}、各样品详细结果（点击样品行展开/收起）</h2>")
    parts.append("<div class='toggle-bar'>"
                 "<button onclick=\"document.querySelectorAll('details.sample').forEach(d=>d.open=true)\">展开全部</button> "
                 "<button onclick=\"document.querySelectorAll('details.sample').forEach(d=>d.open=false)\">收起全部</button>"
                 "</div>")
    for s in samples_data:
        parts.append(_sample_section(s))
    parts.append("</body></html>")
    return "\n".join(parts)


def _sample_section(s: dict) -> str:
    sample = _html_escape(str(s.get("sample", "")))
    v = s.get("verdict", "uncertain")
    label, color = _VERDICT_STYLE.get(v, _VERDICT_STYLE["uncertain"])
    m = s.get("metrics", {})
    dom = s.get("dominant_proportion")
    # rationale 含 < 等字符(如阈值 <0.010%),必须转义否则被浏览器当标签吞掉|
    # rationale contains < (e.g. threshold <0.010%); escape or the browser eats it as a tag.
    rationale = _html_escape(str(s.get("rationale", "")))
    # summary 行(始终可见,点击展开/收起)|always-visible clickable summary row
    sumrow = (
        f"<summary class='sumrow' style='border-left:6px solid {color}'>"
        f"<span class='sname'>{sample}</span>"
        f"<span class='sverdict' style='color:{color}'>● {label}</span>"
        f"<span class='smeta'>杂合率 {_fmt('het_rate', s.get('het_rate'))}</span>"
        f"<span class='smeta'>优势株占比 {_fmt('dominant_proportion', dom)}</span>"
        f"<span class='smeta'>平均深度 {_fmt('mean_depth', s.get('mean_depth'))}</span>"
        f"<span class='conf'>置信度:{s.get('confidence','')}</span>"
        f"</summary>")
    # 展开内容:依据 + 指标表 + 图|expanded body: rationale + metrics table + figures
    order = ["het_rate", "het_sites", "afs_shape", "dominant_proportion", "maf",
             "mean_depth", "heterozygosity", "ploidy"]
    rows = []
    for k in order:
        if k not in m and k not in s:
            continue
        val = m.get(k, s.get(k))
        disp = _shape_cn(val) if k == "afs_shape" else _fmt(k, val)
        rows.append(f"<tr><td class='k'>{k}</td><td class='v'>{disp}</td>"
                    f"<td class='e'>{_METRIC_EXPLAIN.get(k,'')}</td></tr>")
    imgs = []
    for label_cn, key in [("等位频率谱(VAF/AFS)直方图", "afs"),
                          ("GenomeScope k-mer 基因组评估", "genomescope"),
                          ("Smudgeplot 倍性/杂合云图", "smudgeplot")]:
        p = s.get("images", {}).get(key)
        b64 = _embed_image(p) if p else None
        imgs.append(
            f"<figure><figcaption>{label_cn}</figcaption>"
            + (f"<img src='data:image/png;base64,{b64}'>" if b64
               else "<div class='noimg'>（无图:相应步骤未产出或被跳过）</div>")
            + "</figure>")
    return (
        f"<details class='sample' id='{sample}'>{sumrow}"
        f"<div class='sample-body'>"
        f"<p class='rationale'>{rationale}</p>"
        f"<table class='metrics'><tr><th>指标</th><th>数值</th><th>通俗解释</th></tr>{''.join(rows)}</table>"
        f"<div class='figs'>{''.join(imgs)}</div>"
        f"</div></details>")



def _CSS() -> str:
    return """
    body{font-family:'PingFang SC','Microsoft YaHei',sans-serif;max-width:1100px;margin:24px auto;
         padding:0 16px;color:#222;line-height:1.6}
    h1{color:#1a237e;border-bottom:3px solid #1a237e;padding-bottom:8px}
    h2{color:#283593;margin-top:36px;border-left:5px solid #283593;padding-left:10px}
    .muted{color:#666;font-size:0.92em}
    table{border-collapse:collapse;width:100%;margin:12px 0;font-size:0.95em}
    th,td{border:1px solid #ccc;padding:7px 10px;text-align:left;vertical-align:top}
    th{background:#e8eaf6}
    table.sum td:nth-child(2){font-weight:bold}
    .legend dt{font-weight:bold;color:#1a237e;margin-top:10px}
    .legend dd{margin:4px 0 0 0;color:#444}
    details.sample{margin:14px 0;border:1px solid #ddd;border-radius:6px;overflow:hidden}
    summary.sumrow{display:flex;align-items:center;gap:14px;flex-wrap:wrap;padding:12px 16px;
                   cursor:pointer;font-size:1.02em;background:#fafafa;list-style:none}
    summary.sumrow::-webkit-details-marker{display:none}
    details.sample > summary.sumrow::before{content:'▶ ';color:#888;font-size:0.85em}
    details.sample[open] > summary.sumrow::before{content:'▼ '}
    .sumrow .sname{font-weight:bold;font-size:1.15em;min-width:70px}
    .sumrow .sverdict{font-weight:bold}
    .sumrow .smeta{color:#555;font-size:0.92em}
    .sumrow .conf{margin-left:auto;color:#777;font-size:0.88em}
    .sample-body{padding:4px 16px 16px}
    .toggle-bar{margin:12px 0}
    .toggle-bar button{padding:6px 14px;margin-right:8px;border:1px solid #1a237e;background:#1a237e;
                       color:#fff;border-radius:4px;cursor:pointer;font-size:0.9em}
    .toggle-bar button:hover{background:#283593}
    .rationale{background:#fafafa;padding:10px 14px;border-radius:4px;margin:10px 0}
    .metrics td.k{font-weight:bold;color:#283593;width:16%}
    .metrics td.v{width:16%;font-family:monospace}
    .metrics td.e{color:#555}
    .figs{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:14px}
    figure{margin:0;text-align:center;background:#fafafa;padding:8px;border-radius:6px;border:1px solid #eee;min-width:0}
    figcaption{font-size:0.85em;color:#555;margin-bottom:6px}
    figure img{max-width:100%;width:100%;height:auto;border:1px solid #ddd;display:block}
    .noimg{color:#aaa;font-size:0.85em;padding:40px 20px}
    .tree{margin:14px 0;text-align:center}
    .tree img{max-width:100%;height:auto;border:1px solid #ddd;border-radius:6px}
    """
