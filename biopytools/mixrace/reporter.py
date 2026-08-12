"""mixrace 报告生成|mixrace report generation.

step08: VAF 直方图 R 脚本生成、per-sample md 报告、跨样品对比表(tsv+html)、
genomescope summary 解析。|VAF histogram R script, per-sample markdown report,
cross-sample summary table (tsv+html), genomescope summary parsing.
"""
from pathlib import Path
from typing import Tuple

_SUMMARY_COLS = ["sample", "verdict", "confidence", "vaf_mid_ratio",
                 "multiallelic_ratio", "fws", "mean_depth"]


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
        metrics: 指标 dict(vaf_mid_ratio/multiallelic_ratio/fws/hw/mean_depth/heterozygosity?)
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
        md.append(f"| {k} | {v} |")
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
        html += "<tr>" + "".join(f"<td>{r.get(c, '')}</td>" for c in _SUMMARY_COLS) + "</tr>"
    html += "</table>"
    return tsv, html
