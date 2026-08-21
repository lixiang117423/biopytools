"""mixrace 全套图(matplotlib)|mixrace figure suite.

从 03_het_eval/ 的 TSV 表产图;单图失败仅告警不阻断(graceful degradation)。
中文字体优先注册 Noto Sans CJK/WenQuanYi,缺失则标签退英文(禁豆腐块)。
|Figures from het_eval tables; per-figure failures degrade gracefully.
"""
from pathlib import Path
from typing import List, Optional

import numpy as np

matplotlib = None
_CJK_OK = None
_VERDICT_COLOR = {"pure": "#2e7d32", "divergent": "#ef6c00",
                  "contaminated": "#c62828", "uncertain": "#9e9e9e"}


def _t(cn: str, en: str) -> str:
    """中文标签(无中文字体时退英文)|CJK label with English fallback."""
    return cn if _CJK_OK else en


def _setup_font() -> bool:
    """注册中文字体(幂等)|register CJK font (idempotent)."""
    global matplotlib, _CJK_OK
    if _CJK_OK is not None:
        return _CJK_OK
    import matplotlib as mpl
    mpl.use("Agg")
    matplotlib = mpl
    import matplotlib.pyplot  # noqa: F401  绑定 matplotlib.pyplot 属性|bind pyplot attr
    from matplotlib import font_manager
    ok = False
    # 字体名含空格(如 "Noto Sans CJK SC"),匹配时双方都去空格;补 macOS PingFang
    # |font names contain spaces ("Noto Sans CJK SC"); strip spaces on both sides
    pats = ("NotoSansCJK", "NotoSansCJKsc", "wqy", "WenQuanYi", "SourceHanSans",
            "PingFangSC", "HeitiSC", "SimHei")
    try:
        for pat in pats:
            hits = [f for f in font_manager.fontManager.ttflist
                    if pat.lower() in f.name.replace(" ", "").lower()]
            if hits:
                matplotlib.rcParams["font.family"] = hits[0].name
                ok = True
                break
    except Exception:
        pass
    matplotlib.rcParams["axes.unicode_minus"] = False
    _CJK_OK = ok
    return ok


def _read_tsv(path) -> List[dict]:
    p = Path(path)
    if not p.exists():
        return []
    lines = [l.split("\t") for l in p.read_text(encoding="utf-8").splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    return [dict(zip(lines[0], f)) for f in lines[1:]]


def _save(fig, out) -> Path:
    fig.savefig(out, dpi=150, bbox_inches="tight")
    matplotlib.pyplot.close(fig)
    return Path(out)


def _verdict_of(sample, rows):
    for r in rows:
        if r.get("sample") == sample:
            return str(r.get("verdict", "uncertain"))
    return "uncertain"


def _window_matrix(eval_dir):
    """L3 长表 → 样本×窗口矩阵|samples x windows matrix from L3 table."""
    rows = _read_tsv(eval_dir / "l3_window_het.tsv")
    if not rows:
        return None
    wins = sorted({(r["chrom"], int(r["start"])) for r in rows},
                  key=lambda k: (k[0], k[1]))
    samples = sorted({r["sample"] for r in rows})
    w2i = {k: i for i, k in enumerate(wins)}
    s2i = {s: i for i, s in enumerate(samples)}
    mat = np.full((len(samples), len(wins)), np.nan)
    for r in rows:
        mat[s2i[r["sample"]], w2i[(r["chrom"], int(r["start"]))]] = float(r["het_rate"])
    return samples, [f"{c}:{s//1000}kb" for c, s in wins], mat


def fig_heatmap(eval_dir, rows, out, exclude_mask=None):
    """杂合热图(排除热点对比=双面板)|het heatmap (optionally before/after panels)."""
    _setup_font()
    plt = matplotlib.pyplot
    wm = _window_matrix(eval_dir)
    if wm is None:
        raise FileNotFoundError("l3_window_het.tsv")
    samples, cols, mat = wm
    if exclude_mask is None:
        fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.18), max(4, len(samples) * 0.35)))
        ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
        ax.set_yticks(range(len(samples)))
        ax.set_yticklabels(samples, fontsize=8)
        ax.set_title(_t("100kb 窗口杂合率热图", "Het-rate heatmap (100kb windows)"), fontsize=11)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(16, max(4, len(samples) * 0.35)),
                                 sharey=True)
        for ax, m, ttl in zip(axes, (mat, np.where(exclude_mask[None, :], np.nan, mat)),
                              (_t("排除前", "before"), _t("排除后", "after"))):
            ax.imshow(m, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
            ax.set_title(ttl, fontsize=10)
        axes[0].set_yticks(range(len(samples)))
        axes[0].set_yticklabels(samples, fontsize=8)
    return _save(fig, out)


def fig_manhattan_grid(eval_dir, rows, out, only_non_pure=False):
    """窗口级 Manhattan 总览(每样本一行散点)|window-level Manhattan overview."""
    _setup_font()
    plt = matplotlib.pyplot
    wm = _window_matrix(eval_dir)
    if wm is None:
        raise FileNotFoundError("l3_window_het.tsv")
    samples, cols, mat = wm
    if only_non_pure:
        keep = [i for i, s in enumerate(samples)
                if _verdict_of(s, rows) != "pure"]
        if not keep:
            raise ValueError("no non-pure samples")
        samples, mat = [samples[i] for i in keep], mat[keep]
    x = np.arange(mat.shape[1])
    fig, axes = plt.subplots(len(samples), 1, figsize=(12, 1.1 * len(samples)),
                             sharex=True, squeeze=False)
    for ax, s, rates in zip(axes[:, 0], samples, mat):
        ax.scatter(x, rates, s=2, color="#4C78A8")
        ax.set_ylabel(s, rotation=0, ha="right", va="center", fontsize=8)
        ax.set_ylim(0, 1)
    axes[-1, 0].set_xlabel(_t("基因组窗口(按染色体排序)", "genomic windows (sorted)"), fontsize=9)
    fig.suptitle(_t("全基因组杂合分布", "Genome-wide het distribution"), fontsize=11)
    return _save(fig, out)


def fig_dist_heatmap(eval_dir, out):
    """距离热图(格中数值)|distance heatmap with cell values."""
    _setup_font()
    plt = matplotlib.pyplot
    rows = _read_tsv(eval_dir / "distance_matrix.tsv")
    if not rows:
        raise FileNotFoundError("distance_matrix.tsv")
    samples = [r["sample"] for r in rows]
    try:
        # 无共享位点的配对为 NaN/空串,统一转 NaN(imshow 可显示,不整图失败)
        # |pairs with no shared sites are NaN; empty string -> NaN, not a crash
        mat = np.array([[float(r[s]) if str(r.get(s, "")).strip() else float("nan")
                         for s in samples] for r in rows])
    except (KeyError, ValueError):
        raise ValueError("距离矩阵格式坏|bad distance matrix")
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=0.7)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=8)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=8)
    for i in range(len(samples)):
        for j in range(len(samples)):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6)
    ax.set_title(_t("样本间 SNP 不匹配率", "Pairwise SNP mismatch rate"), fontsize=11)
    return _save(fig, out)


def fig_pca(eval_dir, rows, out):
    """PCA 三视图(PC1-2/1-3/2-3,判读着色)|PCA 3 views colored by verdict."""
    _setup_font()
    plt = matplotlib.pyplot
    table = _read_tsv(eval_dir / "pca_coords.tsv")
    if not table:
        raise FileNotFoundError("pca_coords.tsv")
    samples = [r["sample"] for r in table]
    pc = {k: np.array([float(r.get(k, 0) or 0) for r in table])
          for k in ("PC1", "PC2", "PC3")}
    colors = [_VERDICT_COLOR.get(_verdict_of(s, rows), "#9e9e9e") for s in samples]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (a, b) in zip(axes, (("PC1", "PC2"), ("PC1", "PC3"), ("PC2", "PC3"))):
        ax.scatter(pc[a], pc[b], c=colors, s=40)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
        for s, x, y in zip(samples, pc[a], pc[b]):
            ax.annotate(s, (x, y), fontsize=7)
    fig.suptitle(_t("PCA(绿=纯 橙=差异 红=混杂)", "PCA (green=pure orange=divergent red=contam)"),
                 fontsize=11)
    return _save(fig, out)


def fig_tree(eval_dir, out):
    """NJ 树图(权值>阈值才标)|NJ tree via Bio.Phylo."""
    _setup_font()
    plt = matplotlib.pyplot
    nwk = eval_dir / "nj_tree.nwk"
    if not nwk.exists():
        raise FileNotFoundError("nj_tree.nwk")
    from Bio import Phylo
    tree = Phylo.read(str(nwk), "newick")
    fig, ax = plt.subplots(figsize=(9, 7))
    Phylo.draw(tree, axes=ax, do_show=False,
               show_confidence=False, label_func=lambda c: c.name or "")
    ax.set_title(_t("NJ 系统发育树(SNP不匹配率)", "NJ tree (SNP mismatch)"), fontsize=11)
    return _save(fig, out)


def fig_altfrac(eval_dir, rows, out):
    """非纯样本杂合位点 altfrac 叠加分布|altfrac histograms (non-pure samples)."""
    _setup_font()
    plt = matplotlib.pyplot
    from .het_eval import load_gt_ad_dp
    tsv = eval_dir / "gt_ad_dp.tsv"
    if not tsv.exists():
        raise FileNotFoundError("gt_ad_dp.tsv")
    d = load_gt_ad_dp(str(tsv))
    non_pure = [s for s in d.samples if _verdict_of(s, rows) != "pure"]
    if not non_pure:
        raise ValueError("no non-pure samples")
    fig, ax = plt.subplots(figsize=(8, 5))
    idx = [d.samples.index(s) for s in non_pure]
    for i, s in zip(idx, non_pure):
        het = d.gt[i] == 1
        tot = d.ref_ad[i][het] + d.alt_ad[i][het]
        af = np.where(tot > 0, d.alt_ad[i][het] / np.maximum(tot, 1), 0.0)
        ax.hist(af, bins=40, range=(0, 1), histtype="step", linewidth=1.5,
                label=s, density=True)
    ax.set_xlabel("altfrac")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    ax.set_title(_t("杂合位点 altfrac 分布(非纯样本)", "altfrac of het sites (non-pure)"),
                 fontsize=11)
    return _save(fig, out)


def fig_three_panel(eval_dir, out):
    """三面板:总杂合率/DP50杂合率/稳健杂合率|three-panel: het / DP50 het / robust."""
    _setup_font()
    plt = matplotlib.pyplot
    rows = _read_tsv(eval_dir / "l1_het_stats.tsv")
    if not rows:
        raise FileNotFoundError("l1_het_stats.tsv")
    rows = sorted(rows, key=lambda r: float(r.get("het_rate") or 0))
    samples = [r["sample"] for r in rows]

    def col(key):
        return np.array([float(r.get(key) or 0) for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15, max(4, 0.3 * len(samples))), sharey=True)
    y = np.arange(len(samples))
    for ax, key, ttl in zip(axes, ("het_rate", "het_rate_dp50", "robust_rate"),
                            (_t("总杂合率", "het rate"), _t("DP>=50 杂合率", "het@DP50"),
                             _t("稳健杂合率", "robust rate"))):
        ax.barh(y, col(key), color="#4C78A8")
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("rate")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(samples, fontsize=8)
    return _save(fig, out)


def _hotspot_mask(eval_dir, window: int = 100000) -> Optional[np.ndarray]:
    """热点 BED → 窗口布尔掩码(按 L3 窗口顺序)|hotspot BED -> window mask."""
    wm = _window_matrix(eval_dir)
    bed = eval_dir / "hotspots.bed"
    if wm is None or not bed.exists():
        return None
    samples, cols, mat = wm
    mask = np.zeros(mat.shape[1], dtype=bool)
    for line in bed.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        key = f"{f[0]}:{int(f[1]) // 1000}kb"
        for j, c in enumerate(cols):
            if c == key:
                mask[j] = True
    return mask


def build_figures(config, logger, out_dir, payloads: dict) -> List[Path]:
    """产出全部图;单图失败告警继续|produce all figures; per-figure degrade."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = payloads.get("rows", [])
    eval_dir = payloads.get("het_eval_dir")
    hotspot_mask = _hotspot_mask(eval_dir)
    plan = [
        ("het_heatmap_100kb.png", lambda: fig_heatmap(eval_dir, rows, out / "het_heatmap_100kb.png")),
        ("het_heatmap_excl_hotspots.png", lambda: fig_heatmap(
            eval_dir, rows, out / "het_heatmap_excl_hotspots.png", exclude_mask=hotspot_mask)),
        ("het_genome_overview.png", lambda: fig_manhattan_grid(
            eval_dir, rows, out / "het_genome_overview.png")),
        ("manhattan_grid.png", lambda: fig_manhattan_grid(
            eval_dir, rows, out / "manhattan_grid.png", only_non_pure=True)),
        ("dist_heatmap.png", lambda: fig_dist_heatmap(eval_dir, out / "dist_heatmap.png")),
        ("pca_3view.png", lambda: fig_pca(eval_dir, rows, out / "pca_3view.png")),
        ("nj_tree.png", lambda: fig_tree(eval_dir, out / "nj_tree.png")),
        ("altfrac_dist.png", lambda: fig_altfrac(eval_dir, rows, out / "altfrac_dist.png")),
        ("eval_3panel.png", lambda: fig_three_panel(eval_dir, out / "eval_3panel.png")),
    ]
    done: List[Path] = []
    for name, fn in plan:
        try:
            done.append(fn())
        except Exception as e:
            if logger is not None:
                logger.warning(f"图未生成|figure skipped [{name}]: {e}")
    return done
