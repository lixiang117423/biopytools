"""mixrace 杂合评估引擎|mixrace heterozygosity evaluation engine.

移植导师四层方法学(纯计算层): L1 AD/DP 排错 → L2 shared/private+混合伴侣 →
L3 窗口分布 → L4 热点排除 + 群体结构(距离/PCA/NJ)。
输入为 bcftools query 长表(GT/AD/DP)。
|Mentor's four-layer methodology (pure computation). Input = bcftools query
long table (GT/AD/DP).
"""
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from .utils import build_conda_command

# bcftools query 长表格式: 位点4列 + 每样本 GT/AD/DP 三列
# |long-table format: 4 site columns + per-sample GT/AD/DP
QUERY_FORMAT = r"%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n"


@dataclass
class HetData:
    """双等位 SNP 矩阵|biallelic SNP matrices.

    gt: int8 N×M,0/0→0 0/1→1 1/1→2 缺失→-1;ref_ad/alt_ad/dp: int16 N×M
    """

    samples: List[str]
    gt: np.ndarray
    ref_ad: np.ndarray
    alt_ad: np.ndarray
    dp: np.ndarray
    chrom: List[str]
    pos: np.ndarray


def _enc_gt(gt: str) -> int:
    """GT 文本 → 编码(相位分隔符归一)|encode GT text (phase-insensitive)."""
    a1, a2 = gt.replace("|", "/").split("/")
    if a1 == "." or a2 == ".":
        return -1
    if a1 == a2:
        return 0 if a1 == "0" else 2
    return 1


def load_gt_ad_dp(tsv_path: str, samples: Optional[List[str]] = None) -> HetData:
    """读 bcftools query 长表,仅双等位 SNP|load long table, biallelic SNPs only.

    兼容两种格式:带 `{sample}_GT` 表头行(test fixture/外部生成)或 bcftools 原生
    无表头输出(此时样本名取 samples 参数,缺省 S1..SN)。
    |Accepts both a `{sample}_GT` header line and bcftools' native headerless
    output (names from `samples`, else S1..SN).
    """
    names: List[str] = []
    cols: list = []
    chrom: List[str] = []
    pos: list = []
    first = True
    with open(tsv_path, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if first:
                first = False
                if len(p) > 4 and p[4].endswith("_GT") and p[-1].endswith("_DP"):
                    names = [c[:-3] for c in p[4:] if c.endswith("_GT")]
                    continue   # 表头行|header line
                # 无表头:按列数推断样本数|headerless: infer N from column count
                n = (len(p) - 4) // 3
                names = list(samples[:n]) if samples else []
                names += [f"S{i+1}" for i in range(len(names), n)]
            if len(p[2]) != 1 or len(p[3]) != 1 or "," in p[3]:
                continue   # 非双等位SNP(多等位/indel)|not biallelic SNP
            chrom.append(p[0])
            pos.append(int(p[1]))
            row = []
            for i in range(len(names)):
                try:
                    gt, ad, dp = p[4 + i*3], p[5 + i*3], p[6 + i*3]
                except IndexError:
                    break   # 列数不足的截断行:丢弃|truncated row: drop
                adp = ad.split(",")
                ref = int(adp[0]) if adp[0] not in (".", "") else 0
                alt = int(adp[1]) if len(adp) > 1 and adp[1] not in (".", "") else 0
                try:
                    d = int(dp)
                except ValueError:
                    d = ref + alt   # DP 坏值回退 ref+alt|fallback
                g = _enc_gt(gt) if gt not in (".", "./.", ".|.") else -1
                row.append((g, ref, alt, d))
            if len(row) != len(names):
                continue   # 截断行|truncated row
            cols.append(row)
    samples = names
    arr = (np.array(cols, dtype=np.int32) if cols
           else np.zeros((0, len(samples), 4), np.int32))
    # (位点,样本,4) → 转置为 接口约定 (样本,位点)|(sites,samples,4) -> (samples,sites)
    arr = arr.transpose(1, 0, 2)
    return HetData(samples=samples, gt=arr[:, :, 0].astype(np.int8),
                   ref_ad=arr[:, :, 1].astype(np.int32), alt_ad=arr[:, :, 2].astype(np.int32),
                   dp=arr[:, :, 3].astype(np.int32), chrom=chrom,
                   pos=np.array(pos, dtype=np.int32))


def build_query_cmd(bcftools_path: str, vcf: str, out_tsv: str) -> List[str]:
    """构造长表提取命令(conda 包装)|build the query command (conda-wrapped)."""
    return build_conda_command(
        bcftools_path, ["query", "-f", QUERY_FORMAT, "-o", out_tsv, vcf])


def list_samples(runner, bcftools_path: str, vcf: str) -> List[str]:
    """VCF 样本列表(bcftools query -l)|sample list from VCF header."""
    ok, out, _ = runner.run_conda(bcftools_path, ["query", "-l", vcf],
                                  "列样本|list samples")
    if not ok:
        return []
    return [x for x in out.splitlines() if x.strip()]


def _safe_dp_ratio(dp_row, het_mask, hom_dp) -> float:
    """DP检验:杂合DP中位/纯合DP中位(除零/全空 → nan)|het/hom DP median ratio (nan-safe)."""
    het_med = float(np.median(dp_row[het_mask])) if het_mask.any() else float("nan")
    hom_med = float(np.median(hom_dp)) if hom_dp.size else float("nan")
    if np.isfinite(het_med) and hom_med and np.isfinite(hom_med):
        return het_med / hom_med
    return float("nan")


def compute_l1(d: HetData, dp_layers=(10, 20, 50)) -> List[dict]:
    """L1 逐样本:杂合率+稳健杂合+DP检验(导师 het_mixture_analysis/het_final 口径)。

    robust = altAD>=5 且 altfrac>=0.2;dp_ratio = 杂合DP中位/纯合DP中位(>1.5 提示
    真实混合:杂合位点深度反而更高)。|robust het + DP ratio (real-mixture signal).
    """
    rows = []
    for i, s in enumerate(d.samples):
        called = d.gt[i] >= 0
        het = d.gt[i] == 1
        n_sites = int(called.sum())
        n_het = int(het.sum())
        tot = d.ref_ad[i] + d.alt_ad[i]
        af = np.where(tot > 0, d.alt_ad[i] / np.maximum(tot, 1), 0.0)
        robust = int(((d.alt_ad[i] >= 5) & (af >= 0.2) & het).sum())
        hom_dp = d.dp[i][called & ~het]
        row = {"sample": s, "n_sites": n_sites, "n_het": n_het,
               "het_rate": n_het / n_sites if n_sites else 0.0,
               "robust_het": robust,
               "robust_rate": robust / n_sites if n_sites else 0.0,
               "median_altfrac": float(np.median(af[het])) if n_het else float("nan"),
               "median_altAD": float(np.median(d.alt_ad[i][het])) if n_het else float("nan"),
               "median_hetDP": float(np.median(d.dp[i][het])) if n_het else float("nan"),
               "dp_ratio": _safe_dp_ratio(d.dp[i], het, hom_dp)}
        for layer in dp_layers:
            m = called & (d.dp[i] >= layer)
            n_m = int(m.sum())
            row[f"het_rate_dp{layer}"] = int((het & m).sum()) / n_m if n_m else 0.0
        rows.append(row)
    return rows


def compute_shared_partner(d: HetData) -> dict:
    """L2 shared/private ALT + 混合伴侣矩阵(导师 het_private.py 口径)。

    partner_alt_rate[A][B] = A 的杂合位点上 B 携带ALT 计数 / |A的全部杂合位点|
    (分母=全部杂合数,导师 Pb9=154360/156863);partner_hom_rate = 其中 B 纯合 1/1
    占 B 携带的比例(区分 Pb9-Pb22 型 vs 群内互为杂合型)。B 缺失=不携带。
    |shared/private ALT + partner matrices; missing B counts as non-carrier.
    """
    n = len(d.samples)
    par_alt = np.zeros((n, n))
    par_hom = np.zeros((n, n))
    rows = []
    for i, s in enumerate(d.samples):
        het_mask = d.gt[i] == 1
        n_het = int(het_mask.sum())
        n_sites = int((d.gt[i] >= 0).sum())
        if n_het == 0:
            rows.append({"sample": s, "n_het": 0, "shared": 0, "private": 0,
                         "pct_private": 0.0, "shared_only_rate": 0.0})
            continue
        others = np.delete(np.arange(n), i)
        g = d.gt[others][:, het_mask]
        carry_any = ((g == 1) | (g == 2)).any(axis=0)
        shared = int(carry_any.sum())
        private = n_het - shared
        tot = d.ref_ad[i][het_mask] + d.alt_ad[i][het_mask]
        af = np.where(tot > 0, d.alt_ad[i][het_mask] / np.maximum(tot, 1), 0.0)
        rows.append({"sample": s, "n_het": n_het, "shared": shared, "private": private,
                     "pct_private": private / n_het,
                     "shared_only_rate": shared / n_sites if n_sites else 0.0,
                     "median_altfrac_shared": (float(np.median(af[carry_any]))
                                               if shared else float("nan")),
                     "median_altfrac_private": (float(np.median(af[~carry_any]))
                                                if private else float("nan")),
                     "median_altAD_shared": (float(np.median(d.alt_ad[i][het_mask][carry_any]))
                                             if shared else float("nan")),
                     "median_altAD_private": (float(np.median(d.alt_ad[i][het_mask][~carry_any]))
                                              if private else float("nan"))})
        for k, j in enumerate(others):
            n_carry = int((g[k] >= 1).sum())
            n_hom = int((g[k] == 2).sum())
            par_alt[i, j] = n_carry / n_het
            par_hom[i, j] = n_hom / n_carry if n_carry else 0.0
    return {"rows": rows, "partner_alt_rate": par_alt, "partner_hom_rate": par_hom}


def _site_window_index(d: HetData, window: int):
    """位点 → 窗口索引(键=(chrom, pos//window*window),按 chrom+start 排序)。
    |site -> window index (key sorted by chrom+start)."""
    keys = sorted({(c, int(p) // window * window) for c, p in zip(d.chrom, d.pos)})
    k2i = {k: i for i, k in enumerate(keys)}
    site_idx = np.array([k2i[(c, int(p) // window * window)]
                         for c, p in zip(d.chrom, d.pos)], dtype=np.int64)
    chrom_list = [k[0] for k in keys]
    start = np.array([k[1] for k in keys], dtype=np.int64)
    return chrom_list, start, site_idx


def compute_windows(d: HetData, window: int = 100_000) -> dict:
    """L3 100kb 窗口杂合率矩阵|windowed het-rate matrix (N×W)."""
    chrom_list, start, site_idx = _site_window_index(d, window)
    w = len(chrom_list)
    n = len(d.samples)
    het = np.zeros((n, w))
    called = np.zeros((n, w))
    for i in range(n):
        np.add.at(het[i], site_idx[d.gt[i] == 1], 1.0)
        np.add.at(called[i], site_idx[d.gt[i] >= 0], 1.0)
    rate = np.where(called > 0, het / np.where(called > 0, called, 1), 0.0)
    return {"chrom": chrom_list, "start": start, "het": rate,
            "n_sites": called.astype(int), "site_window": site_idx}


def compute_hotspots(windows: dict, l1_rows: List[dict], d: HetData, params: dict) -> dict:
    """L4 共享热点识别+排除(导师 ≥半数候选 & >2×自身均值 & 中位>10% 规则自动化)。

    候选=het_rate>=pure_threshold 的样本;K=max(3, ceil(候选/2));
    extra_bed(chrom,start,end) 与窗口相交则并入排除集(不进 hot 清单)。
    |auto hotspots + extra-bed merge; excluded rates recomputed at site level.
    """
    import math
    pure = params.get("pure_threshold", 0.001)
    fold = params.get("fold", 2.0)
    min_median = params.get("min_median", 0.10)
    het = windows["het"]
    w = het.shape[1]
    cand = [i for i, r in enumerate(l1_rows) if r["het_rate"] >= pure]
    auto_mask = np.zeros(w, dtype=bool)
    hot = []
    if cand:
        k_min = max(3, math.ceil(len(cand) / 2))
        genome_rate = np.array([l1_rows[i]["het_rate"] for i in cand])
        exceed = het[cand] > fold * genome_rate[:, None]
        n_hit = exceed.sum(axis=0).astype(int)
        med = np.median(het[cand], axis=0)
        auto_mask = (n_hit >= k_min) & (med > min_median)
        for j in np.where(auto_mask)[0]:
            hot.append({"chrom": windows["chrom"][j],
                        "start": int(windows["start"][j]),
                        "n_snps": int(windows["n_sites"][:, j].sum()),
                        "median_het": float(med[j]), "n_hit": int(n_hit[j])})
    mask = auto_mask.copy()
    win = params.get("window_size")
    if not win:
        # 回退:窗口宽度=相邻start差(优先用参数,与 hotspots.bed 输出口径一致)
        # |fallback: window width from adjacent starts (prefer the param)
        starts = windows["start"]
        win = int(starts[1] - starts[0]) if len(starts) > 1 else 1
    for c, s, e in params.get("extra_bed", []):
        for j, (cj, sj) in enumerate(zip(windows["chrom"], windows["start"])):
            if cj == c and sj < e and sj + win > s:
                mask[j] = True
    # 位点级重算排除后杂合率/稳健杂合率|site-level recompute after exclusion
    keep = ~mask[windows["site_window"]]
    excluded = []
    for i, s in enumerate(d.samples):
        called = (d.gt[i] >= 0) & keep
        hetm = (d.gt[i] == 1) & keep
        tot = d.ref_ad[i] + d.alt_ad[i]
        af = np.where(tot > 0, d.alt_ad[i] / np.maximum(tot, 1), 0.0)
        robust_all = ((d.alt_ad[i] >= 5) & (af >= 0.2) & called)
        n_c = int(called.sum())
        excluded.append({
            "sample": s,
            "het_rate_before": l1_rows[i]["het_rate"],
            "robust_rate_before": l1_rows[i]["robust_rate"],
            "n_sites_after": n_c,
            "het_rate_after": int(hetm.sum()) / n_c if n_c else 0.0,
            "robust_rate_after": int((robust_all & hetm).sum()) / n_c if n_c else 0.0})
    return {"hot": hot, "excluded": excluded, "mask": mask}


def distance_matrix(d: HetData) -> np.ndarray:
    """成对 SNP 不匹配率(双方均 call 位点)|pairwise mismatch rate (both-called sites).

    导师 het_distance.py 口径:GT 编码不等即不匹配;n_shared=0 → nan。
    """
    n = len(d.samples)
    dist = np.full((n, n), np.nan)
    np.fill_diagonal(dist, 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            m = (d.gt[i] >= 0) & (d.gt[j] >= 0)
            ns = int(m.sum())
            if ns:
                v = float((d.gt[i, m] != d.gt[j, m]).sum()) / ns
                dist[i, j] = dist[j, i] = v
    return dist


def pca_coords(d: HetData, n_comp: int = 4) -> Tuple[np.ndarray, np.ndarray]:
    """GT 编码(缺失填0)中心化后 SVD|PCA via SVD on centered GT codes."""
    x = d.gt.astype(float)
    x[x < 0] = 0.0
    x = x - x.mean(axis=1, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    k = min(n_comp, len(s))
    var = s[:k] ** 2
    explained = var / var.sum() if var.sum() else np.zeros(k)
    return u[:, :k] * s[:k], explained


def nj_newick(dist: np.ndarray, samples: List[str]) -> str:
    """NJ 树 newick(biopython;nan→成对均值填充)|NJ tree via biopython."""
    import io
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio import Phylo
    n = len(samples)
    if n < 3:
        # 样品不足 3 个无法建树:返回空(调用方跳过写文件,fig_tree 自动降级)
        # |fewer than 3 samples: no tree; caller skips the file
        return ""
    dm = np.array(dist, dtype=float)
    finite = dm[np.isfinite(dm)]
    fill = float(finite.mean()) if finite.size else 0.0
    dm[~np.isfinite(dm)] = fill
    np.fill_diagonal(dm, 0.0)
    lower = [[float(dm[i, j]) for j in range(i + 1)] for i in range(n)]
    matrix = DistanceMatrix(names=list(samples), matrix=lower)
    tree = DistanceTreeConstructor().nj(matrix)
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue()


def cluster_groups(dist: np.ndarray, threshold: float = 0.05) -> List[int]:
    """单链接分群(展示用,不参与判读)|single-linkage groups (display only)."""
    n = dist.shape[0]
    dm = np.where(np.isfinite(dist), dist, np.inf)
    labels = list(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if dm[i, j] < threshold:
                old, new = max(labels[i], labels[j]), min(labels[i], labels[j])
                labels = [new if x == old else x for x in labels]
    uniq = {v: i for i, v in enumerate(sorted(set(labels)))}
    return [uniq[x] for x in labels]


# ---------- 编排层|orchestration ----------

def write_tsv(path, rows: List[dict]) -> None:
    """写 TSV(表头=首行键序;nan→空)|write TSV (header from first row; nan->empty)."""
    from pathlib import Path as _P
    if not rows:
        _P(path).write_text("", encoding="utf-8")
        return
    head = list(rows[0].keys())
    lines = ["\t".join(head)]
    for r in rows:
        vals = []
        for k in head:
            v = r.get(k, "")
            if isinstance(v, float) and v != v:
                v = ""
            vals.append(str(v))
        lines.append("\t".join(vals))
    _P(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_bed_regions(bed: str) -> List[Tuple[str, int, int]]:
    """读 BED(chrom,start,end;跳注释)|read BED regions."""
    out = []
    try:
        text = open(bed, encoding="utf-8").read()
    except OSError:
        return out
    for line in text.splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) >= 3:
            try:
                out.append((f[0], int(f[1]), int(f[2])))
            except ValueError:
                continue
    return out


def _top_partner(i: int, samples: List[str], par_alt, par_hom):
    """A 的最佳伴侣(ALT携带率最高,排除自身)|best partner for sample i."""
    alt_row = np.asarray(par_alt[i], dtype=float).copy()
    if alt_row.size <= 1:
        return {"top_partner": "", "top_partner_alt_rate": 0.0,
                "top_partner_hom_rate": 0.0}
    alt_row[i] = -1.0
    j = int(np.argmax(alt_row))
    return {"top_partner": samples[j],
            "top_partner_alt_rate": float(par_alt[i, j]),
            "top_partner_hom_rate": float(par_hom[i, j])}


def read_verdict_table(config) -> list:
    """读已判读表(step4/5 重跑与断点跳过共用)|read verdict table (reruns + checkpoint)."""
    from pathlib import Path as _P
    path = _P(config.output_dir) / "03_het_eval" / "verdict_table.tsv"
    if not path.exists():
        return []
    lines = [l.split("\t") for l in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 2:
        return []
    head = lines[0]
    rows = []
    for f in lines[1:]:
        r = dict(zip(head, f))
        for k in ("het_rate", "robust_rate", "shared_only_rate", "median_altfrac",
                  "het_rate_after_hotspot", "dp_ratio", "top_partner_alt_rate",
                  "top_partner_hom_rate", "mix_proportion"):
            try:
                r[k] = float(r.get(k, "") or "nan")
            except ValueError:
                r[k] = None
        for k in ("n_sites", "n_het"):
            try:
                r[k] = int(r.get(k, 0))
            except (TypeError, ValueError):
                r[k] = 0
        rows.append(r)
    return rows


def run_het_eval(config, runner, ckpt, vcf: str) -> List[dict]:
    """step3: 长表提取→四层评估→群体结构→三分支判读|query + 4 layers + verdicts.

    Returns: 判读汇总行(空列表=失败)|verdict rows (empty on failure).
    """
    from pathlib import Path
    from .pipeline import _done
    from .verdict import judge
    out_dir = Path(config.output_dir) / "03_het_eval"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "gt_ad_dp.tsv"

    # 断点续传:整步已判读则直接读表返回|checkpoint: reuse verdict table
    if config.enable_checkpoint and _done(ckpt, "het_eval",
                                          out_dir / "verdict_table.tsv"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: het_eval")
        return read_verdict_table(config)

    # ① 长表提取|long-table query
    if config.enable_checkpoint and _done(ckpt, "query", tsv):
        runner.logger.info("跳过已完成步骤|Skipping completed step: query")
    else:
        runner.logger.info("开始步骤|Starting step: bcftools query 长表提取")
        cmd = build_query_cmd(config.bcftools_path, vcf, str(tsv))
        ok, _, _ = runner.run(cmd, "长表提取|query GT/AD/DP")
        if not (ok or tsv.exists()):
            runner.logger.error("长表提取失败,中止评估|query failed, eval aborted")
            return []
        if config.enable_checkpoint and ok:
            ckpt.create("query")
    if not tsv.exists():
        runner.logger.error(f"长表不存在: {tsv}|long table missing")
        return []

    # ② L1/L2+伴侣|layers 1-2 + partner matrices
    vcf_samples = list_samples(runner, config.bcftools_path, vcf)
    d = load_gt_ad_dp(str(tsv), vcf_samples or None)
    if len(d.pos) == 0 or len(d.samples) < 1:
        runner.logger.error("长表为空或样本数为0,中止评估|empty table or zero samples")
        return []
    from .utils import format_number
    runner.logger.info(f"双等位SNP位点|biallelic sites: {format_number(len(d.pos))} "
                       f"× 样本|samples {len(d.samples)}")
    l1 = compute_l1(d)
    sp = compute_shared_partner(d)
    write_tsv(out_dir / "l1_het_stats.tsv", l1)
    write_tsv(out_dir / "L2_shared_private.tsv", sp["rows"])
    shared_only = [{"sample": r["sample"], "n_het": r["n_het"],
                    "shared": r["shared"], "private": r["private"],
                    "shared_only_rate": r["shared_only_rate"],
                    "pct_private": r["pct_private"]} for r in sp["rows"]]
    write_tsv(out_dir / "l2_shared_only.tsv", shared_only)
    # 伴侣矩阵(全矩阵长表)|full partner matrices (long form)
    n = len(d.samples)
    par_rows = []
    for i, a in enumerate(d.samples):
        for j, b in enumerate(d.samples):
            if i == j:
                continue
            par_rows.append({"sample": a, "partner": b,
                             "partner_alt_rate": round(float(sp["partner_alt_rate"][i, j]), 4),
                             "partner_hom_rate": round(float(sp["partner_hom_rate"][i, j]), 4)})
    write_tsv(out_dir / "partner_matrix.tsv", par_rows)
    write_tsv(out_dir / "partner_top.tsv",
              [{**_top_partner(i, d.samples, sp["partner_alt_rate"], sp["partner_hom_rate"]),
                "sample": s} for i, s in enumerate(d.samples)])

    # ③ L3/L4 窗口+热点|windows + hotspots
    windows = compute_windows(d, config.window_size)
    win_rows = []
    for j in range(windows["het"].shape[1]):
        for i, s in enumerate(d.samples):
            win_rows.append({"chrom": windows["chrom"][j],
                             "start": int(windows["start"][j]),
                             "sample": s,
                             "het_rate": round(float(windows["het"][i, j]), 5),
                             "n_sites": int(windows["n_sites"][i, j])})
    write_tsv(out_dir / "l3_window_het.tsv", win_rows)
    hot_params = {"pure_threshold": config.pure_het_threshold,
                  "fold": config.hotspot_fold,
                  "min_median": config.hotspot_min_median,
                  "window_size": config.window_size,
                  "extra_bed": read_bed_regions(config.repeat_bed) if config.repeat_bed else []}
    hot = compute_hotspots(windows, l1, d, hot_params)
    write_tsv(out_dir / "l4_hotspot_windows.tsv", hot["hot"])
    with open(out_dir / "hotspots.bed", "w", encoding="utf-8") as fh:
        for j in np.where(hot["mask"])[0]:
            w = config.window_size
            fh.write(f"{windows['chrom'][j]}\t{int(windows['start'][j])}\t"
                     f"{int(windows['start'][j]) + w}\n")
    write_tsv(out_dir / "l4_hotspot_excluded_compare.tsv", hot["excluded"])
    excluded_by_sample = {r["sample"]: r for r in hot["excluded"]}

    # ④ 群体结构|population structure
    dist = distance_matrix(d)
    coords, explained = pca_coords(d)
    nwk = nj_newick(dist, d.samples)
    groups = cluster_groups(dist)
    if nwk:
        (out_dir / "nj_tree.nwk").write_text(nwk, encoding="utf-8")
    write_tsv(out_dir / "distance_matrix.tsv",
              [{"sample": d.samples[i],
                **{d.samples[j]: round(float(dist[i, j]), 4) for j in range(n)}}
               for i in range(n)])
    write_tsv(out_dir / "pca_coords.tsv",
              [{"sample": d.samples[i], "group": groups[i],
                **{f"PC{k+1}": round(float(coords[i, k]), 4) for k in range(coords.shape[1])},
                "explained": ";".join(f"{e:.4f}" for e in explained)}
               for i in range(n)])

    # ⑤ 三分支判读|three-branch verdicts
    jparams = {"pure_het": config.pure_het_threshold,
               "partner_alt_min": config.partner_alt_min,
               "partner_hom_min": config.partner_hom_min,
               "min_sites": config.min_sites}
    sp_by_sample = {row["sample"]: row for row in sp["rows"]}
    rows = []
    for r in l1:
        v = judge(r, sp["partner_alt_rate"], sp["partner_hom_rate"],
                  d.samples, jparams)
        ex = excluded_by_sample.get(r["sample"], {})
        rows.append({"sample": r["sample"], "verdict": v["verdict"],
                     "advice": v["advice"], "partner": v["partner"] or "",
                     "mix_proportion": v["mix_proportion"],
                     "subtag": v["subtag"],
                     "het_rate": r["het_rate"], "robust_rate": r["robust_rate"],
                     "n_sites": r["n_sites"], "n_het": r["n_het"],
                     "shared_only_rate": sp_by_sample[r["sample"]]["shared_only_rate"],
                     **_top_partner(d.samples.index(r["sample"]), d.samples,
                                    sp["partner_alt_rate"], sp["partner_hom_rate"]),
                     "median_altfrac": r["median_altfrac"],
                     "het_rate_after_hotspot": ex.get("het_rate_after", r["het_rate"]),
                     "dp_ratio": r["dp_ratio"],
                     "rationale": v["rationale"]})
    write_tsv(out_dir / "verdict_table.tsv", rows)
    if config.enable_checkpoint:
        ckpt.create("het_eval")
    for r in rows:
        runner.logger.info(
            f"{r['sample']}: 判读|verdict={r['verdict']} "
            f"杂合率|het={r['het_rate']*100:.4f}% "
            f"伴侣|partner={r['top_partner']}({r['top_partner_alt_rate']*100:.1f}%)")
    return rows
