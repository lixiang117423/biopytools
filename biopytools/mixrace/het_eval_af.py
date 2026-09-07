"""mixrace AF直接判定分支(论文式,不依赖GT)|paper-style GT-independent het branch.

方法学来自 Cao et al. 2026 (Front. Microbiol. 17:1789807):不使用 joint VCF 的
GT 字段(GATK 二倍体先验会把偏离 50% 的真实混合信号误判纯合,系统性低估
非对称混合样本的杂合度),改为直接扫每个位点 ref/alt reads 比例——
从联合 VCF 用 bcftools query 只取 AD/DP(等效逐 BAM mpileup 口径,免重跑 pileup)。

杂合度 = 杂合位点数/(杂合位点数+纯合变异位点数),分母不含纯合参考位点;
与 GT-based het_rate(杂合/全部called位点)是两种口径,并列输出供对比,不混用。
|Paper-style (Cao et al. 2026): AD/DP only from the joint VCF (equivalent to
per-BAM mpileup, without re-piling); het_rate_af = het/(het+hom_alt), a
denominator that EXCLUDES hom-ref sites — a different calibre than the
GT-based het_rate; both are reported side by side, never mixed.
"""
from typing import List, Optional

import numpy as np

from .het_eval import HetData, list_samples, write_tsv
from .utils import build_conda_command

# AF 长表格式:位点4列+每样本 AD/DP 两列,刻意不含 GT|AD/DP only, no GT column
AF_QUERY_FORMAT = r"%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n"


def build_af_query_cmd(bcftools_path: str, vcf: str, out_tsv: str) -> List[str]:
    """构造 AF 长表提取命令(conda 包装,不取GT)|build the AF query (no GT)."""
    return build_conda_command(
        bcftools_path, ["query", "-f", AF_QUERY_FORMAT, "-o", out_tsv, vcf])


def load_ad_dp(tsv_path: str, samples: Optional[List[str]] = None) -> HetData:
    """读 AF 长表(仅双等位 SNP;gt 全 -1)|load AF long table (biallelic only).

    兼容 `{sample}_AD` 表头行(test fixture/外部生成)与 bcftools 原生无表头输出
    (样本名取 samples 参数,缺省 S1..SN)。GT 一律填 -1:本分支不依赖 GT。
    |Accepts `{sample}_AD` header or headerless bcftools output; GT is
    hard-coded -1 — this branch never reads GT.
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
                if len(p) > 4 and p[4].endswith("_AD") and p[-1].endswith("_DP"):
                    names = [c[:-3] for c in p[4:] if c.endswith("_AD")]
                    continue   # 表头行|header line
                # 无表头:按列数推断样本数(每样本2列)|headerless: N from columns
                n = (len(p) - 4) // 2
                names = list(samples[:n]) if samples else []
                names += [f"S{i+1}" for i in range(len(names), n)]
            if len(p[2]) != 1 or len(p[3]) != 1 or "," in p[3]:
                continue   # 非双等位SNP(多等位/indel)|not biallelic SNP
            chrom.append(p[0])
            pos.append(int(p[1]))
            row = []
            for i in range(len(names)):
                try:
                    ad, dp = p[4 + i*2], p[5 + i*2]
                except IndexError:
                    break   # 列数不足的截断行:丢弃|truncated row: drop
                adp = ad.split(",")
                ref = int(adp[0]) if adp[0] not in (".", "") else 0
                alt = int(adp[1]) if len(adp) > 1 and adp[1] not in (".", "") else 0
                try:
                    d = int(dp)
                except ValueError:
                    d = ref + alt   # DP 坏值回退 ref+alt|fallback
                row.append((ref, alt, d))
            if len(row) != len(names):
                continue   # 截断行|truncated row
            cols.append(row)
    samples = names
    if cols:
        arr = np.array(cols, dtype=np.int32).transpose(1, 0, 2)   # (样本,位点,3)
    else:
        arr = np.zeros((len(samples), 0, 3), np.int32)
    return HetData(samples=samples,
                   gt=np.full(arr.shape[:2], -1, dtype=np.int8),
                   ref_ad=arr[:, :, 0], alt_ad=arr[:, :, 1], dp=arr[:, :, 2],
                   chrom=chrom, pos=np.array(pos, dtype=np.int32))


def compute_af_based_het(d: HetData, alt_frac_min: float = 0.05,
                         min_depth: int = 10, min_alt_ad: int = 3) -> List[dict]:
    """论文式(Cao et al. 2026)AF直接判定,不依赖GT字段。
    |paper-style (Cao et al. 2026) AF-based calling, GT-independent.

    对每个 depth>=min_depth 的位点:alt_frac = alt_ad/(ref_ad+alt_ad);
    alt_ad>=min_alt_ad 且 alt_frac 落在 [alt_frac_min, 1-alt_frac_min] 之间 → 杂合;
    否则按主导等位判纯合(ref为主→纯合参考型,alt为主→纯合变异型)。
    返回每样本: n_sites_eval, n_het_af, n_hom_alt_af, het_rate_af
    (= n_het_af/(n_het_af+n_hom_alt_af),对齐论文分母口径,不含纯合参考位点)。
    """
    rows = []
    for i, s in enumerate(d.samples):
        tot = d.ref_ad[i].astype(np.int64) + d.alt_ad[i]
        af = np.where(tot > 0, d.alt_ad[i] / np.maximum(tot, 1), 0.0)
        evaluated = d.dp[i] >= min_depth
        het = (evaluated & (d.alt_ad[i] >= min_alt_ad)
               & (af >= alt_frac_min) & (af <= 1 - alt_frac_min))
        hom_alt = evaluated & (af > 1 - alt_frac_min)
        n_het = int(het.sum())
        n_hom_alt = int(hom_alt.sum())
        denom = n_het + n_hom_alt
        rows.append({"sample": s, "n_sites_eval": int(evaluated.sum()),
                     "n_het_af": n_het, "n_hom_alt_af": n_hom_alt,
                     "het_rate_af": n_het / denom if denom else float("nan")})
    return rows


# ---------- 编排层|orchestration ----------

def read_af_table(path: str) -> List[dict]:
    """读回 l1_het_af_based.tsv(断点复用)|read back the AF table (checkpoint)."""
    from pathlib import Path as _P
    p = _P(path)
    if not p.exists():
        return []
    lines = [l.split("\t") for l in p.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 2:
        return []
    head = lines[0]
    rows = []
    for f in lines[1:]:
        r = dict(zip(head, f))
        try:
            r["het_rate_af"] = float(r.get("het_rate_af", ""))
        except ValueError:
            r["het_rate_af"] = None   # 空→None(照 read_verdict_table 口径)|empty -> None
        for k in ("n_sites_eval", "n_het_af", "n_hom_alt_af"):
            try:
                r[k] = int(r.get(k, 0))
            except (TypeError, ValueError):
                r[k] = 0
        rows.append(r)
    return rows


def merge_af_rows(rows: List[dict], af_rows: List[dict]) -> List[dict]:
    """把 het_rate_af 按 sample 并进判读行(只加不改)|merge het_rate_af (add-only).

    AF 侧无数据的样本不加键(该行 het_rate_af 留空),既有 GT-based 判读字段
    原样保留——新分支是纯新增,不改任何既有输出。
    |Samples absent from af_rows simply lack the key; existing GT-based
    fields are untouched (pure addition).
    """
    af = {r["sample"]: r for r in af_rows}
    out = []
    for r in rows:
        a = af.get(r.get("sample"))
        out.append({**r, "het_rate_af": a["het_rate_af"]} if a is not None else r)
    return out


def run_af_het_eval(config, runner, ckpt, vcf: str) -> List[dict]:
    """AF 分支编排:query(AD/DP,不取GT)→判定→写 l1_het_af_based.tsv。
    |AF branch: query (AD/DP, no GT) -> call -> write l1_het_af_based.tsv.

    独立断点(af_query/af_het_eval),支持对旧输出单独补跑本分支;失败返回空表
    (调用方降级,不影响 GT-based 主流程)。数据取自联合 VCF 的 AD/DP——AD 由
    reads 比对直接产生,不含 GT 先验,等效逐 BAM mpileup 口径且免重跑 pileup。
    |Own checkpoints so this branch can be added to an existing run; returns
    [] on failure (caller degrades, GT-based flow unaffected). AD comes from
    the joint VCF (pileup-derived, GT-prior-free; equivalent to per-BAM
    mpileup without re-piling).
    """
    from pathlib import Path
    from .pipeline import _done
    out_dir = Path(config.output_dir) / "04_het_eval"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "af_ad_dp.tsv"
    l1_af = out_dir / "l1_het_af_based.tsv"

    if config.enable_checkpoint and _done(ckpt, "af_het_eval", l1_af):
        runner.logger.info("跳过已完成步骤|Skipping completed step: af_het_eval")
        return read_af_table(str(l1_af))

    # ① AF 长表(独立断点,不动 GT 长表)|AF long table (own checkpoint)
    if config.enable_checkpoint and _done(ckpt, "af_query", tsv):
        runner.logger.info("跳过已完成步骤|Skipping completed step: af_query")
    else:
        runner.logger.info("开始步骤|Starting step: AF长表提取(query AD/DP,不取GT)")
        cmd = build_af_query_cmd(config.bcftools_path, vcf, str(tsv))
        ok, _, _ = runner.run(cmd, "AF长表提取|query AD/DP")
        if not (ok or tsv.exists()):
            runner.logger.error("AF长表提取失败,本分支跳过|AF query failed, branch skipped")
            return []
        if config.enable_checkpoint and ok:
            ckpt.create("af_query")
    if not tsv.exists():
        runner.logger.error(f"AF长表不存在: {tsv}|AF long table missing")
        return []

    # ② 论文式判定+写表|paper-style calling + write
    vcf_samples = list_samples(runner, config.bcftools_path, vcf)
    d = load_ad_dp(str(tsv), vcf_samples or None)
    if len(d.pos) == 0 or len(d.samples) < 1:
        runner.logger.error("AF长表为空或样本数为0,本分支跳过|empty AF table, branch skipped")
        return []
    af_rows = compute_af_based_het(d, alt_frac_min=config.af_het_min_frac,
                                   min_depth=config.af_het_min_depth,
                                   min_alt_ad=config.af_het_min_alt_ad)
    write_tsv(l1_af, af_rows)
    if config.enable_checkpoint:
        ckpt.create("af_het_eval")
    for r in af_rows:
        rate = r["het_rate_af"]
        shown = f"{rate*100:.4f}%" if rate == rate else "n/a"
        runner.logger.info(f"{r['sample']}: AF口径杂合率|AF-based het="
                           f"{shown} (杂合|het={r['n_het_af']} "
                           f"纯合变异|hom_alt={r['n_hom_alt_af']})")
    return af_rows
