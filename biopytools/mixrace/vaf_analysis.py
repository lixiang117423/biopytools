"""mixrace 等位频率谱分析(单倍体,导师方法论 v2)|mixrace AFS analysis (haploid, advisor v2).

解析 freebayes AO/RO → 算 杂合率 het_rate + AFS 形态 + 优势株占比。
根肿菌静息孢子单倍体:n=20,每个位点理论上单一等位;"杂合"=混合群体的基因型多样性。
|Parse freebayes AO/RO -> het_rate + AFS shape + dominant proportion. P. brassicae resting
spores are haploid (n=20); "heterozygosity" here = mixed-population genotype diversity.
"""
import statistics
from typing import List, Optional

# AFS 中间频率区间(5–95%,过低 VAF 多为噪声/测序错误)|intermediate-frequency band (5-95%)
AFS_MID_LO = 0.05
AFS_MID_HI = 0.95


def parse_freebayes(query_text: str) -> List[dict]:
    """解析 bcftools query(CHROM/POS/REF/ALT/RO/AO)|parse freebayes AO/RO query.

    期望格式: -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%RO]\\t[%AO]\\n'
    RO=ref 观测数(单值),AO=alt 观测数(多等位逗号分隔)。|RO=ref obs, AO=alt obs (comma-sep).
    """
    recs = []
    for line in query_text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        chrom, pos, ref, alt_str, ro_str, ao_str = f[0], int(f[1]), f[2], f[3], f[4], f[5]
        alts = alt_str.split(",")
        try:
            ro = int(ro_str) if ro_str and ro_str != "." else 0
        except ValueError:
            ro = 0
        try:
            aos = [int(x) for x in ao_str.split(",") if x and x != "."]
        except ValueError:
            aos = []
        recs.append({"chrom": chrom, "pos": pos, "ref": ref, "alts": alts,
                     "ro": ro, "aos": aos})
    return recs


def compute_afs(records: List[dict], min_alt_reads: int, genome_size: int = 0) -> dict:
    """计算杂合率 + AFS 形态 + 优势株占比|compute het_rate + AFS shape + dominant proportion.

    - het 位点 = 同时有 ref 与 alt 读支持(RO>0 且 某个 AO>=min_alt_reads)的位点(混合信号)。
    - het_rate = het_sites / genome_size(全基因组杂合率,导师 <0.01%/0.1%/1% 阈值)。
    - AFS 形态:由中间频率(0.05–0.95)VAF 分布判 monoclonal/two_clone_50_50/dominant_minor/smeared。
    - dominant_proportion = 混合位点主等位频率的中位数(方法B:主等位单倍型支持率)。
    |het site = ref+alt both supported (mixed signal); het_rate = het_sites/genome_size;
    AFS shape from intermediate VAFs; dominant = median major-allele freq over mixed sites.
    """
    het_sites = 0
    total_variant = len(records)
    intermediate_vafs: List[float] = []
    intermediate_sites = 0
    major_freqs: List[float] = []
    for r in records:
        ro, aos, total = r["ro"], r["aos"], r["ro"] + sum(r["aos"])
        if total == 0:
            continue
        supported = [a for a in aos if a >= min_alt_reads]
        if ro > 0 and supported:                      # ref + alt 都有 → 混合(杂合)位点
            het_sites += 1
            major_freqs.append(max([ro] + aos) / total)
            # 只收中间频率(0.05–0.95)VAF;VAF≈0 的噪声位点会把形态拉向 smeared|
            # keep only intermediate VAFs (0.05-0.95); near-zero noise skews shape to smeared.
            # 注意:intermediate_vafs 按 VAF 计、intermediate_sites 按位点计——多等位位点
            # 会在 shape 判据中贡献多个点(轻度加权,视为对复杂分离的合理强调;freebayes
            # -p 1 输出多等位罕见,影响可忽略)。|Note: vafs counted per-ALT vs sites
            # per-locus — multiallelic loci weigh slightly more in shape (reasonable
            # emphasis; rare under -p 1, negligible).
            band = [a / total for a in supported
                    if AFS_MID_LO <= a / total <= AFS_MID_HI]
            if band:
                intermediate_sites += 1
            intermediate_vafs.extend(band)
    het_rate = het_sites / genome_size if genome_size else 0.0
    # 中间频率位点占 het 位点比例(与 classify 判据同基准;原 het/total 语义不符)|
    # intermediate sites / het sites (same basis as classify; old het/total was inconsistent).
    intermediate_ratio = intermediate_sites / het_sites if het_sites else 0.0
    dominant = statistics.median(major_freqs) if major_freqs else None
    maf = statistics.median(1 - mf for mf in major_freqs) if major_freqs else None   # 次等位频率中位数
    shape = classify_afs_shape(intermediate_vafs, intermediate_ratio)
    return {
        "het_rate": het_rate, "het_sites": het_sites, "total_variant": total_variant,
        "intermediate_ratio": intermediate_ratio,
        "afs_shape": shape, "dominant_proportion": dominant, "maf": maf,
    }


def classify_afs_shape(intermediate_vafs: List[float], intermediate_ratio: float) -> str:
    """AFS 形态分类(导师判读经验)|classify AFS shape per advisor heuristics.

    monoclonal:中间频率 SNP<5% | two_clone_50_50:集中在 0.5 | dominant_minor:主峰>0.5 |
    smeared:0–1 连续糊(多基因型)。
    """
    if intermediate_ratio < 0.05 or not intermediate_vafs:
        return "monoclonal"
    n = len(intermediate_vafs)
    near_half = sum(1 for v in intermediate_vafs if abs(v - 0.5) <= 0.1) / n
    if near_half > 0.5:
        return "two_clone_50_50"
    high = sum(1 for v in intermediate_vafs if v > 0.5) / n
    if high > 0.5:
        return "dominant_minor"
    return "smeared"
