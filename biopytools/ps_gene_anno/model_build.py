"""
ps-gene-anno 基因模型构建 + 质控(GFF3)|Model construction + QC (GFF3)
miniprot 命中 → GFF3 基因模型(方案A, CDS级) + 质控过滤
|miniprot hits → GFF3 gene models (scheme A, CDS-level) + QC filter
"""

from typing import List, Dict, Tuple, Optional

from .evidence import MiniprotHit
from .gap_analysis import cds_overlap_ratio


def qc_filter(hits: List[MiniprotHit], config,
              repeat_regions: Dict[str, List[Tuple[int, int]]],
              genome: Optional[dict] = None,
              expression: Optional[dict] = None
              ) -> List[MiniprotHit]:
    """
    质控过滤命中|QC-filter hits
    检查: identity/coverage/cds_len + 真实完整ORF(①) + 表达证据(②⑤, 唯一reads深度+广度) + 真TE区(可选)
    |identity/coverage/cds_len + real complete ORF + expression (depth+breadth) + TE (optional)

    Args:
        genome: {chrom: seq}, 给定则做真实 ORF 检查①|if given, do real-ORF check
        expression: {id(hit): (mean_depth, breadth%)}, 给定则做表达过滤②⑤|if given, do expression filter
    """
    from .evidence import has_complete_orf
    require_orf = getattr(config, 'require_real_orf', True) and genome is not None
    use_expr = expression is not None
    passed = []
    for h in hits:
        if h.identity < config.gap_min_identity:
            continue
        if h.coverage < config.gap_min_coverage:
            continue
        cds_len = sum(e - s + 1 for s, e, _ in h.cds_exons)
        if cds_len < config.gap_min_cds_len:
            continue
        # 完整 ORF: coverage 接近完整(≥99)|complete ORF: near-full coverage (miniprot)
        if config.require_complete_orf and h.coverage < 99:
            continue
        # ① 真实完整 ORF: 3倍数 + ATG + 终止密码子(翻译验证, 比 miniprot 覆盖率更严)
        # |real complete ORF: 3×len + ATG + stop (translation check)
        if require_orf and not has_complete_orf(h, genome):
            continue
        # ②⑤ 表达证据: 唯一 reads 平均深度 + 覆盖广度
        # |expression: unique-read mean depth + coverage breadth
        if use_expr:
            depth, breadth = expression.get(id(h), (0.0, 0.0))
            if depth < config.min_expression_depth or breadth < config.min_coverage_breadth:
                continue
        # 真 TE 区排除(可选, 默认不排: 疫霉效应子常在 TE 区)
        # |TE exclusion (optional, default off: oomycete effectors often TE-rich)
        if getattr(config, 'exclude_te_gap', False) and repeat_regions:
            te = repeat_regions.get(h.chrom, [])
            te_intervals = [(s, e) for s, e, *_ in te]   # 去掉 family, 只留区间
            hit_cds = [(s, e) for s, e, _ in h.cds_exons]
            if te_intervals and cds_overlap_ratio(hit_cds, te_intervals) > config.te_overlap_cutoff:
                continue
        passed.append(h)
    return passed


def build_gene_models(hits: List[MiniprotHit], prefix: str) -> List[str]:
    """
    miniprot 命中 → GFF3 行(gene/mRNA/exon/CDS + ID/Parent 层级)|hits → GFF3 lines
    ID: {prefix}_gap_{N} / {prefix}_gap_{N}.t1
    """
    lines = ["##gff-version 3"]
    for i, h in enumerate(hits, start=1):
        gid = f"{prefix}_gap_{i}"
        tid = f"{gid}.t1"
        # gene 行|gene line
        lines.append(
            f'{h.chrom}\tpsfill\tgene\t{h.start}\t{h.end}\t.\t'
            f'{h.strand}\t.\tID={gid};')
        # mRNA 行|mRNA line
        lines.append(
            f'{h.chrom}\tpsfill\tmRNA\t{h.start}\t{h.end}\t.\t'
            f'{h.strand}\t.\tID={tid};Parent={gid};')
        for j, (s, e, ph) in enumerate(h.cds_exons, start=1):
            # exon 行|exon line
            lines.append(
                f'{h.chrom}\tpsfill\texon\t{s}\t{e}\t.\t'
                f'{h.strand}\t.\tID={tid}.exon{j};Parent={tid};')
            # CDS 行|CDS line (第8列 phase)
            lines.append(
                f'{h.chrom}\tpsfill\tCDS\t{s}\t{e}\t.\t'
                f'{h.strand}\t{ph}\tID={tid}.cds{j};Parent={tid};')
    return lines
