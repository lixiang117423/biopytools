"""
annorefine 基因模型构建 + 质控(GFF3)|Model construction + QC (GFF3)
miniprot 命中 → GFF3 基因模型(方案A, CDS级) + 质控过滤
|miniprot hits → GFF3 gene models (scheme A, CDS-level) + QC filter
"""

from typing import List, Dict, Tuple, Optional

from .evidence import MiniprotHit, hit_key
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
        expression: {hit_key(hit): (mean_depth, breadth%)}, 给定则做表达过滤②⑤|if given, do expression filter
    """
    from .evidence import has_complete_orf
    require_orf = getattr(config, 'require_real_orf', True) and genome is not None
    use_expr = bool(expression)   # None(失败/无BAM)或空dict 都跳过表达过滤, 不静默丢|None/{} => skip
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
            depth, breadth = expression.get(hit_key(h), (0.0, 0.0))
            if depth < config.min_expression_depth or breadth < config.min_coverage_breadth:
                continue
        # 真 TE 区排除(可选, 默认不排: 真基因可能在 TE 区)
        # |TE exclusion (optional, default off: real genes may reside in TE regions)
        if getattr(config, 'exclude_te_gap', False) and repeat_regions:
            te = repeat_regions.get(h.chrom, [])
            te_intervals = [(s, e) for s, e, *_ in te]   # 去掉 family, 只留区间
            hit_cds = [(s, e) for s, e, _ in h.cds_exons]
            if te_intervals and cds_overlap_ratio(hit_cds, te_intervals) > config.te_overlap_cutoff:
                continue
        passed.append(h)
    return passed


def qc_filter_small(hits: List[MiniprotHit], config,
                    repeat_regions: Dict[str, List[Tuple[int, int]]],
                    genome: Optional[dict] = None,
                    expression: Optional[dict] = None
                    ) -> List[MiniprotHit]:
    """
    小蛋白通道质控(通用)|small-protein lane QC (general, no functional assumption)
    放宽长度, 守 完整ORF① + 同源 + 表达②⑤ + TE 排除。无表达数据时退化为 ORF+更严同源。
    |Relaxes length; gates on real-ORF + homology + expression + TE. Degrades to
    ORF + stricter homology (normal 70/80) when no expression data.

    ⚠ 设计偏离(已记录)|design deviation (logged):
    不应用常规通道的 require_complete_orf(丢 coverage<99, 会抵消 small_min_coverage 放宽);
    完整性由 require_real_orf(翻译校验①)把守。|Does NOT apply require_complete_orf
    (drops cov<99, negates small_min_coverage); completeness via require_real_orf.

    无信号肽/无物种假设 —— 普适。|No signal-peptide / species assumption — general.
    """
    from .evidence import has_complete_orf
    require_orf = getattr(config, 'require_real_orf', True) and genome is not None
    has_expr = expression is not None      # None = 无 BAM/depth 失败 → 退化模式|degraded
    # 有表达: 放宽同源(small 50/50); 无表达(退化): 常规严同源(gap 70/80)
    # |with expression: relaxed (50/50); degraded: normal-grade (70/80)
    min_id = config.small_min_identity if has_expr else config.gap_min_identity
    min_cov = config.small_min_coverage if has_expr else config.gap_min_coverage
    passed = []
    for h in hits:
        if h.identity < min_id:
            continue
        if h.coverage < min_cov:
            continue
        cds_len = sum(e - s + 1 for s, e, _ in h.cds_exons)
        # 长度区间: < gap_min_cds_len(被常规挡掉) 且 <= small_max_cds_len(上限)
        # |length window: < gap_min_cds_len and <= small_max_cds_len
        if cds_len >= config.gap_min_cds_len or cds_len > config.small_max_cds_len:
            continue
        # ① 真实完整 ORF(不可放宽, 核心 completeness 门)|real complete ORF (not relaxed)
        if require_orf and not has_complete_orf(h, genome):
            continue
        # ②⑤ 表达(有数据时必需)|expression (required when data available)
        if has_expr:
            depth, breadth = expression.get(hit_key(h), (0.0, 0.0))
            if depth < config.small_min_expression_depth \
                    or breadth < config.small_min_coverage_breadth:
                continue
        # 强制排 TE(默认 small_exclude_te=True)|force TE exclusion (default on)
        if getattr(config, 'small_exclude_te', True) and repeat_regions:
            te = repeat_regions.get(h.chrom, [])
            te_intervals = [(s, e) for s, e, *_ in te]
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
