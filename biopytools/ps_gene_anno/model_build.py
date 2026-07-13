"""
ps-gene-anno 基因模型构建 + 质控(GFF3)|Model construction + QC (GFF3)
miniprot 命中 → GFF3 基因模型(方案A, CDS级) + 质控过滤
|miniprot hits → GFF3 gene models (scheme A, CDS-level) + QC filter
"""

from typing import List, Dict, Tuple

from .evidence import MiniprotHit
from .gap_analysis import cds_overlap_ratio


def qc_filter(hits: List[MiniprotHit], config,
              repeat_regions: Dict[str, List[Tuple[int, int]]]
              ) -> List[MiniprotHit]:
    """
    质控过滤命中|QC-filter hits
    检查: identity/coverage/cds_len/complete_orf/真TE区
    |Checks: identity/coverage/cds_len/complete_orf/real-TE overlap
    """
    passed = []
    for h in hits:
        if h.identity < config.gap_min_identity:
            continue
        if h.coverage < config.gap_min_coverage:
            continue
        cds_len = sum(e - s + 1 for s, e, _ in h.cds_exons)
        if cds_len < config.gap_min_cds_len:
            continue
        # 完整 ORF: coverage 接近完整(≥99)|complete ORF: near-full coverage
        if config.require_complete_orf and h.coverage < 99:
            continue
        # 真 TE 区排除|real-TE exclusion
        if repeat_regions:
            te = repeat_regions.get(h.chrom, [])
            hit_cds = [(s, e) for s, e, _ in h.cds_exons]
            if cds_overlap_ratio(hit_cds, te) > config.te_overlap_cutoff:
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
            f'{h.chrom}\tps_gene_anno\tgene\t{h.start}\t{h.end}\t.\t'
            f'{h.strand}\t.\tID={gid};')
        # mRNA 行|mRNA line
        lines.append(
            f'{h.chrom}\tps_gene_anno\tmRNA\t{h.start}\t{h.end}\t.\t'
            f'{h.strand}\t.\tID={tid};Parent={gid};')
        for j, (s, e, ph) in enumerate(h.cds_exons, start=1):
            # exon 行|exon line
            lines.append(
                f'{h.chrom}\tps_gene_anno\texon\t{s}\t{e}\t.\t'
                f'{h.strand}\t.\tID={tid}.exon{j};Parent={tid};')
            # CDS 行|CDS line (第8列 phase)
            lines.append(
                f'{h.chrom}\tps_gene_anno\tCDS\t{s}\t{e}\t.\t'
                f'{h.strand}\t{ph}\tID={tid}.cds{j};Parent={tid};')
    return lines
