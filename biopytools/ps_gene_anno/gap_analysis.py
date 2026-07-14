"""
ps-gene-anno 漏检/合并判定(GFF3 + 全 prot 普适)|Gap & merged-gene detection
对比 miniprot 命中与 braker.gff3, 找漏检(证据有·基因无)和错误合并(1基因盖多拷贝)
|Compare miniprot hits vs braker.gff3: find missing copies & merged genes

全 prot 场景关键:命中去重(同位置多 query 合并) + 合并判定按 query 分组
|Multi-query: dedupe hits at same locus + per-query merged-gene detection
"""

import os
import re
from dataclasses import dataclass, field
from typing import List, Tuple, Dict

from .evidence import MiniprotHit


@dataclass
class BrakerGene:
    """braker.gff3 基因(合并所有 mRNA 的 CDS)|braker gene with merged CDS"""
    gene_id: str
    chrom: str
    start: int                      # gene 边界 1-based inclusive
    end: int
    strand: str
    cds_intervals: List[Tuple[int, int]] = field(default_factory=list)


def _parse_gff3_attr(attr_str: str) -> dict:
    """解析 GFF3 第9列(ID=g1;Parent=g2)|Parse GFF3 column 9 (key=value)"""
    attrs = {}
    for kv in attr_str.split(';'):
        kv = kv.strip()
        if '=' in kv:
            k, v = kv.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_braker_gff3(gff3_path: str) -> Dict[str, BrakerGene]:
    """
    解析 braker.gff3 → {gene_id: BrakerGene}|Parse braker.gff3 to genes
    GFF3 层级: gene(ID) → mRNA(Parent=gene) → CDS(Parent=mRNA)
    """
    genes: Dict[str, BrakerGene] = {}
    mrna_to_gene: Dict[str, str] = {}
    cds_temp: Dict[str, List[Tuple[int, int, str, str]]] = {}

    if not os.path.exists(gff3_path):
        return genes

    with open(gff3_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom, src, feat, start, end, score, strand, phase, attr = cols[:9]
            attrs = _parse_gff3_attr(attr)
            if feat == 'gene':
                gid = attrs.get('ID', '')
                if gid:
                    genes[gid] = BrakerGene(
                        gene_id=gid, chrom=chrom,
                        start=int(start), end=int(end), strand=strand)
            elif feat in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                if mid and parent:
                    mrna_to_gene[mid] = parent
            elif feat == 'CDS':
                parent = attrs.get('Parent', '')
                gene_id = mrna_to_gene.get(parent, parent)
                cds_temp.setdefault(gene_id, []).append(
                    (int(start), int(end), chrom, strand))

    for gid, intervals in cds_temp.items():
        if not gid:
            continue
        if gid not in genes:
            chroms = {iv[2] for iv in intervals}
            strands = {iv[3] for iv in intervals}
            genes[gid] = BrakerGene(
                gene_id=gid, chrom=next(iter(chroms)) if chroms else '',
                start=min(iv[0] for iv in intervals),
                end=max(iv[1] for iv in intervals),
                strand=next(iter(strands)) if strands else '+')
        unique = sorted(set((iv[0], iv[1]) for iv in intervals))
        genes[gid].cds_intervals = unique
    return genes


def cds_overlap_ratio(hit_cds: List[Tuple[int, int]],
                      gene_cds: List[Tuple[int, int]]) -> float:
    """hit CDS 被 gene CDS 覆盖的比例(%)|hit CDS overlap (%)"""
    if not hit_cds:
        return 0.0
    hit_len = sum(e - s + 1 for s, e in hit_cds)
    if hit_len == 0:
        return 0.0
    overlap = 0
    for hs, he in hit_cds:
        for gs, ge in gene_cds:
            ov = min(he, ge) - max(hs, gs) + 1
            if ov > 0:
                overlap += ov
    return overlap / hit_len * 100.0


def _cds_total_len(hit: MiniprotHit) -> int:
    """hit 的 CDS 总长|total CDS length"""
    return sum(e - s + 1 for s, e, _ in hit.cds_exons)


def _cds_overlap_len(a: MiniprotHit, b: MiniprotHit) -> int:
    """两个 hit 的 CDS 重叠碱基数|CDS overlap length between two hits"""
    ov = 0
    for s1, e1, _ in a.cds_exons:
        for s2, e2, _ in b.cds_exons:
            ov += max(0, min(e1, e2) - max(s1, s2) + 1)
    return ov


def pairwise_no_cds_overlap(hits: List[MiniprotHit]) -> bool:
    """所有命中两两 CDS 不重叠 → True|True if no two hits' CDS overlap"""
    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            if _cds_overlap_len(hits[i], hits[j]) > 0:
                return False
    return True


def dedupe_hits(hits: List[MiniprotHit], overlap_ratio: float = 0.5
                ) -> List[MiniprotHit]:
    """
    同位置多 query 命中去重(CDS 重叠>overlap_ratio 合并, 保留 identity/coverage 最高)
    |Dedup multi-query hits at same locus

    全 prot 场景:多个蛋白 query 命中同一基因组位置(如 avr1a 拷贝被 Avr1a +
    Psojae_XP_* + G9540 同时命中), 不去重会导致重复补基因 + 合并判定误判重叠。
    |Multi-query hits at same locus must be merged to one.
    """
    sorted_hits = sorted(hits, key=lambda h: (h.chrom, h.strand, h.start))
    deduped: List[MiniprotHit] = []
    for h in sorted_hits:
        merged_into = False
        if deduped:
            last = deduped[-1]
            if last.chrom == h.chrom and last.strand == h.strand:
                ov = _cds_overlap_len(h, last)
                min_len = min(_cds_total_len(h), _cds_total_len(last))
                if min_len > 0 and ov / min_len >= overlap_ratio:
                    # 同位置, 保留 identity/coverage 更高的|keep best
                    if (h.identity, h.coverage) > (last.identity, last.coverage):
                        deduped[-1] = h
                    merged_into = True
        if not merged_into:
            deduped.append(h)
    return deduped


def detect_gaps(hits: List[MiniprotHit],
                braker_genes: Dict[str, BrakerGene],
                overlap_cutoff: float) -> List[MiniprotHit]:
    """
    找漏检命中: 与所有 braker 基因 CDS 重叠 < overlap_cutoff% 的命中
    |Find missing-copy hits: CDS overlap with all braker genes < cutoff
    """
    gaps = []
    gene_list = list(braker_genes.values())
    for hit in hits:
        is_covered = False
        for gene in gene_list:
            if gene.chrom != hit.chrom:
                continue
            hit_cds = [(s, e) for s, e, _ in hit.cds_exons]
            if cds_overlap_ratio(hit_cds, gene.cds_intervals) >= overlap_cutoff:
                is_covered = True
                break
        if not is_covered:
            gaps.append(hit)
    return gaps


def detect_merged_genes(
        hits: List[MiniprotHit],
        braker_genes: Dict[str, BrakerGene],
        split_min_hits: int,
        split_min_copy_coverage: float
        ) -> List[Tuple[BrakerGene, List[MiniprotHit]]]:
    """
    找错误合并基因: 按 query 分组, 任一 query 在 gene 内含 ≥split_min_hits 个
    互相不重叠的完整拷贝 → 判合并|Find merged genes (per-query detection)

    全 prot 场景:不同 query 命中同一基因不同区, 混合 pairwise 会误判重叠。
    改为按 query 分组: 同一 query 的多个完整独立拷贝才算多拷贝合并。
    |Per-query: only same query's >=N independent full copies count.
    """
    merged = []
    for gene in braker_genes.values():
        hits_in = [h for h in hits
                   if h.chrom == gene.chrom
                   and h.start >= gene.start and h.end <= gene.end]
        if len(hits_in) < split_min_hits:
            continue
        # 按 query 分组|group by query
        by_query: Dict[str, List[MiniprotHit]] = {}
        for h in hits_in:
            by_query.setdefault(h.query_id, []).append(h)
        # 任一 query 在 gene 内 ≥N 完整独立拷贝 → 合并
        is_merged = False
        for _q, qhits in by_query.items():
            full = [h for h in qhits if h.coverage >= split_min_copy_coverage]
            if len(full) >= split_min_hits and pairwise_no_cds_overlap(full):
                is_merged = True
                break
        if is_merged:
            # 拆分用 gene 内所有完整拷贝(去重后)|all full copies for split models
            full_all = [h for h in hits_in if h.coverage >= split_min_copy_coverage]
            if full_all:
                merged.append((gene, full_all))
    return merged


def parse_repeat_out(repeat_out: str) -> Dict[str, List[Tuple[int, int]]]:
    """解析 RepeatMasker .out → {chrom: [(start,end)]}|Parse RepeatMasker .out"""
    regions: Dict[str, List[Tuple[int, int]]] = {}
    if not repeat_out or not os.path.exists(repeat_out):
        return regions
    with open(repeat_out) as f:
        for line in f:
            if line.startswith('SW') or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 7:
                continue
            try:
                chrom = cols[4]
                b, e = int(cols[5]), int(cols[6])
                begin, end = (b, e) if b <= e else (e, b)
                regions.setdefault(chrom, []).append((begin, end))
            except (ValueError, IndexError):
                continue
    return regions
