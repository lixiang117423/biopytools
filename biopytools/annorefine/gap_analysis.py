"""
annorefine 漏检/合并判定(GFF3 + 全 prot 普适)|Gap & merged-gene detection
对比 miniprot 命中与 braker.gff3, 找漏检(证据有·基因无)和错误合并(1基因盖多拷贝)
|Compare miniprot hits vs braker.gff3: find missing copies & merged genes

全 prot 场景关键:命中去重(同位置多 query 合并) + 合并判定按 query 分组
|Multi-query: dedupe hits at same locus + per-query merged-gene detection
"""

import bisect
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Tuple

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
            # 坐标解析(畸形行跳过, 不崩)|parse coords, skip malformed
            try:
                s_i, e_i = int(start), int(end)
            except ValueError:
                continue
            if feat == 'gene':
                gid = attrs.get('ID', '')
                if gid:
                    genes[gid] = BrakerGene(
                        gene_id=gid, chrom=chrom,
                        start=s_i, end=e_i, strand=strand)
            elif feat in ('mRNA', 'transcript'):
                mid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                if mid and parent:
                    mrna_to_gene[mid] = parent
            elif feat == 'CDS':
                parent = attrs.get('Parent', '')
                gene_id = mrna_to_gene.get(parent, parent)
                cds_temp.setdefault(gene_id, []).append(
                    (s_i, e_i, chrom, strand))

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


def _build_gene_span_index(
        braker_genes: Dict[str, BrakerGene]
        ) -> Dict[str, Tuple[List[int], List[BrakerGene]]]:
    """按 chrom 建基因 span 索引(按 start 升序)|per-chrom gene index, sorted by start
    用于 detect_gaps 的 span 相交查询, 把 O(hits×genes) 暴力扫描降到按染色体+前缀 bisect。
    |Turns detect_gaps' O(hits×genes) scan into per-chrom bisect prefix scan.
    假设 CDS ⊆ gene span(有效 GFF3 成立; BRAKER 输出满足)|assumes CDS ⊆ gene span.
    """
    by_chrom: Dict[str, List[BrakerGene]] = {}
    for g in braker_genes.values():
        by_chrom.setdefault(g.chrom, []).append(g)
    idx: Dict[str, Tuple[List[int], List[BrakerGene]]] = {}
    for chrom, gs in by_chrom.items():
        gs_sorted = sorted(gs, key=lambda g: g.start)
        idx[chrom] = ([g.start for g in gs_sorted], gs_sorted)
    return idx


def _genes_intersecting_span(idx, chrom: str, qstart: int, qend: int
                             ) -> Iterator[BrakerGene]:
    """yield span 与 [qstart,qend] 相交的基因(按 start 升序)
    |yield genes whose span intersects [qstart,qend] (inclusive, start-sorted)
    相交条件: gene.start <= qend 且 gene.end >= qstart
    |intersect iff gene.start <= qend and gene.end >= qstart
    """
    if chrom not in idx:
        return
    starts, gs = idx[chrom]
    hi = bisect.bisect_right(starts, qend)   # start<=qend 的前缀|prefix with start<=qend
    for i in range(hi):
        if gs[i].end >= qstart:
            yield gs[i]


def detect_gaps(hits: List[MiniprotHit],
                braker_genes: Dict[str, BrakerGene],
                overlap_cutoff: float,
                coord_zero_overlap: bool = False) -> List[MiniprotHit]:
    """
    找漏检命中|Find missing-copy hits
    - coord_zero_overlap=False(默认): 与所有 braker 基因 CDS 重叠 < overlap_cutoff% 的命中
      |CDS overlap with all braker genes < cutoff
    - coord_zero_overlap=True: 与任一 braker 基因坐标 span 有交集即算已覆盖(更严: 真新基因
      不应与现有基因有任何坐标重叠)|any genomic-span intersection with a braker gene => covered

    实现用 per-chrom 基因 span 索引(span 相交预过滤), cds 模式下再对候选做精确 CDS 重叠判定,
    结果与原 O(hits×genes) 扫描完全一致(假设 CDS ⊆ gene span)。
    |Uses a per-chrom span index for fast intersect prefilter; cds mode then does
    exact CDS-overlap on candidates. Identical to the old brute-force (CDS ⊆ span).
    """
    idx = _build_gene_span_index(braker_genes)
    gaps = []
    for hit in hits:
        is_covered = False
        for gene in _genes_intersecting_span(idx, hit.chrom, hit.start, hit.end):
            if coord_zero_overlap:
                # span 相交即已覆盖(任意交集)|genomic span intersection => covered
                is_covered = True
                break
            else:
                hit_cds = [(s, e) for s, e, _ in hit.cds_exons]
                if cds_overlap_ratio(hit_cds, gene.cds_intervals) > overlap_cutoff:
                    is_covered = True
                    break
        if not is_covered:
            gaps.append(hit)
    return gaps


def _build_hit_start_index(
        hits: List[MiniprotHit]
        ) -> Dict[str, Tuple[List[int], List[MiniprotHit]]]:
    """按 chrom 建命中索引(按 start 升序)|per-chrom hit index, sorted by start
    用于 detect_merged_genes 的包含查询, 避免每基因扫描全量 hits。
    |Avoids scanning all hits per gene in detect_merged_genes.
    """
    by_chrom: Dict[str, List[MiniprotHit]] = {}
    for h in hits:
        by_chrom.setdefault(h.chrom, []).append(h)
    idx: Dict[str, Tuple[List[int], List[MiniprotHit]]] = {}
    for chrom, hs in by_chrom.items():
        hs_sorted = sorted(hs, key=lambda h: h.start)
        idx[chrom] = ([h.start for h in hs_sorted], hs_sorted)
    return idx


def _hits_with_start_in_range(idx, chrom: str, lo_pos: int, hi_pos: int
                              ) -> Iterator[MiniprotHit]:
    """yield chrom 上 start ∈ [lo_pos,hi_pos] 的命中(按 start 升序)
    |yield hits with start in [lo_pos,hi_pos] (start-sorted)"""
    if chrom not in idx:
        return
    starts, hs = idx[chrom]
    lo = bisect.bisect_left(starts, lo_pos)
    hi = bisect.bisect_right(starts, hi_pos)
    for i in range(lo, hi):
        yield hs[i]


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

    实现用 per-chrom 命中 start 索引取 gene 内候选(start ∈ [gene.start,gene.end]),
    再过滤 end<=gene.end; 结果集与原实现一致, 并按输入顺序还原(保证确定性输出)。
    |Uses a per-chrom hit start-index to fetch in-gene candidates, then filters
    end<=gene.end. Same result set as brute-force; input order restored for determinism.
    """
    hit_idx = _build_hit_start_index(hits)
    hit_order = {id(h): i for i, h in enumerate(hits)}   # 还原输入顺序|restore input order
    merged = []
    for gene in braker_genes.values():
        # 包含在 gene 内: start ∈ [gene.start, gene.end] 且 end <= gene.end
        # (start<=gene.end 由 end<=gene.end 隐含; range 限定保证与原实现同集)
        # |contained iff start in [gene.start,gene.end] and end<=gene.end
        contained = [h for h in _hits_with_start_in_range(
                        hit_idx, gene.chrom, gene.start, gene.end)
                     if h.end <= gene.end]
        if len(contained) < split_min_hits:
            continue
        contained.sort(key=lambda h: hit_order[id(h)])   # 与原扫描顺序一致|match original order
        # 按 query 分组|group by query
        by_query: Dict[str, List[MiniprotHit]] = {}
        for h in contained:
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
            full_all = [h for h in contained if h.coverage >= split_min_copy_coverage]
            if full_all:
                merged.append((gene, full_all))
    return merged


def parse_repeat_out(repeat_out: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """解析 RepeatMasker .out → {chrom: [(start, end, family)]}
    |Parse RepeatMasker .out (cols[4]=chrom, [5,6]=begin/end, [10]=class/family)"""
    regions: Dict[str, List[Tuple[int, int, str]]] = {}
    if not repeat_out or not os.path.exists(repeat_out):
        return regions
    with open(repeat_out) as f:
        for line in f:
            if not line.strip():
                continue
            stripped = line.lstrip()
            if stripped.startswith('SW') or stripped.startswith('score'):
                continue   # 跳过 header(2 行)|skip header
            cols = line.split()
            if len(cols) < 7:
                continue
            try:
                chrom = cols[4]
                b, e = int(cols[5]), int(cols[6])
                begin, end = (b, e) if b <= e else (e, b)
                # class/family 在 cols[10](标准 .out, 如 LTR/Gypsy, Simple_repeat)
                family = cols[10] if len(cols) > 10 else ''
                regions.setdefault(chrom, []).append((begin, end, family))
            except (ValueError, IndexError):
                continue
    return regions
