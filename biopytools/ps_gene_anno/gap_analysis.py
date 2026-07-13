"""
ps-gene-anno 漏检/合并判定|Gap & merged-gene detection (GFF3 流)
对比 miniprot 命中与 braker.gff3, 找漏检(证据有·基因无)和错误合并(1基因盖多拷贝)
|Compare miniprot hits vs braker.gff3: find missing copies & merged genes
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
    通过 mRNA→gene 映射链, 把 CDS 归到 gene, 合并同 gene 所有 mRNA 的 CDS。
    |Chain CDS→mRNA→gene via Parent, merge CDS across mRNAs per gene.
    """
    genes: Dict[str, BrakerGene] = {}
    mrna_to_gene: Dict[str, str] = {}   # mRNA ID → gene ID
    cds_temp: Dict[str, List[Tuple[int, int, str, str]]] = {}  # gene_id → [(start,end,chrom,strand)]

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
                parent = attrs.get('Parent', '')   # mRNA ID
                gene_id = mrna_to_gene.get(parent, parent)   # 链到 gene
                cds_temp.setdefault(gene_id, []).append(
                    (int(start), int(end), chrom, strand))

    # 合并 CDS 到 gene|Attach merged CDS to genes
    for gid, intervals in cds_temp.items():
        if not gid:
            continue
        if gid not in genes:
            # 无 gene 行, 从 CDS 推断边界|No gene line, infer bounds from CDS
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
    """
    hit CDS 被 gene CDS 覆盖的比例(%)|hit CDS overlap with gene CDS (%)
    = 交集长度 / hit CDS 总长度 × 100
    """
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


def pairwise_no_cds_overlap(hits: List[MiniprotHit]) -> bool:
    """所有命中两两 CDS 不重叠 → True|True if no two hits' CDS overlap"""
    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            for s1, e1, _ in hits[i].cds_exons:
                for s2, e2, _ in hits[j].cds_exons:
                    if min(e1, e2) - max(s1, s2) + 1 > 0:
                        return False
    return True


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
    找错误合并基因: 一个 braker 基因含 ≥split_min_hits 个互相不重叠的完整蛋白拷贝
    |Find merged genes: a braker gene covering >= N independent full copies

    保守条件: 每个命中 coverage >= split_min_copy_coverage(完整独立拷贝)
    |Conservative: each hit must cover its protein >= cutoff (independent full copy)
    """
    merged = []
    for gene in braker_genes.values():
        hits_in = [h for h in hits
                   if h.chrom == gene.chrom
                   and h.start >= gene.start and h.end <= gene.end]
        if len(hits_in) < split_min_hits:
            continue
        if (all(h.coverage >= split_min_copy_coverage for h in hits_in)
                and pairwise_no_cds_overlap(hits_in)):
            merged.append((gene, hits_in))
    return merged


def parse_repeat_out(repeat_out: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    解析 RepeatMasker .out → {chrom: [(start,end)]}(1-based inclusive)
    |Parse RepeatMasker .out to TE regions
    RepeatMasker .out 列: score div del ins sequence begin end ...
    """
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
