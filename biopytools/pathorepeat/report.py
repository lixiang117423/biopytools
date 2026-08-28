"""pathorepeat 汇总报告模块|pathorepeat Report Module

repeat_summary(总体+逐类 bp/%/GC%) + families_classified(家族分类表)
|repeat_summary (overall + per-class bp/%/GC%) + families_classified
"""

import os
from typing import Dict, List, Optional

from .utils import (RepeatHit, fasta_ids, find_overlaps,
                    nearest_repeat_distance, overlap_bp)

SUMMARY_HEADER = 'category\tclass\tbp\tpct_genome\tn_families\tgc_pct'
FAMILIES_HEADER = 'family\torder\tsuperfamily\tclade\tstatus'


def compute_gc(seq: str) -> float:
    """GC 占比(N 计入分母)|GC fraction (N counts in denominator)"""
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq)


def _read_genome(genome_fasta: str) -> Dict[str, str]:
    """读入基因组序列(name→seq)|Load genome sequences (name -> seq)"""
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open(genome_fasta, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def gc_per_class(genome_fasta: str,
                 hits: List[RepeatHit]) -> Dict[str, Dict[str, float]]:
    """每类 repeat 的 bp 与 GC(服务 two-speed AT-rich 重复区识别)
    |Per-class bp and GC (serves two-speed AT-rich region finding)"""
    seqs = _read_genome(genome_fasta)
    acc: Dict[str, List[int]] = {}
    for hit in hits:
        seq = seqs.get(hit.seqid, '')
        seg = seq[hit.start - 1:hit.end].upper()
        slot = acc.setdefault(hit.family_class, [0, 0])  # [bp, gc_count]
        slot[0] += len(seg)
        slot[1] += seg.count('G') + seg.count('C')
    return {cls: {'bp': int(bp), 'gc': (gc / bp if bp else 0.0)}
            for cls, (bp, gc) in acc.items()}


def write_repeat_summary(out_path: str, tbl: Dict, hits: List[RepeatHit],
                         cls_map: Dict[str, Dict[str, str]],
                         lib_fasta: str, genome_fasta: str = '') -> None:
    """写 repeat_summary.tsv(总体行+逐类行)|Write repeat_summary.tsv

    cls_map 为空(TEsorter 降级)时逐类行仍按 .out family_class 统计,
    n_families 取该类 .out 家族与 cls.tsv 归组家族的并集
    |With degraded cls_map, per-class rows still use .out family_class;
    n_families = union of .out families and cls.tsv-grouped families
    """
    total_length = tbl.get('total_length') or 0
    masked_bp = tbl.get('masked_bp')
    if masked_bp is None:
        masked_bp = sum(h.end - h.start + 1 for h in hits)
    masked_pct = tbl.get('masked_pct')
    if masked_pct is None and total_length:
        masked_pct = masked_bp / total_length * 100

    lib_names = set(fasta_ids(lib_fasta))
    per_class_bp: Dict[str, int] = {}
    per_class_fams: Dict[str, set] = {}
    for h in hits:
        per_class_bp[h.family_class] = per_class_bp.get(
            h.family_class, 0) + (h.end - h.start + 1)
        per_class_fams.setdefault(h.family_class, set()).add(h.family)
    # cls.tsv 的 Order/Superfamily 归组(降级时为空)|Grouping from cls.tsv
    classified_by_class: Dict[str, set] = {}
    for family, info in cls_map.items():
        key = f"{info['order']}/{info['superfamily']}"
        classified_by_class.setdefault(key, set()).add(family)

    cls_gc = gc_per_class(genome_fasta, hits) if genome_fasta else {}

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(SUMMARY_HEADER + '\n')
        fh.write('\t'.join([
            'overall', 'all', str(masked_bp),
            f"{masked_pct:.2f}" if masked_pct is not None else 'NA',
            str(len(lib_names & set(h.family for h in hits))), 'NA']) + '\n')
        for cls in sorted(per_class_bp):
            bp = per_class_bp[cls]
            pct = bp / total_length * 100 if total_length else 0.0
            fams = (per_class_fams.get(cls, set())
                    | classified_by_class.get(cls, set()))
            gc = cls_gc.get(cls, {}).get('gc')
            fh.write('\t'.join([
                'per_class', cls, str(bp), f"{pct:.2f}", str(len(fams)),
                f"{gc:.4f}" if gc is not None else 'NA']) + '\n')


def write_families_classified(out_path: str, lib_fasta: str,
                              cls_map: Dict[str, Dict[str, str]]) -> None:
    """写家族分类表(family→Order/Superfamily/Clade,未分类标 unknown)
    |Write family classification table (unclassified marked unknown)"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(FAMILIES_HEADER + '\n')
        for family in fasta_ids(lib_fasta):
            info = cls_map.get(family)
            if info:
                fh.write('\t'.join([family, info['order'], info['superfamily'],
                                    info['clade'], 'classified']) + '\n')
            else:
                fh.write('\t'.join([family, 'NA', 'NA', 'NA', 'unknown']) + '\n')


EFFECTOR_HEADER = ('effector_id\tseqid\tstart\tend\tlength_bp\trepeat_family'
                   '\trepeat_class\toverlap_bp\teffector_coverage'
                   '\tdistance_to_nearest_repeat')
BATCH_HEADER = ('sample\tstatus\tmasked_pct\tn_families'
                '\tclassified_pct\treason')


def write_effector_overlap(out_path: str, regions: List[Dict],
                           hits: List[RepeatHit]) -> None:
    """写 effector↔repeat 交叉检查表(人工核查 edge case 用,不做自动过滤)
    |Write effector↔repeat overlap table (for manual edge-case review)"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    lines = [EFFECTOR_HEADER]
    overlapped = 0
    for region in regions:
        length = int(region['end']) - int(region['start']) + 1
        hits_here = find_overlaps(region, hits)
        hits_here.sort(key=lambda h: overlap_bp(region, h), reverse=True)
        if hits_here:
            overlapped += 1
        for hit in hits_here:
            obp = overlap_bp(region, hit)
            # 退化区(BED start==end 转1-based后 length=0;GFF3 点特征)除零防护,
            # 否则 ZeroDivisionError 逃出样品隔离杀死整批
            # |Guard degenerate regions (length=0) against obp/0, whose
            # ZeroDivisionError would escape sample isolation and kill the batch
            coverage = f"{obp / length:.4f}" if length > 0 else 'NA'
            lines.append('\t'.join([
                str(region['id']), str(region['seqid']),
                str(region['start']), str(region['end']), str(length),
                hit.family, hit.family_class, str(obp),
                coverage, '0']))
        if not hits_here:
            dist = nearest_repeat_distance(region, hits)
            lines.append('\t'.join([
                str(region['id']), str(region['seqid']),
                str(region['start']), str(region['end']), str(length),
                'NA', 'NA', '0', '0.0000',
                str(dist) if dist is not None else 'NA']))
    total = len(regions)
    lines.append('#summary\ttotal={0}\toverlapped={1}\tpct={2:.2f}'.format(
        total, overlapped, (overlapped / total * 100) if total else 0.0))
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')


def write_batch_summary(out_path: str, results: List[Dict]) -> None:
    """写批量模式每样品状态表|Write per-sample batch status table"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(BATCH_HEADER + '\n')
        for r in results:
            fh.write('\t'.join([
                r['sample'], r['status'],
                f"{r['masked_pct']:.2f}" if r.get('masked_pct') is not None else 'NA',
                str(r['n_families']) if r.get('n_families') is not None else 'NA',
                f"{r['classified_pct']:.2f}" if r.get('classified_pct') is not None else 'NA',
                r.get('reason', '')]) + '\n')
