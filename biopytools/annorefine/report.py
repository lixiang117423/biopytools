"""
annorefine gap 验证报告|Gap evidence report
每个 gap 基因: 蛋白证据(miniprot) + RNA-seq raw depth + StringTie FPKM/TPM + TE family
|Per-gap: protein + raw depth + FPKM/TPM + TE family → TSV
"""

import os
import re
from typing import List, Optional, Dict, Tuple

from .evidence import MiniprotHit
from .gap_analysis import parse_repeat_out, cds_overlap_ratio


def _parse_stringtie_fpkm(stringtie_gtf: str) -> Dict[str, Tuple[float, float]]:
    """
    解析 StringTie 输出 gtf → {transcript_id: (fpkm, tpm)}
    |Parse StringTie gtf to per-transcript FPKM/TPM
    """
    fpkm_tpm: Dict[str, Tuple[float, float]] = {}
    if not os.path.exists(stringtie_gtf):
        return fpkm_tpm
    with open(stringtie_gtf) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] != 'transcript':
                continue
            attr = cols[8]
            # transcript_id("X" 或 X)|transcript_id
            tid_m = (re.search(r'transcript_id "([^"]+)"', attr)
                     or re.search(r'transcript_id=([^;]+)', attr))
            fpkm_m = (re.search(r'FPKM "([^"]+)"', attr)
                      or re.search(r'FPKM=([^;]+)', attr))
            tpm_m = (re.search(r'TPM "([^"]+)"', attr)
                     or re.search(r'TPM=([^;]+)', attr))
            if tid_m:
                tid = tid_m.group(1)
                fpkm = float(fpkm_m.group(1)) if fpkm_m else 0.0
                tpm = float(tpm_m.group(1)) if tpm_m else 0.0
                fpkm_tpm[tid] = (fpkm, tpm)
    return fpkm_tpm


def build_gap_report(gap_hits: List[MiniprotHit], prefix: str,
                     rnaseq_bam: Optional[List[str]],
                     repeat_out: Optional[str],
                     out_tsv: str,
                     gap_filled_gff3: str,
                     config, cmd_runner, logger,
                     expression: Optional[dict] = None,
                     unique_bam: Optional[str] = None) -> str:
    """
    生成 gap 验证报告 TSV|Build gap evidence report

    每行: 坐标 + 蛋白证据 + depth + coverage_breadth + FPKM/TPM + TE family
    |Per gap: coords + protein + depth + breadth + FPKM/TPM + TE family

    Args:
        expression: {id(hit): (mean_depth, breadth%)} 预计算表达证据, 给定则复用
                    |precomputed expression; if given, reused (no recompute)
        unique_bam: 唯一比对 BAM 路径, StringTie 定量用(优先于 rnaseq_bam)
                    |unique BAM for StringTie quant (priority over rnaseq_bam)
    """
    logger.info("=" * 60)
    logger.info("gap 验证报告(蛋白+RNA-seq depth+breadth+FPKM/TPM+TE)|gap report")
    logger.info("=" * 60)

    # 1. 表达证据(depth+breadth): 优先复用预计算 expression, 否则现算 depth(breadth 置0)
    # |expression: prefer precomputed; else compute depth (breadth=0)
    gap_depth_breadth: Dict[int, Tuple[float, float]] = {}
    if expression is not None:
        gap_depth_breadth = expression
        logger.info(f"复用预计算表达证据|reuse precomputed expression: "
                    f"{len(gap_depth_breadth)}")
    elif rnaseq_bam:
        try:
            from ..braker.repeat_refine import compute_region_mean_depth
            regions: List = []
            seg_map: List[int] = []
            for h in gap_hits:
                for s, e, _ in h.cds_exons:
                    regions.append((h.chrom, s, e))
                    seg_map.append(id(h))
            regions_bed = out_tsv + '.regions.bed'
            with open(regions_bed, 'w') as f:
                for chrom, s, e in regions:
                    f.write(f"{chrom}\t{s - 1}\t{e}\n")
            depth = compute_region_mean_depth(
                rnaseq_bam, regions, regions_bed,
                config.samtools_bin, cmd_runner, logger)
            depth_by_hit: Dict[int, list] = {}
            for idx, d in depth.items():
                if idx < len(seg_map):
                    depth_by_hit.setdefault(seg_map[idx], []).append(d)
            gap_depth_breadth = {hid: (sum(ds) / len(ds), 0.0)
                                 for hid, ds in depth_by_hit.items()}
            logger.info(f"depth 现算|depth computed: "
                        f"{len(gap_depth_breadth)}/{len(gap_hits)} gaps")
            if os.path.exists(regions_bed):
                os.remove(regions_bed)
        except Exception as e:
            logger.warning(f"depth 计算失败|depth failed: {e}")

    # 2. StringTie FPKM/TPM(优先用 unique_bam)|StringTie quant (unique BAM preferred)
    fpkm_tpm: Dict[str, Tuple[float, float]] = {}
    bam_for_quant = unique_bam or (rnaseq_bam[0] if rnaseq_bam else None)
    if bam_for_quant and gap_filled_gff3 and os.path.exists(gap_filled_gff3):
        stringtie_gtf = out_tsv + '.stringtie.gtf'
        cmd = (f"{config.stringtie_bin} {bam_for_quant} -e "
               f"-G {gap_filled_gff3} -o {stringtie_gtf}")
        if cmd_runner.run_command(cmd, "StringTie FPKM/TPM 定量|StringTie quant"):
            fpkm_tpm = _parse_stringtie_fpkm(stringtie_gtf)
            logger.info(f"StringTie 定量|quantified: {len(fpkm_tpm)} transcripts")
        if os.path.exists(stringtie_gtf):
            os.remove(stringtie_gtf)
    else:
        logger.info("无 BAM 或 gap_filled, FPKM/TPM 列置 0"
                    "|no BAM/gap_filled, FPKM/TPM=0")

    # 3. TE 区(含 family)|TE regions with family
    te_regions = parse_repeat_out(repeat_out) if repeat_out else {}

    # 4. 写报告 TSV(含 coverage_breadth)|write report
    with open(out_tsv, 'w') as f:
        f.write("gap_id\tchrom\tstart\tend\tstrand\tcds_bp\tcds_aa\t"
                "prot_query\tprot_identity\tprot_coverage\t"
                "rnaseq_mean_depth\tcoverage_breadth\tfpkm\ttpm\t"
                "te_overlap_pct\tte_family\n")
        for i, h in enumerate(gap_hits, 1):
            gid = f"{prefix}_gap_{i}"
            tid = f"{gid}.t1"
            cds_bp = sum(e - s + 1 for s, e, _ in h.cds_exons)
            cds_list = [(s, e) for s, e, _ in h.cds_exons]
            te_list = te_regions.get(h.chrom, [])
            te_intervals = [(s, e) for s, e, *_ in te_list]
            te_ov = cds_overlap_ratio(cds_list, te_intervals) if te_intervals else 0.0
            # 最大重叠 TE family|top-overlap TE family
            best_family = ''
            best_ov_len = 0
            for s, e, fam in te_list:
                ov_len = sum(max(0, min(e, ge) - max(s, gs) + 1)
                             for gs, ge in cds_list)
                if ov_len > best_ov_len:
                    best_ov_len = ov_len
                    best_family = fam
            depth, breadth = gap_depth_breadth.get(id(h), (0.0, 0.0))
            fpkm, tpm = fpkm_tpm.get(tid, (0.0, 0.0))
            f.write(
                f"{gid}\t{h.chrom}\t{h.start}\t{h.end}\t{h.strand}\t"
                f"{cds_bp}\t{cds_bp // 3}\t{h.query_id}\t"
                f"{h.identity:.1f}\t{h.coverage:.1f}\t"
                f"{depth:.1f}\t{breadth:.1f}\t"
                f"{fpkm:.4f}\t{tpm:.4f}\t"
                f"{te_ov:.1f}\t{best_family}\n")

    logger.info(f"报告写出|report written: {out_tsv} ({len(gap_hits)} gaps)")
    return out_tsv
