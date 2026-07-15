"""
ps-gene-anno gap 验证报告|Gap evidence report
每个 gap 基因的蛋白证据(miniprot) + RNA-seq 表达(samtools depth) + TE 重叠
|Per-gap: protein evidence + RNA-seq depth + TE overlap → TSV
保持 gap 结果不变, 仅输出报告供用户判断真假|Report-only, no change to gap results
"""

import os
from typing import List, Optional, Dict

from .evidence import MiniprotHit
from .gap_analysis import parse_repeat_out, cds_overlap_ratio


def build_gap_report(gap_hits: List[MiniprotHit], prefix: str,
                     rnaseq_bam: Optional[List[str]], repeat_out: Optional[str],
                     out_tsv: str, config, cmd_runner, logger) -> str:
    """
    生成 gap 验证报告 TSV|Build gap evidence report

    每行一个 gap 基因:坐标 + 蛋白证据 + RNA-seq mean depth + TE 重叠%
    |Per gap: coords + protein evidence + RNA-seq depth + TE overlap

    Args:
        gap_hits: 通过质控的 gap 命中(已建模型)|passed gap hits
        prefix: gap ID 前缀|ID prefix
        rnaseq_bam: RNA-seq BAM 路径列表(可选)|RNA-seq BAMs
        repeat_out: RepeatMasker .out(可选, TE 重叠)|RepeatMasker out
        out_tsv: 报告输出路径|output TSV path
        config/samtools_bin: 提供|provides samtools_bin
        cmd_runner/logger: 命令执行/日志|runner/logger

    Returns:
        报告路径|report path
    """
    logger.info("=" * 60)
    logger.info("gap 验证报告(蛋白+RNA-seq+TE)|gap evidence report")
    logger.info("=" * 60)

    # 1. RNA-seq mean depth(每 gap 的 CDS 区, 复用 braker 的 depth 计算)
    # |RNA-seq depth per gap CDS (reuse braker's compute_region_mean_depth)
    gap_mean_depth: Dict[str, float] = {}
    if rnaseq_bam:
        try:
            from ..braker.repeat_refine import compute_region_mean_depth
            # 每 gap 的 CDS 段作为 region, 映射回 gap_id
            # |each gap CDS segment as region, map back to gap_id
            regions: List = []
            gap_seg_map: List[str] = []   # region_idx → gap_id
            for i, h in enumerate(gap_hits, 1):
                gid = f"{prefix}_gap_{i}"
                for s, e, _ in h.cds_exons:
                    regions.append((h.chrom, s, e))
                    gap_seg_map.append(gid)
            regions_bed = out_tsv + '.regions.bed'
            with open(regions_bed, 'w') as f:
                for chrom, s, e in regions:
                    f.write(f"{chrom}\t{s - 1}\t{e}\n")
            depth = compute_region_mean_depth(
                rnaseq_bam, regions, regions_bed,
                config.samtools_bin, cmd_runner, logger)
            # 汇总 per gap(多 CDS 段取平均)|aggregate per gap
            gap_depth_list: Dict[str, list] = {}
            for idx, d in depth.items():
                if idx < len(gap_seg_map):
                    gap_depth_list.setdefault(gap_seg_map[idx], []).append(d)
            gap_mean_depth = {gid: sum(ds) / len(ds)
                              for gid, ds in gap_depth_list.items()}
            logger.info(f"RNA-seq depth 计算|depth computed: "
                        f"{len(gap_mean_depth)}/{len(gap_hits)} gaps")
            # 清理临时 bed
            if os.path.exists(regions_bed):
                os.remove(regions_bed)
        except Exception as e:
            logger.warning(f"RNA-seq depth 计算失败, 该列置 0|depth failed: {e}")
    else:
        logger.info("未提供 RNA-seq BAM, depth 列置 0|no RNA-seq BAM, depth=0")

    # 2. TE 区(可选, 含 family)|TE regions with family
    te_regions = parse_repeat_out(repeat_out) if repeat_out else {}
    if not repeat_out:
        logger.info("未提供 repeat_out, TE 列置空|no repeat_out, TE empty")

    # 3. 写报告 TSV(含 te_family)|write report
    with open(out_tsv, 'w') as f:
        f.write("gap_id\tchrom\tstart\tend\tstrand\tcds_bp\tcds_aa\t"
                "prot_query\tprot_identity\tprot_coverage\t"
                "rnaseq_mean_depth\tte_overlap_pct\tte_family\n")
        for i, h in enumerate(gap_hits, 1):
            gid = f"{prefix}_gap_{i}"
            cds_bp = sum(e - s + 1 for s, e, _ in h.cds_exons)
            cds_list = [(s, e) for s, e, _ in h.cds_exons]
            te_list = te_regions.get(h.chrom, [])
            te_intervals = [(s, e) for s, e, *_ in te_list]
            te_ov = cds_overlap_ratio(cds_list, te_intervals) if te_intervals else 0.0
            # 最大重叠的 TE family(gap 主要落在哪个 TE 家族)
            # |top-overlap TE family
            best_family = ''
            best_ov_len = 0
            for s, e, fam in te_list:
                ov_len = sum(max(0, min(e, ge) - max(s, gs) + 1)
                             for gs, ge in cds_list)
                if ov_len > best_ov_len:
                    best_ov_len = ov_len
                    best_family = fam
            f.write(
                f"{gid}\t{h.chrom}\t{h.start}\t{h.end}\t{h.strand}\t"
                f"{cds_bp}\t{cds_bp // 3}\t{h.query_id}\t"
                f"{h.identity:.1f}\t{h.coverage:.1f}\t"
                f"{gap_mean_depth.get(gid, 0):.1f}\t{te_ov:.1f}\t{best_family}\n")

    logger.info(f"报告写出|report written: {out_tsv} ({len(gap_hits)} gaps)")
    return out_tsv
