"""
annorefine 表达证据(唯一比对 reads)|Expression evidence (unique-mapping reads)
- prepare_unique_bam: BAM 过滤为唯一比对(NH==1 / MAPQ 兜底)|filter to unique reads
- compute_hit_depth_breadth: 每 hit 的平均深度 + 覆盖广度|mean depth + coverage breadth per hit

用唯一比对 reads 的目的: 多比对 reads 多落在 TE/重复区, 会让假基因看起来"有表达"。
|unique reads only: multi-mappers inflate TE/repeat regions spuriously
"""

import bisect
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .evidence import hit_key


def prepare_unique_bam(rnaseq_bam: Optional[List[str]], config,
                       cmd_runner, logger) -> Optional[str]:
    """
    过滤 BAM 为唯一比对 reads|Filter BAM to unique-mapping reads
    优先 samtools view -e '[NH]==1'(samtools≥1.18); 失败回退 -q <min_unique_mapq>
    |Prefer NH==1 (samtools≥1.18); fallback to MAPQ threshold
    缓存到 output_dir/tmp/|cache under output_dir/tmp/

    Returns:
        唯一比对 BAM 路径; 无 BAM 返回 None; unique_reads_only=False 返回原 BAM
        |unique BAM path; None if no BAM; raw BAM if filtering disabled
    """
    if not rnaseq_bam:
        return None
    src = rnaseq_bam[0]
    if not getattr(config, 'unique_reads_only', True):
        logger.info("未启用唯一比对过滤, 用原 BAM|unique filter off, using raw BAM")
        return src
    tmp_dir = os.path.join(config.output_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    unique_bam = os.path.join(tmp_dir, 'rnaseq_unique.bam')
    # 优先 NH==1 | prefer NH==1
    cmd_nh = f"{config.samtools_bin} view -b -e '[NH]==1' {src} -o {unique_bam}"
    if cmd_runner.run_command(cmd_nh, "唯一比对过滤(NH==1)|unique filter (NH==1)"):
        cmd_runner.run_command(
            f"{config.samtools_bin} index {unique_bam}",
            "索引 unique BAM|index unique BAM")
        logger.info(f"唯一比对 BAM(NH==1)|unique BAM (NH==1): {unique_bam}")
        return unique_bam
    # 回退 MAPQ | fallback MAPQ
    logger.warning(
        f"NH 过滤失败(旧版 samtools?), 回退 MAPQ≥{config.min_unique_mapq}"
        "|NH filter failed (old samtools?), fallback to MAPQ")
    cmd_q = (f"{config.samtools_bin} view -b -q {config.min_unique_mapq} "
             f"{src} -o {unique_bam}")
    if cmd_runner.run_command(
            cmd_q,
            f"唯一比对过滤(MAPQ≥{config.min_unique_mapq})|unique filter (MAPQ)"):
        cmd_runner.run_command(
            f"{config.samtools_bin} index {unique_bam}",
            "索引 unique BAM|index unique BAM")
        logger.info(f"唯一比对 BAM(MAPQ)|unique BAM (MAPQ): {unique_bam}")
        return unique_bam
    logger.error("唯一比对过滤失败, 回退原 BAM|unique filter failed, using raw BAM")
    return src


def compute_hit_depth_breadth(hits, bam: Optional[str], config,
                              cmd_runner, logger
                              ) -> Dict[tuple, Tuple[float, float]]:
    """
    算每 hit 的 (mean_depth, breadth%)|Compute mean depth + coverage breadth per hit
    - mean_depth = CDS 各位置深度之和 / CDS 总长
    - breadth    = depth≥1 的位置数 / CDS 总长 × 100
    用 samtools depth -a -b(逐碱基含0); 按 hit_key(hit) 键返回
    |samtools depth -a -b (per-position incl 0); keyed by hit_key(hit)
    """
    if not hits:
        return {}
    if not bam:
        return None     # 有候选但无 BAM → 无法算表达, 返回 None 让 qc_filter 跳过表达过滤|no BAM => skip expr filter

    # 1. 收集每 hit 的 CDS 总长 + BED 行|collect per-hit CDS length + BED rows
    hit_total = {}                       # hit_key -> total CDS length
    bed_rows = []
    for h in hits:
        key = hit_key(h)
        total = sum(e - s + 1 for s, e, _ in h.cds_exons)
        hit_total[key] = total
        for s, e, _ in h.cds_exons:
            bed_rows.append((h.chrom, s, e))
    if not bed_rows:
        return {}

    # 2. 写 BED(0-based half-open)|write BED
    bed = os.path.join(config.output_dir, 'tmp', 'expression_regions.bed')
    os.makedirs(os.path.dirname(bed), exist_ok=True)
    with open(bed, 'w') as f:
        for chrom, s, e in bed_rows:
            f.write(f"{chrom}\t{s - 1}\t{e}\n")
    depth_tsv = bed.rsplit('.', 1)[0] + '_depth.tsv'

    # 3. samtools depth -a -b|run depth
    cmd = f"{config.samtools_bin} depth -a -b {bed} {bam} > {depth_tsv}"
    if not cmd_runner.run_command(cmd, "表达深度计算|expression depth"):
        logger.warning("depth 计算失败, 表达过滤将不生效|depth failed, expression filter disabled")
        for p in (bed, depth_tsv):
            if os.path.exists(p):
                os.remove(p)
        return None     # None(非{})→qc_filter 跳过表达过滤, 避免"全判无表达"静默丢光|None => skip

    # 4. 解析 depth → 每染色体的有序 (pos, depth)|parse to per-chrom sorted positions+depths
    #    samtools depth -a -b 输出 BED 区并集的逐位置(含0), 已按 chrom,pos 排序
    pos_by_chrom = defaultdict(list)
    depth_by_chrom = defaultdict(list)
    with open(depth_tsv) as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < 3:
                continue
            try:
                pos_by_chrom[c[0]].append(int(c[1]))
                depth_by_chrom[c[0]].append(float(c[2]))
            except ValueError:
                continue

    # 清理临时|cleanup tmp
    for p in (bed, depth_tsv):
        if os.path.exists(p):
            os.remove(p)

    # 5. 每 hit 在【自己的 CDS 区间】独立求和|per-hit sum over its OWN CDS exons
    #    关键: 重叠 CDS 的位置对每个覆盖它的 hit 都算(不归属单一 hit),
    #    消除排序依赖/不确定性|overlapping positions count for each covering hit (deterministic)
    result = {}
    for h in hits:
        key = hit_key(h)
        total = hit_total[key]
        if total <= 0:
            result[key] = (0.0, 0.0)
            continue
        positions = pos_by_chrom.get(h.chrom)
        depths = depth_by_chrom.get(h.chrom)
        sum_d = 0.0
        covered = 0
        if positions:
            for s, e, _ in h.cds_exons:
                lo = bisect.bisect_left(positions, s)
                hi = bisect.bisect_right(positions, e)
                for i in range(lo, hi):
                    sum_d += depths[i]
                    if depths[i] >= 1:
                        covered += 1
        result[key] = (sum_d / total, covered / total * 100.0)
    return result
