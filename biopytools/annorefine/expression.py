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
    unique_bam = os.path.join(tmp_dir, 'rnaseq.unique.bam')
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
                              ) -> Dict[int, Tuple[float, float]]:
    """
    算每 hit 的 (mean_depth, breadth%)|Compute mean depth + coverage breadth per hit
    - mean_depth = CDS 各位置深度之和 / CDS 总长
    - breadth    = depth≥1 的位置数 / CDS 总长 × 100
    用 samtools depth -a -b(逐碱基含0); 按 id(hit) 键返回
    |samtools depth -a -b (per-position incl 0); keyed by id(hit)
    """
    if not hits or not bam:
        return {}

    # 1. 收集每 hit 的 CDS exon + 建 position→hit 查找|collect CDS exons + position->hit map
    exon_by_chrom = defaultdict(list)   # chrom -> [(start, end, hit_key)]
    hit_stats = {}                      # hit_key -> [sum_depth, covered_bases, total_bases]
    bed_rows = []
    for h in hits:
        key = id(h)
        total = 0
        for s, e, _ in h.cds_exons:
            exon_by_chrom[h.chrom].append((s, e, key))
            bed_rows.append((h.chrom, s, e))
            total += e - s + 1
        hit_stats[key] = [0.0, 0, total]
    if not bed_rows:
        return {}

    # 2. 写 BED(0-based half-open)|write BED
    bed = os.path.join(config.output_dir, 'tmp', 'expression_regions.bed')
    os.makedirs(os.path.dirname(bed), exist_ok=True)
    with open(bed, 'w') as f:
        for chrom, s, e in bed_rows:
            f.write(f"{chrom}\t{s - 1}\t{e}\n")
    depth_tsv = bed + '.depth.tsv'

    # 3. samtools depth -a -b|run depth
    cmd = f"{config.samtools_bin} depth -a -b {bed} {bam} > {depth_tsv}"
    if not cmd_runner.run_command(cmd, "表达深度计算|expression depth"):
        logger.warning("depth 计算失败|depth computation failed")
        for p in (bed, depth_tsv):
            if os.path.exists(p):
                os.remove(p)
        return {}

    # 排序 exon 便于 bisect 查找|sort exons per chrom for bisect
    sorted_exons = {}
    for chrom, lst in exon_by_chrom.items():
        lst.sort()
        sorted_exons[chrom] = (lst, [x[0] for x in lst])

    # 4. 解析 depth, 累加到 hit|parse depth, accumulate per hit
    with open(depth_tsv) as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < 3:
                continue
            try:
                pos = int(c[1])
                d = float(c[2])
            except ValueError:
                continue
            info = sorted_exons.get(c[0])
            if not info:
                continue
            lst, starts = info
            idx = bisect.bisect_right(starts, pos) - 1   # start≤pos 的最后一个|last exon start<=pos
            if idx < 0:
                continue
            s, e, key = lst[idx]
            if s <= pos <= e:
                st = hit_stats[key]
                st[0] += d
                if d >= 1:
                    st[1] += 1

    # 清理临时|cleanup tmp
    for p in (bed, depth_tsv):
        if os.path.exists(p):
            os.remove(p)

    # 5. 汇总|summarize
    result = {}
    for h in hits:
        s_depth, covered, total = hit_stats[id(h)]
        if total > 0:
            result[id(h)] = (s_depth / total, covered / total * 100.0)
        else:
            result[id(h)] = (0.0, 0.0)
    return result
