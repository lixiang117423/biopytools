"""isoseq3 cluster2 转录本聚类执行器|isoseq3 cluster2 transcript clustering executor

isoseq3 26.2: cluster2(官方推荐,比cluster快) 输入FLNC BAM或FOFN,输出transcripts.bam
|isoseq3 26.2: cluster2 (recommended, faster than cluster) takes FLNC BAM/FOFN
"""

import os
from pathlib import Path
from typing import List

from .utils import run_command


def _is_done(path: str) -> bool:
    return os.path.exists(path) and os.path.getsize(path) > 0


def count_fasta_records(fasta_path: str) -> int:
    """纯Python统计fasta序列数(§11.A登录节点安全)|Count fasta records in pure Python"""
    if not os.path.exists(fasta_path):
        return 0
    n = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def run_cluster2(config, logger, flnc_bams: List[str]) -> bool:
    """聚类FLNC得转录本|Cluster FLNC reads into transcripts

    cluster2输入:单文件直接传,多文件写fofn(原生支持)|Single file passed
    directly, multiple via fofn (native support)
    """
    out_bam = os.path.join(config.isoseq3_dir, f"{config.prefix}.transcripts.bam")
    if _is_done(out_bam):
        logger.info(f"跳过已完成cluster2|Skipping completed cluster2: {out_bam}")
        return True

    if len(flnc_bams) == 1:
        cluster_input = flnc_bams[0]
    else:
        Path(config.tmp_dir).mkdir(parents=True, exist_ok=True)
        cluster_input = os.path.join(config.tmp_dir, "flnc.fofn")
        with open(cluster_input, "w") as f:
            f.write("\n".join(flnc_bams) + "\n")

    ok = run_command(
        [config.isoseq3_path, "cluster2", cluster_input, out_bam,
         "-j", str(config.threads)],
        logger
    )
    if not ok or not _is_done(out_bam):
        logger.error("cluster2失败或无输出|cluster2 failed or no output")
        return False

    # fasta保底:cluster2不产fasta时用samtools转|fasta fallback via samtools if absent
    out_fa = out_bam.replace(".bam", ".fasta")
    if not _is_done(out_fa):
        run_command(
            [config.samtools_path, "fasta", "-@", str(config.threads), out_bam, "-o", out_fa],
            logger
        )

    n_transcripts = count_fasta_records(out_fa)
    logger.info(f"转录本数量|Transcript count: {n_transcripts}")
    return True
