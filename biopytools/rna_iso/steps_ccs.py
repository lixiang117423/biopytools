"""PacBio前处理执行器|PacBio preprocessing executor

ccs(subreads→HiFi) + refine(去引物/polyA/concatemer→FLNC) + FLNC转fasta
|ccs (subreads→HiFi) + refine (primers/polyA/concatemer removal→FLNC) + FLNC→fasta

isoseq3 26.2 CLI:refine 直接吃 primer FASTA,无 lima/polish 环节
|isoseq3 26.2 CLI: refine takes primer FASTA directly; no lima/polish steps
"""

import os
from typing import List

from .utils import run_command


def _is_done(path: str) -> bool:
    """断点续传marker:文件存在且非空(§10.2)|Resume marker: exists and non-empty"""
    return os.path.exists(path) and os.path.getsize(path) > 0


def _movie_name(bam_path: str, suffix: str) -> str:
    """提取movie名(去.bam后缀链)|Extract movie name (strip .bam suffix chain)"""
    base = os.path.basename(bam_path)
    for cut in (".subreads.bam.pbi", ".subreads.bam", ".ccs.bam", ".flnc.bam", suffix):
        if base.endswith(cut):
            return base[: -len(cut)]
    return base


def run_ccs(config, logger) -> List[str]:
    """逐subreads跑ccs生成HiFi reads BAM|Run ccs per subreads to HiFi BAM

    Returns:
        ccs.bam路径列表|List of ccs.bam paths
    """
    outs = []
    for sub in config.reads:
        if sub.endswith(".pbi"):
            continue
        movie = _movie_name(sub, "")
        out = os.path.join(config.ccs_dir, f"{movie}.ccs.bam")
        if _is_done(out):
            logger.info(f"跳过已完成ccs|Skipping completed ccs: {movie}")
            outs.append(out)
            continue
        # pbi缺失仅警告(pbccs可无pbi运行但速度下降)|pbi missing: warn only
        if not os.path.exists(sub + ".pbi"):
            logger.warning(f"缺少索引文件(降速运行)|Missing index (slower run): {sub}.pbi")
        ok = run_command(
            [config.ccs_path, sub, out,
             "--min-passes", str(config.min_passes),
             "--num-threads", str(config.threads)],
            logger
        )
        if not ok or not _is_done(out):
            logger.error(f"ccs失败或无输出|ccs failed or no output: {movie}")
            return []
        outs.append(out)
    return outs


def run_refine(config, logger, ccs_bams: List[str]) -> List[str]:
    """逐ccs.bam跑去引物得FLNC BAM|Run refine per ccs.bam to FLNC BAMs

    Returns:
        flnc.bam路径列表|List of flnc.bam paths
    """
    primers = config.write_primers_fasta()
    flncs = []
    for ccs_bam in ccs_bams:
        movie = _movie_name(ccs_bam, "")
        out = os.path.join(config.refine_dir, f"{movie}.flnc.bam")
        if _is_done(out):
            logger.info(f"跳过已完成refine|Skipping completed refine: {movie}")
            flncs.append(out)
            continue
        ok = run_command(
            [config.isoseq3_path, "refine", ccs_bam, primers, out,
             "-j", str(config.threads), "--verbose"],
            logger
        )
        if not ok or not _is_done(out):
            logger.error(f"refine失败或无输出|refine failed or no output: {movie}")
            return []
        flncs.append(out)
    return flncs


def flnc_to_fasta(config, logger, flnc_bams: List[str]) -> List[str]:
    """FLNC BAM转fasta(IsoQuant输入)|Convert FLNC BAMs to fasta (IsoQuant input)

    Returns:
        fasta路径列表|List of fasta paths
    """
    outs = []
    for flnc in flnc_bams:
        movie = _movie_name(flnc, "")
        out = os.path.join(config.refine_dir, f"{movie}.flnc.fa")
        if _is_done(out):
            logger.info(f"跳过已完成fasta转换|Skipping completed fasta: {movie}")
            outs.append(out)
            continue
        ok = run_command(
            [config.samtools_path, "fasta", "-@", str(config.threads), flnc, "-o", out],
            logger
        )
        if not ok or not _is_done(out):
            logger.error(f"fasta转换失败|fasta conversion failed: {movie}")
            return []
        outs.append(out)
    return outs
