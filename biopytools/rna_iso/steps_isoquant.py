"""IsoQuant 参考引导转录本重建执行器|IsoQuant reference-guided transcript reconstruction executor

IsoQuant 4.0:长reads比对(minimap2内置)+转录本模型构建+定量,输出GTF/GFF3+表达量表
|IsoQuant 4.0: long-read mapping (built-in minimap2) + transcript model
construction + quantification; outputs GTF/GFF3 and expression tables
"""

import glob
import os
from typing import List

from .utils import run_command


def _isoquant_done(config) -> bool:
    """断点marker:核心产物transcript_models.gtf存在|Resume marker: core GTF exists

    注意:IsoQuant 4.0 不在样本目录产 isoquant.log(README为老版行为),故以
    确定性核心产物为marker(e2e实测2026-09-03)
    """
    core = os.path.join(config.isoquant_dir, config.prefix,
                        f"{config.prefix}.transcript_models.gtf")
    return os.path.exists(core) and os.path.getsize(core) > 0


def run_isoquant(config, logger, reads_inputs: List[str]) -> bool:
    """运行IsoQuant|Run IsoQuant

    Args:
        reads_inputs: 传给--fastq的reads文件|Read files passed to --fastq
    """
    if _isoquant_done(config):
        logger.info(f"跳过已完成IsoQuant|Skipping completed IsoQuant: {config.prefix}")
        return True

    cmd = [
        config.isoquant_path,
        "--reference", config.reference,
        "--fastq"] + reads_inputs + [
        "--data_type", config.isoquant_data_type,
        "-t", str(config.threads),
        "-o", config.isoquant_dir,
        "--prefix", config.prefix,
    ]
    if config.genedb:
        cmd += ["--genedb", config.genedb]

    ok = run_command(cmd, logger)
    if not ok:
        logger.error("IsoQuant执行失败|IsoQuant failed")
        return False

    if not _isoquant_done(config):
        logger.error("IsoQuant完成但未找到预期产物(GTF)|IsoQuant finished but no GTF output found")
        return False

    # 摘要:列出关键产物|Summary: list key outputs
    sample_dir = os.path.join(config.isoquant_dir, config.prefix)
    for pattern in ("*.gtf", "*counts*.tsv", "*tpm*.tsv"):
        for f in sorted(glob.glob(os.path.join(sample_dir, pattern))):
            logger.info(f"产物|Output: {os.path.basename(f)}")
    return True
