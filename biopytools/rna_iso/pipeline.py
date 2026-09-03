"""全长转录本分析流程编排|Full-length RNA pipeline orchestration

分派逻辑|Dispatch:
- pacbio_subreads: ccs → refine → fasta → 引擎
- pacbio_ccs:      refine → fasta → 引擎
- reads(fq/fa):    直接引擎(isoquant; isoseq3不可用)
"""

import os
import shutil
from pathlib import Path

from .steps_ccs import run_ccs, run_refine, flnc_to_fasta
from .steps_isoquant import run_isoquant
from .steps_isoseq3 import run_cluster2
from .utils import generate_software_versions_yml


class RnaIsoPipeline:
    """流程编排器|Pipeline orchestrator"""

    def __init__(self, config):
        self.config = config

    def _create_dirs(self, logger):
        """条件创建目录(缺省步不建)|Conditionally create directories"""
        dirs = [self.config.output_dir, self.config.stat_dir,
                self.config.log_dir, self.config.tmp_dir]
        if self.config.needs_ccs:
            dirs.append(self.config.ccs_dir)
        if self.config.needs_refine:
            dirs.append(self.config.refine_dir)
        if self.config.run_isoquant:
            dirs.append(self.config.isoquant_dir)
        if self.config.run_isoseq3:
            dirs.append(self.config.isoseq3_dir)
        for d in dirs:
            Path(d).mkdir(parents=True, exist_ok=True)
        logger.info(f"输出目录|Output directory: {self.config.output_dir}")

    def _cleanup_tmp(self, logger):
        """清理临时目录(§12.4)|Clean tmp directory"""
        if os.path.isdir(self.config.tmp_dir):
            shutil.rmtree(self.config.tmp_dir, ignore_errors=True)
            logger.info("临时文件已清理|Tmp files cleaned")

    def run(self, logger) -> bool:
        """执行全流程|Run the full pipeline"""
        config = self.config
        self._create_dirs(logger)

        generate_software_versions_yml(
            config, os.path.join(config.stat_dir, "software_versions.yml"), logger)

        logger.info("=" * 60)
        logger.info("开始全长转录本分析|Starting full-length RNA analysis")
        logger.info(f"  输入形态|Input kind: {config.input_kind}")
        logger.info(f"  引擎|Engine: {config.engine}")
        logger.info("=" * 60)

        # PacBio前处理链|PacBio preprocessing chain
        flnc_bams = None
        if config.needs_refine:
            if config.needs_ccs:
                ccs_bams = run_ccs(config, logger)
                if not ccs_bams:
                    logger.error("ccs步骤失败,流程终止|ccs step failed, aborting")
                    return False
            else:
                # ccs.bam直接作为refine输入|ccs.bam feeds refine directly
                ccs_bams = list(config.reads)

            flnc_bams = run_refine(config, logger, ccs_bams)
            if not flnc_bams:
                logger.error("refine步骤失败,流程终止|refine step failed, aborting")
                return False

        # IsoQuant引擎输入:PacBio链用FLNC fasta,reads输入原样直传
        # |IsoQuant input: FLNC fasta for PacBio chain, raw reads otherwise
        if config.run_isoquant:
            if flnc_bams is not None:
                flnc_fas = flnc_to_fasta(config, logger, flnc_bams)
                if not flnc_fas:
                    logger.error("FLNC转fasta失败,流程终止|FLNC to fasta failed, aborting")
                    return False
                reads_inputs = flnc_fas
            else:
                reads_inputs = config.reads
            if not run_isoquant(config, logger, reads_inputs):
                logger.error("IsoQuant步骤失败,流程终止|IsoQuant step failed, aborting")
                return False

        # isoseq3引擎(cluster2)|isoseq3 engine (cluster2)
        if config.run_isoseq3:
            if not run_cluster2(config, logger, flnc_bams):
                logger.error("cluster2步骤失败,流程终止|cluster2 step failed, aborting")
                return False

        self._cleanup_tmp(logger)
        logger.info("全长转录本分析完成|Full-length RNA analysis completed")
        return True
