"""KMERIA完整流程处理器|KMERIA Complete Pipeline Processor"""

import os
from ..utils import CommandRunner
from .count import CountProcessor
from .kctm import KctmProcessor
from .filter import FilterProcessor
from .m2b import M2bProcessor
from .asso import AssoProcessor
from .qc import QCProcessor
from .post_gwas import PostGwasProcessor


class PipelineProcessor:
    """KMERIA完整流程处理器|KMERIA Complete Pipeline Processor"""

    def __init__(self, config, logger, cmd_runner=None):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner or CommandRunner(logger, config.output_dir)

        # 初始化各步骤处理器|Initialize step processors
        self.processors = {}

    def run(self) -> bool:
        """运行完整流程|Run complete pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("开始KMERIA完整分析流程|Starting KMERIA complete analysis pipeline")
        self.logger.info("=" * 60)

        # 确定起始步骤|Determine starting step
        start_step = self.config.step
        steps_to_run = self.config.steps

        # 步骤顺序|Step order
        step_order = ['count', 'kctm', 'filter', 'm2b', 'asso']

        if start_step:
            start_idx = step_order.index(start_step)
            steps_to_run = step_order[start_idx:]
            self.logger.info(f"从步骤开始|Starting from step: {start_step}")

        # 执行各步骤|Execute each step
        step_results = {}

        for step in steps_to_run:
            if step not in step_order:
                self.logger.warning(f"跳过未知步骤|Skipping unknown step: {step}")
                continue

            self.logger.info(f"{'='*60}")
            self.logger.info(f"执行步骤|Executing step: {step.upper()}")
            self.logger.info(f"{'='*60}")

            success = self._run_step(step)
            step_results[step] = success

            if not success:
                self.logger.error(f"步骤{step}失败，流程终止|Step {step} failed, pipeline terminated")
                return False

            self.logger.info(f"步骤{step}完成|Step {step} completed")

        # 执行可选功能|Execute optional features
        if self.config.enable_qc:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("执行质控统计|Executing QC statistics")
            self.logger.info("=" * 60)
            self._run_qc()

        if self.config.genome_file:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("执行Post-GWAS分析|Executing Post-GWAS analysis")
            self.logger.info("=" * 60)
            self._run_post_gwas()

        # 完成|Complete
        self.logger.info(f"\n{'='*60}")
        self.logger.info("KMERIA完整分析流程全部完成|KMERIA complete analysis pipeline finished")
        self.logger.info(f"{'='*60}\n")

        # 保存命令日志|Save commands log
        log_file = os.path.join(self.config.output_dir, 'commands_log.json')
        self.cmd_runner.save_commands_log(log_file)

        return True

    def _run_step(self, step: str) -> bool:
        """运行单个步骤|Run single step"""
        try:
            if step == 'count':
                return self._run_count()
            elif step == 'kctm':
                return self._run_kctm()
            elif step == 'filter':
                return self._run_filter()
            elif step == 'm2b':
                return self._run_m2b()
            elif step == 'asso':
                return self._run_asso()
            else:
                self.logger.warning(f"未知步骤|Unknown step: {step}")
                return False
        except Exception as e:
            self.logger.error(f"步骤{step}执行出错|Step {step} execution error: {e}")
            return False

    def _run_count(self) -> bool:
        """运行k-mer计数|Run k-mer counting"""
        from ..config import CountConfig

        count_config = CountConfig(
            fastq_dir=self.config.fastq_dir,
            samples_file=self.config.samples_file,
            output_dir=str(self.config.dirs['counts']),
            kmer_size=self.config.kmer_size,
            threads=self.config.threads,
            batch_size=self.config.batch_size
        )
        count_config.validate()

        processor = CountProcessor(count_config, self.logger, self.cmd_runner)
        return processor.run()

    def _run_kctm(self) -> bool:
        """运行矩阵构建|Run matrix construction"""
        from ..config import KctmConfig

        kctm_config = KctmConfig(
            input_dir=str(self.config.dirs['counts']),
            output_dir=str(self.config.dirs['matrices']),
            threads=self.config.threads
        )
        kctm_config.validate()

        processor = KctmProcessor(kctm_config, self.logger, self.cmd_runner)
        return processor.run()

    def _run_filter(self) -> bool:
        """运行过滤|Run filtering"""
        from ..config import FilterConfig

        filter_config = FilterConfig(
            input_dir=str(self.config.dirs['matrices']),
            output_dir=str(self.config.dirs['filtered']),
            depth_file=self.config.depth_file,
            max_abund=self.config.max_abund,
            missing_ratio=self.config.missing_ratio,
            ploidy=self.config.ploidy,
            threads=self.config.threads
        )
        filter_config.validate()

        processor = FilterProcessor(filter_config, self.logger, self.cmd_runner)
        return processor.run()

    def _run_m2b(self) -> bool:
        """运行格式转换|Run format conversion"""
        from ..config import M2bConfig

        m2b_config = M2bConfig(
            input_dir=str(self.config.dirs['filtered']),
            output_dir=str(self.config.dirs['bimbam']),
            threads=self.config.threads
        )
        m2b_config.validate()

        processor = M2bProcessor(m2b_config, self.logger, self.cmd_runner)
        return processor.run()

    def _run_asso(self) -> bool:
        """运行关联分析|Run association analysis"""
        from ..config import AssoConfig

        asso_config = AssoConfig(
            input_dir=str(self.config.dirs['bimbam']),
            pheno_file=self.config.pheno_file,
            output_dir=str(self.config.dirs['association']),
            pheno_col=self.config.pheno_col,
            covar_file=self.config.covar_file,
            kinship_file=self.config.kinship_file,
            use_bimbam_tools=True,
            threads=self.config.threads
        )
        asso_config.validate()

        processor = AssoProcessor(asso_config, self.logger, self.cmd_runner)
        return processor.run()

    def _run_qc(self):
        """运行质控统计|Run QC statistics"""
        qc_processor = QCProcessor(self.config, self.logger, self.cmd_runner)
        return qc_processor.run()

    def _run_post_gwas(self):
        """运行Post-GWAS分析|Run Post-GWAS analysis"""
        post_gwas_processor = PostGwasProcessor(self.config, self.logger, self.cmd_runner)
        return post_gwas_processor.run()
