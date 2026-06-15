"""
CPhasing流程执行模块|CPhasing Pipeline Execution Module

负责构建和执行CPhasing pipeline命令
Responsible for building and executing CPhasing pipeline commands
"""

import os
from typing import List
from .config import CPhasingConfig
from .utils import CPhasingLogger, CommandRunner, get_conda_env, get_cphasing_env


class CPhasingPipeline:
    """
    CPhasing流程执行器|CPhasing Pipeline Executor

    封装CPhasing pipeline命令，支持Hi-C数据分相和挂载
    Wraps CPhasing pipeline command for phasing and scaffolding with Hi-C data
    """

    def __init__(self, config: CPhasingConfig):
        self.config = config

        log_file = str(config.output_path / "99_logs" / "cphasing_pipeline.log")
        self.logger_manager = CPhasingLogger(log_file)
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger)

    def build_command(self) -> List[str]:
        """
        构建CPhasing pipeline命令|Build CPhasing pipeline command

        Returns:
            命令列表|Command list
        """
        cmd = ['cphasing', 'pipeline']

        # 必需参数|Required parameters
        cmd.extend(['-f', self.config.fasta])
        cmd.extend(['-hic1', self.config.hic1])
        cmd.extend(['-hic2', self.config.hic2])

        # 核心参数|Core parameters
        cmd.extend(['-t', str(self.config.threads)])
        cmd.extend(['-n', self.config.groups])

        # 模式参数|Mode parameter
        if self.config.mode != "phasing":
            cmd.extend(['--mode', self.config.mode])

        # 预设参数|Preset parameter
        if self.config.preset != "precision":
            cmd.extend(['--preset', self.config.preset])

        # 输出目录|Output directory
        if self.config.output_dir != "./cphasing_output":
            cmd.extend(['--outdir', self.config.output_dir])

        # 步骤控制|Step control
        if self.config.steps:
            cmd.extend(['--steps', self.config.steps])

        if self.config.skip_steps:
            cmd.extend(['--skip-steps', self.config.skip_steps])

        # Hi-C Mapper参数|Hi-C Mapper parameters
        if self.config.hic_aligner != "_chromap":
            cmd.extend(['--hic-aligner', self.config.hic_aligner])

        if self.config.hic_mapper_k is not None:
            cmd.extend(['--hic-mapper-k', str(self.config.hic_mapper_k)])

        if self.config.hic_mapper_w is not None:
            cmd.extend(['--hic-mapper-w', str(self.config.hic_mapper_w)])

        if self.config.mapping_quality != 0:
            cmd.extend(['--mapping-quality', str(self.config.mapping_quality)])

        # HCR参数|HCR parameters
        if self.config.hcr:
            cmd.append('--hcr')

        if self.config.pattern:
            cmd.extend(['--pattern', self.config.pattern])

        # 高级选项|Advanced options
        if self.config.low_memory:
            cmd.append('--low-memory')

        return cmd

    def run(self) -> bool:
        """
        执行CPhasing pipeline|Execute CPhasing pipeline

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始CPhasing流程|Starting CPhasing pipeline")
        self.logger.info("=" * 60)

        self.logger.info(f"基因组|Genome: {self.config.fasta}")
        self.logger.info(f"Hi-C R1|Hi-C R1: {self.config.hic1}")
        self.logger.info(f"Hi-C R2|Hi-C R2: {self.config.hic2}")
        self.logger.info(f"分组数|Groups: {self.config.groups}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"模式|Mode: {self.config.mode}")
        self.logger.info(f"预设|Preset: {self.config.preset}")
        self.logger.info(f"CPhasing目录|CPhasing dir: {self.config.cphasing_dir}")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info("-" * 60)

        cmd = self.build_command()

        # 获取conda环境名|Get conda environment name
        conda_env = get_conda_env('cphasing')
        if conda_env:
            self.logger.info(f"检测到conda环境|Detected conda env: {conda_env}")
            cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output'] + cmd
        else:
            self.logger.warning("未检测到cphasing conda环境|CPhasing conda env not detected")

        # 获取CPhasing运行环境（含bin目录PATH）|Get CPhasing runtime env (with bin PATH)
        extra_env = get_cphasing_env()

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        success, stdout, stderr = self.cmd_runner.run_command(
            cmd,
            description="CPhasing pipeline分析|CPhasing pipeline analysis",
            extra_env=extra_env,
        )

        if success:
            self.logger.info("=" * 60)
            self.logger.info("CPhasing流程完成|CPhasing pipeline completed successfully")
            self.logger.info("=" * 60)
        else:
            self.logger.error("=" * 60)
            self.logger.error("CPhasing流程失败|CPhasing pipeline failed")
            self.logger.error("=" * 60)

        return success
