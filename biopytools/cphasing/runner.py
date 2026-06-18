"""
CPhasing命令执行模块|CPhasing Command Execution Module

负责构建和执行CPhasing命令，支持所有子命令
Responsible for building and executing CPhasing commands, supports all subcommands
"""

import os
from typing import List, Optional
from .config import CPhasingConfig
from .utils import CPhasingLogger, CommandRunner, get_conda_env, get_cphasing_env


class CPhasingRunner:
    """
    CPhasing命令执行器|CPhasing Command Runner

    支持两种模式|Supports two modes:
    1. pipeline模式：完整参数封装|Full parameter wrapping
    2. 通用模式：透传子命令和参数|Pass-through subcommand and arguments
    """

    def __init__(self, config: CPhasingConfig):
        self.config = config

        log_file = str(config.output_path / "99_logs" / "cphasing.log")
        self.logger_manager = CPhasingLogger(log_file)
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger)

    def build_pipeline_command(self) -> List[str]:
        """构建pipeline子命令|Build pipeline subcommand"""
        cmd = ['cphasing', 'pipeline']
        cfg = self.config

        if cfg.fasta:
            cmd.extend(['-f', cfg.fasta])
        if cfg.hic1:
            cmd.extend(['-hic1', cfg.hic1])
        if cfg.hic2:
            cmd.extend(['-hic2', cfg.hic2])

        cmd.extend(['-t', str(cfg.threads)])
        cmd.extend(['-n', cfg.groups])

        if cfg.mode != "phasing":
            cmd.extend(['--mode', cfg.mode])
        if cfg.preset != "precision":
            cmd.extend(['--preset', cfg.preset])
        if cfg.output_dir != "./cphasing_output":
            cmd.extend(['--outdir', cfg.output_dir])
        if cfg.steps:
            cmd.extend(['--steps', cfg.steps])
        if cfg.skip_steps:
            cmd.extend(['--skip-steps', cfg.skip_steps])
        if cfg.hic_aligner != "_chromap":
            cmd.extend(['--hic-aligner', cfg.hic_aligner])
        if cfg.hic_mapper_k is not None:
            cmd.extend(['--hic-mapper-k', str(cfg.hic_mapper_k)])
        if cfg.hic_mapper_w is not None:
            cmd.extend(['--hic-mapper-w', str(cfg.hic_mapper_w)])
        if cfg.mapping_quality != 0:
            cmd.extend(['--mapping-quality', str(cfg.mapping_quality)])
        if cfg.hcr:
            cmd.append('--hcr')
        if cfg.pattern:
            cmd.extend(['--pattern', cfg.pattern])
        if cfg.low_memory:
            cmd.append('--low-memory')

        return cmd

    def build_command(self) -> List[str]:
        """
        构建CPhasing命令|Build CPhasing command

        根据subcommand自动选择构建方式|Auto-select build method based on subcommand
        """
        if self.config.subcommand == "pipeline":
            cmd = self.build_pipeline_command()
        else:
            # 通用模式：cphasing <subcommand>|Generic: cphasing <subcommand>
            cmd = ['cphasing', self.config.subcommand]

        # 透传参数对所有子命令都生效（包括pipeline）
        # |Pass-through args apply to all subcommands (including pipeline)
        cmd.extend(self.config.extra_args)

        return cmd

    def run(self) -> bool:
        """
        执行CPhasing命令|Execute CPhasing command

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        subcommand = self.config.subcommand
        self.logger.info(f"开始CPhasing流程|Starting CPhasing: {subcommand}")
        self.logger.info("=" * 60)

        if subcommand == "pipeline":
            self._log_pipeline_config()
        else:
            self._log_generic_config()

        cmd = self.build_command()

        # conda环境包装|conda env wrapper
        conda_env = get_conda_env('cphasing')
        if conda_env:
            self.logger.info(f"conda环境|Conda env: {conda_env}")
            cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output'] + cmd
        else:
            self.logger.warning("未检测到cphasing conda环境|CPhasing conda env not detected")

        extra_env = get_cphasing_env()

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        success, stdout, stderr = self.cmd_runner.run_command(
            cmd,
            description=f"CPhasing {subcommand}",
            extra_env=extra_env,
        )

        if success:
            self.logger.info("=" * 60)
            self.logger.info(f"CPhasing {subcommand} 完成|CPhasing {subcommand} completed successfully")
            self.logger.info("=" * 60)
        else:
            self.logger.error("=" * 60)
            self.logger.error(f"CPhasing {subcommand} 失败|CPhasing {subcommand} failed")
            self.logger.error("=" * 60)

        return success

    def _log_pipeline_config(self):
        """记录pipeline配置|Log pipeline configuration"""
        cfg = self.config
        self.logger.info(f"基因组|Genome: {cfg.fasta}")
        self.logger.info(f"Hi-C R1: {cfg.hic1}")
        self.logger.info(f"Hi-C R2: {cfg.hic2}")
        self.logger.info(f"分组数|Groups: {cfg.groups}")
        self.logger.info(f"线程数|Threads: {cfg.threads}")
        self.logger.info(f"模式|Mode: {cfg.mode}")
        self.logger.info(f"预设|Preset: {cfg.preset}")
        self.logger.info(f"输出目录|Output: {cfg.output_dir}")
        self.logger.info("-" * 60)

    def _log_generic_config(self):
        """记录通用配置|Log generic configuration"""
        self.logger.info(f"子命令|Subcommand: {self.config.subcommand}")
        if self.config.extra_args:
            self.logger.info(f"额外参数|Extra args: {' '.join(self.config.extra_args)}")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info("-" * 60)


# 保持向后兼容|Keep backward compatibility
CPhasingPipeline = CPhasingRunner
