"""
GCTB分析运行模块|GCTB Analysis Runner Module

负责运行GCTB的贝叶斯分析
Responsible for running GCTB Bayesian analysis
"""

import os
import logging
from pathlib import Path
from typing import Optional
from .config import GCTBConfig
from .utils import CommandRunner


class GCTBAnalyzer:
    """GCTB分析器|GCTB Analyzer"""

    def __init__(self, config: GCTBConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_individual_analysis(self) -> bool:
        """
        运行个体水平分析|Run individual-level analysis

        使用个体基因型和表型数据进行贝叶斯分析
        Use individual genotype and phenotype data for Bayesian analysis

        Returns:
            bool: 分析是否成功|Whether analysis succeeded
        """
        self.logger.info("开始个体水平贝叶斯分析|Starting individual-level Bayesian analysis")

        output_prefix = str(self.config.analysis_dir / f"bayes{self.config.bayes_type.lower()}")

        # 构建GCTB命令|Build GCTB command
        args = [
            '--bfile', str(self.config.qc_dir / "input_qc"),
            '--pheno', self.config.pheno_file,
            '--bayes', self.config.bayes_type,
            '--out', output_prefix
        ]

        # 添加可选参数|Add optional parameters
        if self.config.seed is not None:
            args.extend(['--seed', str(self.config.seed)])
        if self.config.pi is not None:
            args.extend(['--pi', str(self.config.pi)])
        if self.config.sigma_g is not None:
            args.extend(['--sigma-g', str(self.config.sigma_g)])
        if self.config.rho is not None:
            args.extend(['--rho', self.config.rho])

        cmd = self._build_gctb_command(args)

        success = self.cmd_runner.run(
            cmd,
            f"个体水平{self.config.bayes_type}贝叶斯分析|Individual-level Bayesian analysis (Bayes{self.config.bayes_type})"
        )

        if success:
            self.logger.info(f"分析完成|Analysis completed: {output_prefix}")

        return success

    def run_summary_analysis(self, ld_matrix_path: Optional[str] = None) -> bool:
        """
        运行汇总水平分析|Run summary-level analysis

        使用GWAS汇总统计和LD矩阵进行贝叶斯分析
        Use GWAS summary statistics and LD matrix for Bayesian analysis

        Args:
            ld_matrix_path: LD矩阵路径，如果为None则自动检测|LD matrix path, auto-detect if None

        Returns:
            bool: 分析是否成功|Whether analysis succeeded
        """
        self.logger.info("开始汇总水平贝叶斯分析|Starting summary-level Bayesian analysis")

        # 确定LD矩阵路径|Determine LD matrix path
        if ld_matrix_path is None:
            ld_matrix_path = self._get_ld_matrix_path()

        if ld_matrix_path is None or not os.path.exists(ld_matrix_path):
            self.logger.error(f"LD矩阵不存在|LD matrix not found: {ld_matrix_path}")
            return False

        output_prefix = str(self.config.analysis_dir / f"sbayes{self.config.bayes_type.lower()}")

        # 构建GCTB命令|Build GCTB command
        args = [
            '--gwas-summary', self._prepare_summary_stats(),
            '--ldm', ld_matrix_path,
            '--sbayes', self.config.bayes_type,
            '--out', output_prefix
        ]

        # 添加可选参数|Add optional parameters
        if self.config.seed is not None:
            args.extend(['--seed', str(self.config.seed)])
        if self.config.pi is not None:
            args.extend(['--pi', str(self.config.pi)])
        if self.config.sigma_g is not None:
            args.extend(['--sigma-g', str(self.config.sigma_g)])
        if self.config.rho is not None:
            args.extend(['--rho', self.config.rho])

        cmd = self._build_gctb_command(args)

        success = self.cmd_runner.run(
            cmd,
            f"汇总水平{self.config.bayes_type}贝叶斯分析|Summary-level Bayesian analysis (SBayes{self.config.bayes_type})"
        )

        if success:
            self.logger.info(f"分析完成|Analysis completed: {output_prefix}")

        return success

    def _get_ld_matrix_path(self) -> Optional[str]:
        """
        获取LD矩阵路径|Get LD matrix path

        Returns:
            LD矩阵路径或None|LD matrix path or None
        """
        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        if self.config.ld_matrix_type == "sparse":
            return f"{ld_prefix}.ldm.sparse"
        elif self.config.ld_matrix_type == "block":
            return f"{ld_prefix}.ldm.block"
        elif self.config.ld_matrix_type == "eigen":
            return ld_prefix
        else:
            return None

    def _prepare_summary_stats(self) -> str:
        """
        准备GWAS汇总统计数据|Prepare GWAS summary statistics

        如果用户没有提供汇总统计文件，可以从个体数据生成
        If user didn't provide summary stats, generate from individual data

        Returns:
            汇总统计文件路径|Summary statistics file path
        """
        # TODO: 实现从个体数据生成汇总统计的功能
        # TODO: Implement generating summary stats from individual data

        # 暂时假设用户提供了汇总统计文件
        # Temporarily assume user provided summary stats file
        return self.config.pheno_file  # 这里应该修改为实际的汇总统计文件路径

    def _build_gctb_command(self, args: list) -> list:
        """构建GCTB命令|Build GCTB command"""
        from .utils import build_conda_command
        return build_conda_command(self.config.gctb_path, args)

    def get_output_files(self) -> dict:
        """
        获取分析输出文件列表|Get list of analysis output files

        Returns:
            dict: 文件类型到路径的映射|Mapping from file type to path
        """
        if self.config.analysis_mode == "individual":
            output_prefix = str(self.config.analysis_dir / f"bayes{self.config.bayes_type.lower()}")
        else:
            output_prefix = str(self.config.analysis_dir / f"sbayes{self.config.bayes_type.lower()}")

        files = {
            'prefix': output_prefix,
            'snp_res': f"{output_prefix}.snpRes",
            'par_res': f"{output_prefix}.parRes",
            'log': f"{output_prefix}.log"
        }

        return files
