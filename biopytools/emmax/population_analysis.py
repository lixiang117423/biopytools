"""
👥 群体结构分析模块 | Population Structure Analysis Module
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from .utils import CommandRunner
except ImportError:
    from utils import CommandRunner


class PopulationAnalyzer:
    """群体结构分析器 | Population Structure Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.temp_dir = None

    def analyze_population_structure(self, plink_prefix: Path) -> bool:
        """分析群体结构 | Analyze population structure"""
        self.logger.info("🧬 开始群体结构分析 | Starting population structure analysis")

        # 创建临时目录 | Create temporary directory
        try:
            from .utils import create_temp_dir
        except ImportError:
            from utils import create_temp_dir
        self.temp_dir = create_temp_dir(self.config.working_dir, "pop_analysis")

        # 步骤1: LD修剪 | Step 1: LD pruning
        if not self._ld_pruning(plink_prefix):
            return False

        # 步骤2: PCA分析 | Step 2: PCA analysis
        if not self._pca_analysis():
            return False

        # 步骤3: Admixture分析 | Step 3: Admixture analysis
        if not self._admixture_analysis():
            return False

        # 步骤4: 计算亲缘关系矩阵 | Step 4: Calculate kinship matrix
        if not self._calculate_kinship():
            return False

        self.logger.info(
            "✅ 群体结构分析完成 | Population structure analysis completed"
        )
        return True

    def _ld_pruning(self, plink_prefix: Path) -> bool:
        """LD修剪 | LD pruning"""
        self.logger.info("✂️ 执行LD修剪 | Performing LD pruning")

        pruned_prefix = self.temp_dir / "pruned"

        cmd = (
            f"{self.config.plink_path} "
            f"--bfile {plink_prefix} "
            f"--indep-pairwise {self.config.ld_window} {self.config.ld_step} {self.config.ld_r2} "
            f"--make-bed "
            f"--out {pruned_prefix} "
            f"--allow-extra-chr"
        )

        success = self.cmd_runner.run(cmd, "LD修剪 | LD pruning")

        if success:
            self.pruned_prefix = pruned_prefix
            self.logger.info(f"✅ LD修剪完成 | LD pruning completed")

        return success

    def _pca_analysis(self) -> bool:
        """PCA分析 | PCA analysis"""
        self.logger.info("📊 执行PCA分析 | Performing PCA analysis")

        pca_prefix = self.temp_dir / "pca"

        cmd = (
            f"{self.config.plink_path} "
            f"--bfile {self.pruned_prefix} "
            f"--pca {self.config.pca_components} "
            f"--out {pca_prefix} "
            f"--allow-extra-chr"
        )

        success = self.cmd_runner.run(cmd, "PCA分析 | PCA analysis")

        if success:
            self.pca_prefix = pca_prefix
            self.logger.info(f"✅ PCA分析完成 | PCA analysis completed")

        return success

    def _admixture_analysis(self) -> bool:
        """Admixture分析 | Admixture analysis"""
        self.logger.info("🧬 执行Admixture分析 | Performing Admixture analysis")

        k_start, k_end = self.config.admixture_k_range

        for k in range(k_start, min(k_end + 1, 6)):  # 限制最大K值以避免计算时间过长
            self.logger.info(
                f"🔍 执行K={k}的Admixture分析 | Performing Admixture analysis for K={k}"
            )

            admixture_prefix = self.temp_dir / f"admixture_k{k}"

            cmd = f"{self.config.admixture_path} -j4 {self.pruned_prefix}.bed {k}"

            success = self.cmd_runner.run(
                cmd, f"Admixture分析 K={k} | Admixture analysis K={k}"
            )

            if not success:
                self.logger.warning(
                    f"⚠️ K={k}的Admixture分析失败 | Admixture analysis for K={k} failed"
                )
                continue

        self.logger.info(f"✅ Admixture分析完成 | Admixture analysis completed")
        return True

    def _calculate_kinship(self) -> bool:
        """计算亲缘关系矩阵 | Calculate kinship matrix"""
        self.logger.info("🧬 计算亲缘关系矩阵 | Calculating kinship matrix")

        kinship_prefix = self.temp_dir / "kinship"

        cmd = (
            f"{self.config.emmax_path} "
            f"-d -v -h {self.pruned_prefix}.bed "
            f"-p {kinship_prefix}.phen "
            f"-k {kinship_prefix}.kin "
            f"-o {kinship_prefix}"
        )

        success = self.cmd_runner.run(
            cmd, "亲缘关系矩阵计算 | Kinship matrix calculation"
        )

        if success:
            self.kinship_prefix = kinship_prefix
            self.logger.info(
                f"✅ 亲缘关系矩阵计算完成 | Kinship matrix calculation completed"
            )

        return success

    def get_analysis_files(self) -> Dict[str, Path]:
        """获取分析文件 | Get analysis files"""
        return {
            "pruned_prefix": self.pruned_prefix,
            "pca_prefix": self.pca_prefix,
            "kinship_prefix": self.kinship_prefix,
            "temp_dir": self.temp_dir,
        }
