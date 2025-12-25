"""
🧬 GWAS分析引擎模块 | GWAS Analysis Engine Module
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


class GWASEngine:
    """GWAS分析引擎 | GWAS Analysis Engine"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.temp_dir = None

    def run_gwas_analysis(
        self, plink_prefix: Path, phenotype_file: Path, kinship_prefix: Path
    ) -> bool:
        """运行GWAS分析 | Run GWAS analysis"""
        self.logger.info("🧬 开始GWAS关联分析 | Starting GWAS association analysis")

        # 创建临时目录 | Create temporary directory
        try:
            from .utils import create_temp_dir
        except ImportError:
            from utils import create_temp_dir
        self.temp_dir = create_temp_dir(self.config.working_dir, "gwas")

        # 准备表型文件 | Prepare phenotype file
        phenotype_file_prepared = self._prepare_phenotype_file(phenotype_file)

        # 运行EMMAX分析 | Run EMMAX analysis
        if not self._run_emmax_analysis(
            plink_prefix, phenotype_file_prepared, kinship_prefix
        ):
            return False

        # 处理结果文件 | Process result files
        if not self._process_results():
            return False

        self.logger.info("✅ GWAS分析完成 | GWAS analysis completed")
        return True

    def _prepare_phenotype_file(self, phenotype_file: Path) -> Path:
        """准备表型文件 | Prepare phenotype file"""
        self.logger.info("📝 准备表型文件 | Preparing phenotype file")

        # 读取表型文件 | Read phenotype file
        phenotype_df = pd.read_csv(
            phenotype_file, sep="\t", header=None, names=["sample_id", "phenotype"]
        )

        # 保存为EMMAX格式 | Save as EMMAX format
        emmax_phenotype_file = self.temp_dir / "phenotype.txt"
        phenotype_df["phenotype"].to_csv(
            emmax_phenotype_file, header=False, index=False
        )

        self.logger.info(
            f"✅ 表型文件准备完成 | Phenotype file prepared: {emmax_phenotype_file}"
        )
        return emmax_phenotype_file

    def _run_emmax_analysis(
        self, plink_prefix: Path, phenotype_file: Path, kinship_prefix: Path
    ) -> bool:
        """运行EMMAX分析 | Run EMMAX analysis"""
        self.logger.info("🔬 执行EMMAX关联分析 | Performing EMMAX association analysis")

        emmax_prefix = self.temp_dir / "emmax_results"

        cmd = (
            f"{self.config.emmax_path} "
            f"-d -v -h {plink_prefix}.bed "
            f"-p {phenotype_file} "
            f"-k {kinship_prefix}.kin "
            f"-o {emmax_prefix}"
        )

        success = self.cmd_runner.run(cmd, "EMMAX关联分析 | EMMAX association analysis")

        if success:
            self.emmax_prefix = emmax_prefix
            self.logger.info(f"✅ EMMAX分析完成 | EMMAX analysis completed")

        return success

    def _process_results(self) -> bool:
        """处理结果文件 | Process result files"""
        self.logger.info("📊 处理GWAS结果 | Processing GWAS results")

        try:
            # 读取EMMAX结果 | Read EMMAX results
            results_file = f"{self.emmax_prefix}.ps"
            if not os.path.exists(results_file):
                self.logger.error(
                    f"❌ 结果文件不存在 | Result file does not exist: {results_file}"
                )
                return False

            # 读取BIM文件获取SNP信息 | Read BIM file for SNP information
            bim_file = f"{self.emmax_prefix}.bim"
            if not os.path.exists(bim_file):
                self.logger.error(
                    f"❌ BIM文件不存在 | BIM file does not exist: {bim_file}"
                )
                return False

            # 读取结果 | Read results
            results_df = pd.read_csv(
                results_file,
                sep="\t",
                header=None,
                names=["snp_id", "beta", "se", "p_value"],
            )
            bim_df = pd.read_csv(
                bim_file,
                sep="\t",
                header=None,
                names=["chrom", "snp_id", "pos_cm", "bp", "allele1", "allele2"],
            )

            # 合并结果 | Merge results
            merged_df = pd.merge(bim_df, results_df, on="snp_id")

            # 计算-log10(p值) | Calculate -log10(p-value)
            merged_df["neg_log10_p"] = -np.log10(merged_df["p_value"])

            # 保存结果 | Save results
            self.results_file = self.temp_dir / "gwas_results.txt"
            merged_df.to_csv(self.results_file, sep="\t", index=False)

            self.logger.info(
                f"✅ 结果处理完成 | Results processing completed: {self.results_file}"
            )
            return True

        except Exception as e:
            self.logger.error(f"❌ 结果处理失败 | Results processing failed: {e}")
            return False

    def get_results_file(self) -> Path:
        """获取结果文件 | Get results file"""
        return self.results_file
