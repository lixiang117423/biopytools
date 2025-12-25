"""
🧬 VCF文件处理模块 | VCF File Processing Module
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from .utils import CommandRunner, get_vcf_stats
except ImportError:
    from utils import CommandRunner, get_vcf_stats


class VCFProcessor:
    """VCF文件处理器 | VCF File Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.temp_dir = None

    def process_vcf(self) -> bool:
        """处理VCF文件 | Process VCF file"""
        self.logger.info("🔬 开始VCF文件处理 | Starting VCF file processing")

        # 获取VCF统计信息 | Get VCF statistics
        stats = get_vcf_stats(self.config.vcf_file, self.logger)
        if stats:
            self.logger.info(f"📊 VCF文件统计 | VCF file statistics:")
            self.logger.info(f"  - 样本数 | Samples: {stats.get('samples', 'unknown')}")
            self.logger.info(
                f"  - 变异数 | Variants: {stats.get('variants', 'unknown')}"
            )
            self.logger.info(
                f"  - 染色体数 | Chromosomes: {len(stats.get('chromosomes', []))}"
            )

        # 步骤1: VCF质量过滤 | Step 1: VCF quality filtering
        if not self._filter_vcf():
            return False

        # 步骤2: 格式转换 | Step 2: Format conversion
        if not self._convert_format():
            return False

        # 步骤3: 样本匹配 | Step 3: Sample matching
        if not self._match_samples():
            return False

        self.logger.info("✅ VCF文件处理完成 | VCF file processing completed")
        return True

    def _filter_vcf(self) -> bool:
        """过滤VCF文件 | Filter VCF file"""
        self.logger.info("🔍 执行VCF质量过滤 | Performing VCF quality filtering")

        # 创建临时目录 | Create temporary directory
        try:
            from .utils import create_temp_dir
        except ImportError:
            from utils import create_temp_dir
        self.temp_dir = create_temp_dir(self.config.working_dir, "vcf_filter")

        # 过滤VCF文件 | Filter VCF file
        filtered_vcf = self.temp_dir / "filtered.vcf.gz"

        cmd = (
            f"{self.config.bcftools_path} filter "
            f"-i 'QUAL>={self.config.qual_min} && DP>={self.config.depth_min} && DP<={self.config.depth_max}' "
            f"{self.config.vcf_file} | "
            f"{self.config.bcftools_path} view "
            f"-i 'F_MISSING<={self.config.missing_threshold}' | "
            f"{self.config.bcftools_path} view "
            f"-i 'MAF>={self.config.maf_threshold}' "
            f"-O z -o {filtered_vcf}"
        )

        success = self.cmd_runner.run(cmd, "VCF质量过滤 | VCF quality filtering")

        if success:
            self.filtered_vcf = filtered_vcf
            self.logger.info(
                f"✅ VCF过滤完成 | VCF filtering completed: {filtered_vcf}"
            )

        return success

    def _convert_format(self) -> bool:
        """转换格式为PLINK | Convert format to PLINK"""
        self.logger.info("🔄 转换VCF格式为PLINK | Converting VCF format to PLINK")

        # 使用PLINK转换格式 | Use PLINK to convert format
        plink_prefix = self.temp_dir / "gwas_data"

        cmd = (
            f"{self.config.plink_path} "
            f"--vcf {self.filtered_vcf} "
            f"--make-bed "
            f"--out {plink_prefix} "
            f"--allow-extra-chr"
        )

        success = self.cmd_runner.run(cmd, "PLINK格式转换 | PLINK format conversion")

        if success:
            self.plink_prefix = plink_prefix
            self.logger.info(
                f"✅ PLINK格式转换完成 | PLINK format conversion completed"
            )

        return success

    def _match_samples(self) -> bool:
        """匹配样本 | Match samples"""
        self.logger.info(
            "🔗 匹配VCF和表型文件样本 | Matching VCF and phenotype file samples"
        )

        try:
            # 读取PLINK样本文件 | Read PLINK sample file
            plink_samples = pd.read_csv(
                f"{self.plink_prefix}.fam",
                sep="\t",
                header=None,
                names=[
                    "family_id",
                    "individual_id",
                    "paternal_id",
                    "maternal_id",
                    "sex",
                    "phenotype",
                ],
            )

            # 读取表型文件 | Read phenotype file
            phenotype_df = pd.read_csv(
                self.config.phenotype_file,
                sep="\t",
                header=None,
                names=["sample_id", "phenotype"],
            )

            # 找到共同样本 | Find common samples
            common_samples = set(plink_samples["individual_id"]) & set(
                phenotype_df["sample_id"]
            )

            if not common_samples:
                self.logger.error(
                    "❌ VCF和表型文件没有共同样本 | No common samples between VCF and phenotype files"
                )
                return False

            self.logger.info(
                f"✅ 找到 {len(common_samples)} 个共同样本 | Found {len(common_samples)} common samples"
            )

            # 过滤PLINK数据 | Filter PLINK data
            plink_samples_filtered = plink_samples[
                plink_samples["individual_id"].isin(common_samples)
            ]
            phenotype_df_filtered = phenotype_df[
                phenotype_df["sample_id"].isin(common_samples)
            ]

            # 保存过滤后的文件 | Save filtered files
            self.plink_prefix_filtered = self.temp_dir / "gwas_data_filtered"
            plink_samples_filtered.to_csv(
                f"{self.plink_prefix_filtered}.fam", sep="\t", header=False, index=False
            )

            # 复制其他PLINK文件 | Copy other PLINK files
            import shutil

            shutil.copy(f"{self.plink_prefix}.bim", f"{self.plink_prefix_filtered}.bim")
            shutil.copy(f"{self.plink_prefix}.bed", f"{self.plink_prefix_filtered}.bed")

            # 保存表型文件 | Save phenotype file
            self.matched_phenotype_file = self.temp_dir / "matched_phenotypes.txt"
            phenotype_df_filtered.to_csv(
                self.matched_phenotype_file, sep="\t", header=False, index=False
            )

            self.logger.info(f"✅ 样本匹配完成 | Sample matching completed")
            return True

        except Exception as e:
            self.logger.error(f"❌ 样本匹配失败 | Sample matching failed: {e}")
            return False

    def get_processed_files(self) -> Dict[str, Path]:
        """获取处理后的文件 | Get processed files"""
        return {
            "plink_prefix": self.plink_prefix_filtered,
            "phenotype_file": self.matched_phenotype_file,
            "temp_dir": self.temp_dir,
        }
