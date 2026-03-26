"""
Purge_Dups去冗余封装模块|Purge_Dups Deduplication Wrapper Module

在hifi_hic流程中集成Purge_Dups去冗余功能|Integrate Purge_Dups deduplication in hifi_hic pipeline
"""

import os
import sys
from pathlib import Path


class PurgeDupsWrapper:
    """Purge_Dups去冗余封装类|Purge_Dups Deduplication Wrapper Class"""

    def __init__(self, config, logger):
        """初始化封装器|Initialize wrapper

        Args:
            config: AssemblyConfig对象|AssemblyConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 动态导入purge_dups模块|Dynamically import purge_dups module
        try:
            # 使用绝对导入从biopytools.purge_dups包导入|Use absolute import from biopytools.purge_dups package
            from biopytools.purge_dups.main import PurgeDupsRunner
            from biopytools.purge_dups.config import PurgeDupsConfig
            self.PurgeDupsRunner = PurgeDupsRunner
            self.PurgeDupsConfig = PurgeDupsConfig
            self.logger.info("Purge_Dups模块导入成功|Purge_Dups module imported successfully")
        except ImportError as e:
            self.logger.error(f"无法导入Purge_Dups模块|Failed to import Purge_Dups module: {e}")
            raise

    def run_purge_dups(self, input_fa: str) -> str:
        """运行Purge_Dups去冗余流程|Run Purge_Dups deduplication pipeline

        Args:
            input_fa: 输入基因组FASTA文件|Input genome FASTA file

        Returns:
            str: 去冗余后的FASTA文件路径|Path to purged FASTA file
        """
        self.logger.info("=" * 80)
        self.logger.info("开始Purge_Dups去冗余流程|Starting Purge_Dups Deduplication Pipeline")
        self.logger.info("=" * 80)

        try:
            # 确定输入文件|Determine input file
            input_fa = os.path.abspath(input_fa)
            if not os.path.exists(input_fa):
                self.logger.error(f"输入文件不存在|Input file not found: {input_fa}")
                return None

            # 确定输入reads文件|Determine input reads file
            # 优先使用NGS polish后的reads，否则使用原始HiFi reads
            # Prefer NGS polished reads, otherwise use original HiFi reads
            if self.config.has_ngs:
                # 检查是否有筛选后的reads|Check if there are filtered reads
                filtered_reads = os.path.join(self.config.ngs_polish_dir, "03.filtered_reads",
                                            f"{self.config.prefix}_high_quality_reads.fq.gz")
                if os.path.exists(filtered_reads):
                    input_reads = filtered_reads
                    self.logger.info(f"使用NGS筛选后的reads|Using NGS filtered reads: {input_reads}")
                else:
                    input_reads = self.config.hifi_data
                    self.logger.info(f"使用原始HiFi reads|Using original HiFi reads: {input_reads}")
            else:
                input_reads = self.config.hifi_data
                self.logger.info(f"使用HiFi reads|Using HiFi reads: {input_reads}")

            # 创建Purge_Dups配置|Create Purge_Dups configuration
            purge_config = self.PurgeDupsConfig(
                input=input_fa,
                reads=input_reads,
                purge_dups_path=self.config.purge_dups_path,
                output_dir=self.config.purge_dups_dir,
                threads=self.config.purge_dups_threads,
                read_type=self.config.purge_dups_read_type
            )

            # 创建并运行Purge_Dups|Create and run Purge_Dups
            purge_runner = self.PurgeDupsRunner(
                input=input_fa,
                reads=input_reads,
                output_dir=self.config.purge_dups_dir,
                threads=self.config.purge_dups_threads,
                read_type=self.config.purge_dups_read_type,
                purge_dups_path=self.config.purge_dups_path
            )

            # 运行完整的去冗余流程|Run complete deduplication pipeline
            success = purge_runner.run_full_pipeline()

            if not success:
                self.logger.error("Purge_Dups去冗余流程失败|Purge_Dups deduplication pipeline failed")
                return None

            # 获取输出文件|Get output files
            # Purge_Dups输出目录可能是sequences或seqs|Purge_Dups output directory could be sequences or seqs
            sequences_dir = os.path.join(self.config.purge_dups_dir, "sequences")
            if not os.path.exists(sequences_dir):
                sequences_dir = os.path.join(self.config.purge_dups_dir, "seqs")

            # 尝试多种可能的文件名格式|Try multiple possible filename formats
            input_stem = Path(input_fa).stem
            possible_filenames = [
                f"{input_stem}_purged.purge.fa",  # 原始格式|Original format
                f"{Path(input_fa).name}_purged.purge.fa",  # 完整文件名|Full filename
                f"{input_stem.replace('.primary', '')}_purged.purge.fa",  # 去掉.primary|Remove .primary
                f"{input_stem.replace('.hap', '')}_purged.purge.fa",  # 去掉.hap|Remove .hap
            ]

            purged_fa = None
            for filename in possible_filenames:
                candidate = os.path.join(sequences_dir, filename)
                if os.path.exists(candidate):
                    purged_fa = candidate
                    break

            # 如果还是找不到，扫描目录查找任何*_purged.purge.fa文件|If still not found, scan directory for any *_purged.purge.fa
            if not purged_fa and os.path.exists(sequences_dir):
                import glob
                candidates = glob.glob(os.path.join(sequences_dir, "*_purged.purge.fa"))
                if candidates:
                    # 使用找到的第一个文件|Use the first found file
                    purged_fa = candidates[0]
                    self.logger.info(f"使用自动检测的去冗余文件|Using auto-detected purged file: {purged_fa}")

            if not purged_fa:
                self.logger.error(f"去冗余输出文件未找到|Purged output file not found in {sequences_dir}")
                self.logger.error(f"  尝试的文件名|Tried filenames: {possible_filenames}")
                return None

            self.logger.info("=" * 80)
            self.logger.info("Purge_Dups去冗余流程完成|Purge_Dups Deduplication Pipeline Completed")
            self.logger.info("=" * 80)
            self.logger.info(f"去冗余基因组|Purged genome: {purged_fa}")

            return purged_fa

        except Exception as e:
            self.logger.error(f"Purge_Dups去冗余流程出错|Error in Purge_Dups deduplication pipeline: {e}")
            import traceback
            traceback.print_exc()
            return None
