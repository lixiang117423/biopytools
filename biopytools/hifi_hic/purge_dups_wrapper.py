"""
Purge_Dups去冗余封装模块|Purge_Dups Deduplication Wrapper Module

在hifi_hic流程中集成Purge_Dups去冗余功能|Integrate Purge_Dups deduplication in hifi_hic pipeline
"""

import os
import glob
from pathlib import Path

from ..common.paths import resolve_legacy_path


class PurgeDupsWrapper:
    """Purge_Dups去冗余封装类|Purge_Dups Deduplication Wrapper Class"""

    def __init__(self, config, logger):
        """
        初始化封装器|Initialize wrapper

        Args:
            config: AssemblyConfig对象|AssemblyConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 动态导入purge_dups模块(该模块已conda包装)|Import purge_dups module (already conda-wrapped)
        try:
            from biopytools.purge_dups.main import PurgeDupsRunner
            from biopytools.purge_dups.config import PurgeDupsConfig
            self.PurgeDupsRunner = PurgeDupsRunner
            self.PurgeDupsConfig = PurgeDupsConfig
            self.logger.info("Purge_Dups模块导入成功|Purge_Dups module imported successfully")
        except ImportError as e:
            self.logger.error(f"无法导入Purge_Dups模块|Failed to import Purge_Dups module: {e}")
            raise

    def _find_purged_output(self, input_fa: str, output_subdir: str = None) -> str:
        """
        查找去冗余输出文件|Find purged output file

        Purge_Dups输出目录可能是sequences或seqs,文件名格式多样,逐一尝试+glob兜底
        |Purge_Dups output dir may be sequences/seqs; filename varies; try each + glob fallback

        Args:
            input_fa: 输入基因组FASTA文件|Input genome FASTA file
            output_subdir: 单倍型等子目录|Subdirectory for haplotypes etc.

        Returns:
            找到的文件路径或None|Found file path or None
        """
        # hap{i}/子目录隔离输出,固定名中间文件(cutoffs/dups.bed)不互踩
        # |hap{i}/ subdirs isolate outputs; fixed-name intermediates (cutoffs/dups.bed) never collide
        purge_base = self.config.purge_dups_dir
        if output_subdir:
            purge_base = os.path.join(purge_base, output_subdir)

        sequences_dir = os.path.join(purge_base, "sequences")
        if not os.path.exists(sequences_dir):
            sequences_dir = os.path.join(purge_base, "seqs")
        if not os.path.exists(sequences_dir):
            return None

        input_stem = Path(input_fa).stem
        possible_filenames = [
            f"{input_stem}_purged.purged.fa",  # get_seqs默认输出|get_seqs default output
            f"{input_stem}_purged.purge.fa",   # 旧版本格式|Legacy format
            f"{Path(input_fa).name}_purged.purged.fa",
            f"{input_stem.replace('_primary', '')}_purged.purged.fa",
        ]
        for filename in possible_filenames:
            candidate = os.path.join(sequences_dir, filename)
            if os.path.exists(candidate):
                return candidate

        # glob兜底:按输入stem过滤,同目录多个purged产物时不拿错
        # |glob fallback filtered by input stem; never grabs the wrong file among multiple outputs
        candidates = sorted(glob.glob(os.path.join(sequences_dir, f"{input_stem}_purged.purged.fa")))
        if not candidates:
            candidates = sorted(glob.glob(os.path.join(sequences_dir, f"{input_stem}_purged.purge.fa")))
        if candidates:
            self.logger.info(f"使用自动检测的去冗余文件|Using auto-detected purged file: {candidates[0]}")
            return candidates[0]
        return None

    def run_purge_dups(self, input_fa: str, output_subdir: str = None,
                       force_hifi_reads: bool = False) -> str:
        """
        运行Purge_Dups去冗余流程(支持断点续传,§10.2)|Run Purge_Dups deduplication (resume support)

        Args:
            input_fa: 输入基因组FASTA文件|Input genome FASTA file
            output_subdir: 输出子目录(单倍型用hap{i})|Output subdirectory (hap{i} for haplotypes)
            force_hifi_reads: 强制使用全量HiFi reads(单倍型用,避免NGS筛选子集带偏覆盖度阈值)
                |Force full HiFi reads (for haplotypes; NGS-filtered subset biases haploid cutoffs)

        Returns:
            str: 去冗余后的FASTA文件路径|Path to purged FASTA file
        """
        input_label = output_subdir or Path(input_fa).stem
        self.logger.info("=" * 80)
        self.logger.info(f"开始Purge_Dups去冗余流程[{input_label}]|Starting Purge_Dups Deduplication Pipeline [{input_label}]")
        self.logger.info("=" * 80)

        # 单倍型输出到独立子目录|Haplotypes write into isolated subdirectories
        purge_output_dir = self.config.purge_dups_dir
        if output_subdir:
            purge_output_dir = os.path.join(purge_output_dir, output_subdir)

        try:
            # 确定输入文件|Determine input file
            input_fa = os.path.abspath(input_fa)
            if not os.path.exists(input_fa):
                self.logger.error(f"输入文件不存在|Input file not found: {input_fa}")
                return None

            # 断点续传:输出已存在且未禁用续传则跳过去冗余(§10.2)
            # |Resume: skip dedup if output exists and resume enabled (§10.2)
            existing_purged = self._find_purged_output(input_fa, output_subdir)
            if existing_purged and os.path.getsize(existing_purged) > 0 and self.config.resume:
                self.logger.info(f"检测到去冗余输出已存在，跳过|Purge_dups output exists, skipping: {existing_purged}")
                return existing_purged

            # 确定输入reads文件|Determine input reads file
            # 单倍型固定全量HiFi reads:NGS筛选子集按primary高覆盖contig筛出,
            # 拿去算单套hap的覆盖度直方图会带偏calcuts阈值
            # |Haplotypes always use full HiFi reads: the NGS-filtered subset is selected
            # by primary high-coverage contigs and would bias calcuts cutoffs on a hap
            if force_hifi_reads:
                input_reads = self.config.hifi_data
                self.logger.info(f"使用全量HiFi reads|Using full HiFi reads: {input_reads}")
            elif self.config.has_ngs:
                # 优先使用NGS polish后的reads，否则使用原始HiFi reads
                # Prefer NGS polished reads, otherwise use original HiFi reads
                filtered_reads = os.path.join(resolve_legacy_path(self.config.ngs_polish_dir, "03_filtered_reads"),
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
                output_dir=purge_output_dir,
                threads=self.config.purge_dups_threads,
                read_type=self.config.purge_dups_read_type
            )

            # 创建并运行Purge_Dups|Create and run Purge_Dups
            purge_runner = self.PurgeDupsRunner(
                input=input_fa,
                reads=input_reads,
                output_dir=purge_output_dir,
                threads=self.config.purge_dups_threads,
                read_type=self.config.purge_dups_read_type,
                purge_dups_path=self.config.purge_dups_path
            )

            # 运行完整的去冗余流程|Run complete deduplication pipeline
            success = purge_runner.run_full_pipeline()

            if not success:
                self.logger.error("Purge_Dups去冗余流程失败|Purge_Dups deduplication pipeline failed")
                return None

            # 查找输出文件|Find output files
            purged_fa = self._find_purged_output(input_fa, output_subdir)

            if not purged_fa:
                self.logger.error(f"去冗余输出文件未找到|Purged output file not found in {purge_output_dir}")
                return None

            self.logger.info("=" * 80)
            self.logger.info(f"Purge_Dups去冗余流程完成[{input_label}]|Purge_Dups Deduplication Pipeline Completed [{input_label}]")
            self.logger.info("=" * 80)
            self.logger.info(f"去冗余基因组|Purged genome: {purged_fa}")

            return purged_fa

        except Exception as e:
            self.logger.error(f"Purge_Dups去冗余流程出错|Error in Purge_Dups deduplication pipeline: {e}", exc_info=True)
            return None
