"""
甲基化分析模块 | Methylation Analysis Module
包含甲基化比对、提取、差异分析等核心功能
"""

import glob
import os
from pathlib import Path

from .utils import CommandRunner


class MethylationMapper:
    """甲基化比对器 | Methylation Mapper"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_methylation_mapping(self) -> bool:
        """运行甲基化比对分析 | Run methylation mapping analysis"""
        self.logger.info("步骤4: 开始甲基化比对分析...")

        # 获取clean样品
        clean_samples = []
        if os.path.exists(self.config.clean_dir):
            for filename in os.listdir(self.config.clean_dir):
                if filename.endswith("_1_clean.fq.gz"):
                    sample = filename.replace("_1_clean.fq.gz", "")
                    clean_samples.append(sample)

        if not clean_samples:
            self.logger.error("没有找到clean数据文件")
            return False

        self.logger.info(f"发现clean样品: {clean_samples}")

        bismark_index_dir = os.path.join(self.config.mapping_dir, "bismark_index")

        for sample in clean_samples:
            group = self.config.sample_groups.get(sample, "Unknown")
            sample_output_dir = os.path.join(self.config.mapping_dir, sample)

            # 检查该样品是否已经完成甲基化分析
            cx_report_files = glob.glob(
                os.path.join(sample_output_dir, "*CX_report.txt")
            )
            if cx_report_files:
                self.logger.info(
                    f"🚀 跳过甲基化分析: {sample} (分组: {group}) - 已完成"
                )
                continue

            self.logger.info(f"甲基化分析: {sample} (分组: {group})")

            clean_r1 = os.path.join(self.config.clean_dir, f"{sample}_1_clean.fq.gz")
            clean_r2 = os.path.join(self.config.clean_dir, f"{sample}_2_clean.fq.gz")

            os.makedirs(sample_output_dir, exist_ok=True)

            if os.path.exists(clean_r1) and os.path.exists(clean_r2):
                # Step 4.1: Bismark比对
                if not self._run_bismark_alignment(
                    sample, clean_r1, clean_r2, bismark_index_dir, sample_output_dir
                ):
                    continue

                # Step 4.2: 去重复
                if not self._run_deduplication(sample, sample_output_dir):
                    continue

                # Step 4.3: 甲基化提取
                if not self._run_methylation_extraction(
                    sample, sample_output_dir, bismark_index_dir
                ):
                    continue

                # Step 4.4: 生成报告
                self._generate_reports(sample, sample_output_dir)

                self.logger.info(f"样品分析完成: {sample}")
            else:
                self.logger.warning(f"找不到clean文件 - {sample}")

        self.logger.info("步骤4: 甲基化比对分析完成")
        return True

    def _run_bismark_alignment(
        self, sample, clean_r1, clean_r2, bismark_index_dir, output_dir
    ) -> bool:
        """运行Bismark比对 | Run Bismark alignment"""
        self.logger.info(f"  Bismark比对: {sample}")

        # 设置环境变量确保bismark能找到bowtie2
        bowtie2_path = shutil.which(self.config.bowtie2_path)
        if bowtie2_path:
            bowtie2_dir = os.path.dirname(bowtie2_path)
            os.environ["PATH"] = f"{bowtie2_dir}:{os.environ.get('PATH', '')}"

        cmd = (
            f"{self.config.bismark_path} "
            f"--genome {bismark_index_dir} "
            f"--multicore {self.config.threads // 4} "
            f"-1 {clean_r1} -2 {clean_r2} "
            f"--output_dir {output_dir} "
            f"--temp_dir {output_dir} "
            f"--bowtie2 --non_directional"
        )

        return self.cmd_runner.run(cmd, f"Bismark比对: {sample}")

    def _run_deduplication(self, sample, output_dir) -> bool:
        """运行去重复 | Run deduplication"""
        self.logger.info(f"  去重复: {sample}")

        bam_files = glob.glob(os.path.join(output_dir, "*_pe.bam"))
        if not bam_files:
            self.logger.error(f"找不到BAM文件 - {sample}")
            return False

        bam_file = bam_files[0]
        cmd = (
            f"{self.config.deduplicate_bismark_path} "
            f"--paired --output_dir {output_dir} {bam_file}"
        )

        return self.cmd_runner.run(cmd, f"去重复: {sample}")

    def _run_methylation_extraction(
        self, sample, output_dir, bismark_index_dir
    ) -> bool:
        """运行甲基化提取 | Run methylation extraction"""
        self.logger.info(f"  甲基化提取: {sample}")

        dedup_bam_files = glob.glob(os.path.join(output_dir, "*deduplicated.bam"))
        if not dedup_bam_files:
            self.logger.error(f"找不到去重复BAM文件 - {sample}")
            return False

        dedup_bam = dedup_bam_files[0]
        cmd = (
            f"{self.config.bismark_methylation_extractor_path} "
            f"--paired-end --comprehensive --merge_non_CpG --report "
            f"--cytosine_report --genome_folder {bismark_index_dir} "
            f"--output {output_dir} --multicore {self.config.threads} {dedup_bam}"
        )

        return self.cmd_runner.run(cmd, f"甲基化提取: {sample}")

    def _generate_reports(self, sample, output_dir):
        """生成HTML报告 | Generate HTML reports"""
        self.logger.info(f"  生成HTML报告: {sample}")

        # 生成HTML报告
        try:
            alignment_reports = glob.glob(os.path.join(output_dir, "*_PE_report.txt"))
            dedup_reports = glob.glob(
                os.path.join(output_dir, "*_deduplication_report.txt")
            )
            splitting_reports = glob.glob(
                os.path.join(output_dir, "*_splitting_report.txt")
            )
            mbias_reports = glob.glob(os.path.join(output_dir, "*_M-bias.txt"))

            if alignment_reports:
                cmd = (
                    f"{self.config.bismark2report_path} --dir {output_dir} "
                    f"--alignment_report {alignment_reports[0]}"
                )
                if dedup_reports:
                    cmd += f" --dedup_report {dedup_reports[0]}"
                if splitting_reports:
                    cmd += f" --splitting_report {splitting_reports[0]}"
                if mbias_reports:
                    cmd += f" --mbias_report {mbias_reports[0]}"

                self.cmd_runner.run(cmd, f"生成HTML报告: {sample}")
        except Exception as e:
            self.logger.warning(f"HTML报告生成可能有问题: {e}")

        # 生成汇总报告
        try:
            cmd = f"{self.config.bismark2summary_path} --dir {output_dir}"
            self.cmd_runner.run(cmd, f"生成汇总报告: {sample}")
        except Exception as e:
            self.logger.warning(f"汇总报告生成可能有问题: {e}")
