"""
质量控制模块|Quality Control Module

集成fastp质控功能，支持从原始FASTQ到清洁FASTQ的处理
"""

import os
from pathlib import Path
from .utils import CommandRunner, check_file_exists


class QualityController:
    """质控控制器|Quality Control Controller"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_quality_control(self) -> bool:
        """
        运行质量控制流程|Run quality control pipeline

        Returns:
            bool: 是否成功|Success status
        """
        step_name = "quality_control"

        # 检查是否跳过质控|Check if skipping QC
        if self.config.skip_qc:
            self.logger.info("跳过质控步骤（用户指定）|Skipping QC step (user specified)")
            self.config.clean_fastq_dir = self.config.input_path
            return True

        # 检查断点续传|Check for resume
        if not self.config.force_restart and self._check_qc_completed():
            self.logger.info("质控已完成（文件已存在）|QC completed (files exist)")
            return True

        self.logger.info("=" * 80)
        self.logger.info("步骤1: 质量控制|Step 1: Quality Control")
        self.logger.info("=" * 80)

        # 创建输出目录|Create output directories
        Path(self.config.clean_fastq_dir).mkdir(parents=True, exist_ok=True)

        # 构建fastp命令|Build fastp command
        cmd = self._build_fastp_command()

        # 执行质控|Execute QC
        success = self.cmd_runner.run(cmd, "执行fastp质控|Running fastp QC")

        if success:
            # 验证输出文件|Verify output files
            if self._verify_qc_output():
                self.logger.info("质控完成|Quality control completed successfully")
            else:
                self.logger.error("质控输出验证失败|QC output verification failed")
                return False
        else:
            self.logger.error("质控执行失败|QC execution failed")
            return False

        return True

    def _build_fastp_command(self) -> str:
        """构建fastp命令|Build fastp command"""
        # 基础命令
        cmd_parts = [
            "biopytools", "fastp",
            "-i", self.config.input_path,
            "-o", self.config.clean_fastq_dir,
            "-t", str(self.config.qc_threads)
        ]

        # 添加质控参数|Add QC parameters
        if self.config.qc_quality_threshold != 20:
            cmd_parts.extend(["-q", str(self.config.qc_quality_threshold)])

        if self.config.qc_min_length != 50:
            cmd_parts.extend(["-l", str(self.config.qc_min_length)])

        if self.config.qc_unqualified_percent != 40:
            cmd_parts.extend(["-u", str(self.config.qc_unqualified_percent)])

        if self.config.qc_n_base_limit != 10:
            cmd_parts.extend(["-n", str(self.config.qc_n_base_limit)])

        # 添加文件模式（如果指定）|Add file patterns if specified
        if self.config.qc_read1_suffix:
            cmd_parts.extend(["--read1-suffix", self.config.qc_read1_suffix])

        if self.config.qc_read2_suffix:
            cmd_parts.extend(["--read2-suffix", self.config.qc_read2_suffix])

        # 添加单末端模式（如果指定）|Add single-end mode if specified
        if self.config.qc_single_end:
            cmd_parts.append("--single-end")

        # 添加fastp路径（如果指定）|Add fastp path if specified
        if self.config.fastp_path != "fastp":
            cmd_parts.extend(["--fastp-path", self.config.fastp_path])

        # 添加强制覆盖|Add force overwrite
        if self.config.force_restart:
            cmd_parts.append("--force")

        return " ".join(cmd_parts)

    def _check_qc_completed(self) -> bool:
        """检查质控是否已完成|Check if QC is completed"""
        # 检查输出目录中是否有清洁的FASTQ文件
        clean_dir = Path(self.config.clean_fastq_dir)

        if not clean_dir.exists():
            return False

        # 检查是否有输出文件
        fq_files = list(clean_dir.glob("*.fq.gz")) + list(clean_dir.glob("*.fastq.gz"))

        return len(fq_files) > 0

    def _verify_qc_output(self) -> bool:
        """验证质控输出|Verify QC output"""
        import time

        clean_dir = Path(self.config.clean_fastq_dir)

        self.logger.info(f"验证质控输出目录|Verifying QC output directory: {clean_dir}")
        self.logger.info(f"目录存在|Directory exists: {clean_dir.exists()}")

        if not clean_dir.exists():
            self.logger.error(f"输出目录不存在|Output directory does not exist: {clean_dir}")
            return False

        # 等待一小段时间确保文件系统同步|Wait a bit to ensure filesystem sync
        time.sleep(1)

        # 统计输出文件数量
        clean_files = list(clean_dir.glob("*.fq.gz")) + list(clean_dir.glob("*.fastq.gz"))

        self.logger.info(f"搜索模式|Search pattern: *.fq.gz, *.fastq.gz")
        self.logger.info(f"找到的文件|Files found: {len(clean_files)}")

        if len(clean_files) == 0:
            self.logger.error("未找到质控后的清洁文件|No clean files found after QC")
            self.logger.error(f"请检查输出目录内容|Please check output directory contents: {clean_dir}")
            # 列出目录中的所有文件以便调试
            all_files = list(clean_dir.glob("*"))
            self.logger.error(f"目录中的所有文件|All files in directory:")
            for f in all_files:
                if f.is_file():
                    self.logger.error(f"   - {f.name}")
            return False

        self.logger.info(f"找到 {len(clean_files)} 个清洁FASTQ文件|Found {len(clean_files)} clean FASTQ files")
        for f in clean_files[:5]:  # 只显示前5个文件
            self.logger.info(f"   - {f.name}")

        # 验证文件是否可读且非空|Verify files are readable and non-empty
        for f in clean_files:
            if f.stat().st_size == 0:
                self.logger.error(f"文件为空|File is empty: {f.name}")
                return False

        return True
