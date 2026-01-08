"""
FASTP处理核心模块|FASTP Processing Core Module
"""

import os
from pathlib import Path
from .utils import CommandRunner


class FastpCore:
    """FASTP核心处理器|FASTP Core Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """
        初始化FASTP核心处理器|Initialize FASTP core processor

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def create_output_directories(self):
        """创建输出目录|Create output directories"""
        self.config.output_path.mkdir(parents=True, exist_ok=True)
        self.config.report_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"创建清洁数据输出目录|Created clean data output directory: {self.config.output_path}")
        self.logger.info(f"创建报告文件输出目录|Created report output directory: {self.config.report_path}")

    def validate_fastp(self) -> bool:
        """
        验证fastp可执行性|Validate fastp executable

        Returns:
            验证是否成功|Whether validation succeeded
        """
        if not self.cmd_runner.check_executable(self.config.fastp_path):
            self.logger.error(
                f"Fastp未找到或无执行权限|"
                f"Fastp not found or no execute permission: {self.config.fastp_path}"
            )
            return False

        self.logger.info(f"Fastp可执行文件验证成功|Fastp executable validated: {self.config.fastp_path}")
        return True

    def process_sample(self, sample_name: str, read1_file: Path, read2_file: Path = None) -> bool:
        """
        处理单个样本|Process single sample

        Args:
            sample_name: 样本名称|Sample name
            read1_file: Read1文件路径|Read1 file path
            read2_file: Read2文件路径（单末端模式为None）|Read2 file path (None for single-end mode)

        Returns:
            处理是否成功|Whether processing succeeded
        """
        # 单末端模式|Single-end mode
        if self.config.single_end:
            # 构建输出文件路径|Construct output file paths
            output_read1 = self.config.output_path / f"{sample_name}.clean.fq.gz"
            html_report = self.config.report_path / f"{sample_name}.html"
            json_report = self.config.report_path / f"{sample_name}.json"

            # 检查输出文件是否已存在|Check if output files already exist
            if output_read1.exists():
                if not self.config.force:
                    self.logger.info(f"样本 {sample_name} 的输出文件已存在，跳过处理|Output files exist for sample {sample_name}, skipping")
                    return True
                else:
                    self.logger.info(f"强制模式：覆盖已存在的文件|Force mode: overwriting existing files for sample {sample_name}")

            # 打印处理信息|Print processing information
            self.logger.info("---")
            self.logger.info(f"正在处理样本（单末端模式）|Processing sample (single-end mode): {sample_name}")
            self.logger.info(f"  Input (输入): {read1_file}")
            self.logger.info(f"  Output (输出): {output_read1}")
            self.logger.info(f"  HTML report (报告): {html_report}")
            self.logger.info("---")

            # 构建fastp命令（单末端）|Build fastp command (single-end)
            cmd = [
                self.config.fastp_path,
                "-i", str(read1_file),
                "-o", str(output_read1),
                "-h", str(html_report),
                "-j", str(json_report),
                "-w", str(self.config.threads),
                "-q", str(self.config.quality_threshold),
                "-l", str(self.config.min_length),
                "-u", str(self.config.unqualified_percent),
                "-n", str(self.config.n_base_limit)
            ]

        else:
            # 双末端模式|Paired-end mode
            # 构建输出文件路径|Construct output file paths
            output_read1 = self.config.output_path / f"{sample_name}_1.clean.fq.gz"
            output_read2 = self.config.output_path / f"{sample_name}_2.clean.fq.gz"
            html_report = self.config.report_path / f"{sample_name}.html"
            json_report = self.config.report_path / f"{sample_name}.json"

            # 检查输出文件是否已存在|Check if output files already exist
            if output_read1.exists() and output_read2.exists():
                if not self.config.force:
                    self.logger.info(f"样本 {sample_name} 的输出文件已存在，跳过处理|Output files exist for sample {sample_name}, skipping")
                    return True
                else:
                    self.logger.info(f"强制模式：覆盖已存在的文件|Force mode: overwriting existing files for sample {sample_name}")

            # 打印处理信息|Print processing information
            self.logger.info("---")
            self.logger.info(f"正在处理样本（双末端模式）|Processing sample (paired-end mode): {sample_name}")
            self.logger.info(f"  Read 1 (输入|input): {read1_file}")
            self.logger.info(f"  Read 2 (输入|input): {read2_file}")
            self.logger.info(f"  Read 1 (输出|output): {output_read1}")
            self.logger.info(f"  Read 2 (输出|output): {output_read2}")
            self.logger.info(f"  HTML report (报告): {html_report}")
            self.logger.info("---")

            # 构建fastp命令|Build fastp command
            cmd = [
                self.config.fastp_path,
                "-i", str(read1_file),
                "-I", str(read2_file),
                "-o", str(output_read1),
                "-O", str(output_read2),
                "-h", str(html_report),
                "-j", str(json_report),
                "-w", str(self.config.threads),
                "-q", str(self.config.quality_threshold),
                "-l", str(self.config.min_length),
                "-u", str(self.config.unqualified_percent),
                "-n", str(self.config.n_base_limit)
            ]

        # 检查是否为dry-run模式|Check if dry-run mode
        if self.config.dry_run:
            self.logger.info(f"[DRY-RUN] 将执行命令|Would execute command for: {sample_name}")
            return True

        # 执行fastp命令|Execute fastp command
        success = self.cmd_runner.run(cmd, f"FASTP质控处理|FASTP quality control -> {sample_name}")

        if success:
            self.logger.info(f"样本 {sample_name} 处理完成|Sample {sample_name} processing completed")
        else:
            self.logger.error(f"样本 {sample_name} 处理失败|Sample {sample_name} processing failed")

        return success
