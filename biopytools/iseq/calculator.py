"""
iSeq下载核心计算逻辑|iSeq Download Core Calculation Logic
"""

import os
from pathlib import Path
from typing import Dict, List, Optional


class ISeqCalculator:
    """iSeq下载计算器|iSeq Download Calculator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 输出文件路径|Output file paths
        self.metadata_file = os.path.join(self.config.output_dir, "metadata.txt")
        self.download_log = os.path.join(self.config.output_dir, "download_summary.txt")

    def build_command(self) -> str:
        """构建iSeq命令|Build iSeq command"""
        cmd_parts = [self.config.iseq_path]

        # 必需参数|Required parameter
        cmd_parts.extend(["-i", self.config.accession])

        # 输出目录|Output directory
        cmd_parts.extend(["-o", self.config.output_dir])

        # 元数据选项|Metadata option
        if self.config.metadata_only:
            cmd_parts.append("-m")

        # 下载格式选项|Download format option
        if self.config.gzip:
            cmd_parts.append("-g")
        elif self.config.fastq:
            cmd_parts.append("-q")

        # 合并选项|Merge option
        if self.config.merge:
            cmd_parts.extend(["-e", self.config.merge])

        # 性能参数|Performance parameters
        cmd_parts.extend(["-t", str(self.config.threads)])
        cmd_parts.extend(["-p", str(self.config.parallel)])

        # 速度限制|Speed limit
        if self.config.speed:
            cmd_parts.extend(["-s", str(self.config.speed)])

        # 数据库选项|Database option
        if self.config.database != 'ena':
            cmd_parts.extend(["-d", self.config.database])

        # Aspera选项|Aspera option
        if self.config.use_aspera:
            cmd_parts.append("-a")

        # 跳过MD5|Skip MD5
        if self.config.skip_md5:
            cmd_parts.append("-k")

        # 协议选项|Protocol option
        if self.config.protocol != 'ftp':
            cmd_parts.extend(["-r", self.config.protocol])

        # 静默模式|Quiet mode
        if self.config.quiet:
            cmd_parts.append("-Q")

        return " ".join(cmd_parts)

    def validate_accession(self) -> bool:
        """验证accession格式是否有效|Validate if accession format is valid"""
        acc = self.config.accession.upper()

        # 支持的accession前缀|Supported accession prefixes
        valid_prefixes = {
            'Projects': ['PRJEB', 'PRJNA', 'PRJDB', 'PRJC', 'GSE'],
            'Studies': ['ERP', 'DRP', 'SRP', 'CRA'],
            'BioSamples': ['SAMD', 'SAME', 'SAMN', 'SAMC'],
            'Samples': ['ERS', 'DRS', 'SRS', 'GSM'],
            'Experiments': ['ERX', 'DRX', 'SRX', 'CRX'],
            'Runs': ['ERR', 'DRR', 'SRR', 'CRR']
        }

        all_prefixes = []
        for prefixes in valid_prefixes.values():
            all_prefixes.extend(prefixes)

        for prefix in all_prefixes:
            if acc.startswith(prefix):
                self.logger.info(f"有效accession格式|Valid accession format: {acc} (类型|Type: {prefix})")
                return True

        self.logger.warning(
            f"无法验证accession格式|Cannot verify accession format: {acc}. "
            f"将尝试下载|Will attempt download anyway."
        )
        return True

    def run_download(self) -> bool:
        """运行下载|Run download"""
        self.logger.info("开始iSeq下载|Starting iSeq download")
        self.logger.info(f"Accession|Accession: {self.config.accession}")

        # 验证accession|Validate accession
        if not self.validate_accession():
            self.logger.error("无效的accession格式|Invalid accession format")
            return False

        # 构建命令|Build command
        cmd = self.build_command()
        self.logger.info(f"构建iSeq命令|iSeq command built successfully")

        # 检查是否仅获取元数据|Check if metadata only
        if self.config.metadata_only:
            self.logger.info("仅获取元数据模式|Metadata only mode")
        else:
            mode_desc = []
            if self.config.gzip:
                mode_desc.append("下载gzip格式|Download gzip format")
            elif self.config.fastq:
                mode_desc.append("转换为FASTQ|Convert to FASTQ")

            if self.config.use_aspera:
                mode_desc.append("使用Aspera高速下载|Use Aspera high-speed download")

            self.logger.info(f"下载模式|Download mode: {', '.join(mode_desc)}")

        # 执行下载|Execute download
        success = self.cmd_runner.run(cmd, "iSeq下载|iSeq download")

        if success:
            self.logger.info("iSeq下载完成|iSeq download completed successfully")
            self._write_summary()
        else:
            self.logger.error("iSeq下载失败|iSeq download failed")
            return False

        return success

    def _write_summary(self):
        """写入下载汇总信息|Write download summary"""
        self.logger.info(f"写入下载汇总|Writing download summary: {self.download_log}")

        try:
            with open(self.download_log, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("iSeq下载汇总报告|iSeq Download Summary Report\n")
                f.write("=" * 80 + "\n\n")

                # 基本信息|Basic information
                f.write("1. 基本信息|Basic Information\n")
                f.write("-" * 80 + "\n")
                f.write(f"Accession|Accession: {self.config.accession}\n")
                f.write(f"输出目录|Output directory: {self.config.output_dir}\n")
                f.write(f"仅元数据|Metadata only: {self.config.metadata_only}\n\n")

                # 下载选项|Download options
                f.write("2. 下载选项|Download Options\n")
                f.write("-" * 80 + "\n")
                f.write(f"Gzip格式|Gzip format: {self.config.gzip}\n")
                f.write(f"FASTQ格式|FASTQ format: {self.config.fastq}\n")
                f.write(f"合并选项|Merge option: {self.config.merge}\n")
                f.write(f"数据库|Database: {self.config.database}\n")
                f.write(f"协议|Protocol: {self.config.protocol}\n\n")

                # 性能参数|Performance parameters
                f.write("3. 性能参数|Performance Parameters\n")
                f.write("-" * 80 + "\n")
                f.write(f"线程数|Threads: {self.config.threads}\n")
                f.write(f"并行连接|Parallel connections: {self.config.parallel}\n")
                if self.config.speed:
                    f.write(f"速度限制|Speed limit: {self.config.speed} MB/s\n")
                f.write(f"使用Aspera|Use Aspera: {self.config.use_aspera}\n\n")

                f.write("=" * 80 + "\n")
                f.write("报告生成完成|Report generation completed\n")
                f.write("=" * 80 + "\n")

            self.logger.info("下载汇总写入完成|Download summary written successfully")

        except Exception as e:
            self.logger.error(f"写入下载汇总失败|Failed to write download summary: {e}")

    def run_analysis(self) -> bool:
        """运行完整分析流程|Run complete analysis pipeline"""
        return self.run_download()
