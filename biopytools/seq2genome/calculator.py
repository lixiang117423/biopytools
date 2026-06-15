"""
序列比对核心逻辑模块|Sequence Alignment Core Logic Module
支持Minimap2（DNA）和Miniprot（蛋白质）|Support Minimap2 (DNA) and Miniprot (Protein)
"""

import os
from typing import Optional
from .utils import CommandRunner


class MiniprotAligner:
    """Miniprot对齐器|Miniprot Aligner"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_alignment(self, output_paf: Optional[str] = None) -> str:
        """运行miniprot对齐|Run miniprot alignment

        Args:
            output_paf: 输出PAF文件路径|Output PAF file path
                       如果为None，使用默认路径|If None, use default path

        Returns:
            str: PAF文件路径|PAF file path
        """
        # 确定输出文件路径|Determine output file path
        if output_paf is None:
            output_paf = os.path.join(self.config.output_dir, "alignment.paf")

        self.logger.info("=" * 60)
        self.logger.info("运行Miniprot对齐|Running Miniprot alignment")
        self.logger.info("=" * 60)
        self.logger.info(f"基因组|Genome: {self.config.genome_fa}")
        self.logger.info(f"蛋白质|Protein: {self.config.query_fa}")
        self.logger.info(f"输出|Output: {output_paf}")
        self.logger.info(f"线程|Threads: {self.config.threads}")

        # 构建miniprot命令|Build miniprot command
        cmd = self._build_miniprot_command(output_paf)

        # 执行命令|Execute command
        timeout = 86400  # 24小时超时|24 hours timeout
        success = self.cmd_runner.run(cmd, "Miniprot对齐|Miniprot alignment", timeout=timeout)

        if not success:
            raise RuntimeError("Miniprot对齐失败|Miniprot alignment failed")

        # 验证输出文件|Verify output file
        if not os.path.exists(output_paf):
            raise RuntimeError(f"PAF文件未生成|PAF file not generated: {output_paf}")

        file_size = os.path.getsize(output_paf)
        self.logger.info(f"PAF文件已生成|PAF file generated: {output_paf} ({file_size} bytes)")

        return output_paf

    def _build_miniprot_command(self, output_paf: str) -> str:
        """构建miniprot命令|Build miniprot command

        Args:
            output_paf: 输出PAF文件|Output PAF file

        Returns:
            str: miniprot命令|miniprot command
        """
        cmd = (
            f"{self.config.miniprot_path} "
            f"{self.config.genome_fa} "
            f"{self.config.query_fa} "
            f"-t {self.config.threads} "
            f"> {output_paf}"
        )

        return cmd


class Minimap2Aligner:
    """Minimap2对齐器（用于DNA序列）|Minimap2 Aligner (for DNA sequences)"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_alignment(self, output_paf: Optional[str] = None) -> str:
        """运行minimap2对齐|Run minimap2 alignment

        Args:
            output_paf: 输出PAF文件路径|Output PAF file path
                       如果为None，使用默认路径|If None, use default path

        Returns:
            str: PAF文件路径|PAF file path
        """
        # 确定输出文件路径|Determine output file path
        if output_paf is None:
            output_paf = os.path.join(self.config.output_dir, "alignment.paf")

        self.logger.info("=" * 60)
        self.logger.info("运行Minimap2对齐|Running Minimap2 alignment")
        self.logger.info("=" * 60)
        self.logger.info(f"基因组|Genome: {self.config.genome_fa}")
        self.logger.info(f"查询序列|Query sequences: {self.config.query_fa}")
        self.logger.info(f"输出|Output: {output_paf}")
        self.logger.info(f"线程|Threads: {self.config.threads}")

        # 构建minimap2命令|Build minimap2 command
        cmd = self._build_minimap2_command(output_paf)

        # 执行命令|Execute command
        timeout = 86400  # 24小时超时|24 hours timeout
        success = self.cmd_runner.run(cmd, "Minimap2对齐|Minimap2 alignment", timeout=timeout)

        if not success:
            raise RuntimeError("Minimap2对齐失败|Minimap2 alignment failed")

        # 验证输出文件|Verify output file
        if not os.path.exists(output_paf):
            raise RuntimeError(f"PAF文件未生成|PAF file not generated: {output_paf}")

        file_size = os.path.getsize(output_paf)
        self.logger.info(f"PAF文件已生成|PAF file generated: {output_paf} ({file_size} bytes)")

        return output_paf

    def _build_minimap2_command(self, output_paf: str) -> str:
        """构建minimap2命令|Build minimap2 command

        Args:
            output_paf: 输出PAF文件|Output PAF file

        Returns:
            str: minimap2命令|minimap2 command
        """
        # 使用asm5 preset (for mapping similar sequences)
        # 也可以使用asm10, asm20等preset depending on divergence
        cmd = (
            f"{self.config.minimap2_path} "
            f"-x asm5 "  # preset for mapping sequences with ~5% divergence
            f"-t {self.config.threads} "
            f"--secondary=no "  # 只输出primary alignment|Only output primary alignments
            f"{self.config.genome_fa} "
            f"{self.config.query_fa} "
            f"> {output_paf}"
        )

        return cmd
