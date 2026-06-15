"""
Miniprot对齐核心逻辑模块|Miniprot Alignment Core Logic Module
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
        self.logger.info(f"蛋白质|Protein: {self.config.protein_fa}")
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
            f"{self.config.protein_fa} "
            f"-t {self.config.threads} "
            f"> {output_paf}"
        )

        return cmd
