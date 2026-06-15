"""
序列提取模块|Sequence Extraction Module
根据BED文件从基因组提取对应的序列|Extract sequences from genome based on BED file
"""

import os
import subprocess
from typing import List
from .parser import PAFRecord


class SequenceExtractor:
    """序列提取器|Sequence Extractor"""

    def __init__(self, config, logger, cmd_runner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_sequences(self, bed_file: str, output_file: str = None):
        """根据BED文件从基因组提取序列|Extract sequences from genome based on BED file

        Args:
            bed_file: BED文件路径|BED file path
            output_file: 输出FASTA文件路径|Output FASTA file path
        """
        if output_file is None:
            output_file = os.path.join(self.config.output_dir, "alignment.fa")

        self.logger.info(f"提取基因组序列|Extracting genome sequences")
        self.logger.info(f"  BED文件|BED file: {bed_file}")
        self.logger.info(f"  基因组|Genome: {self.config.genome_fa}")
        self.logger.info(f"  输出文件|Output file: {output_file}")

        # 构建seqkit命令|Build seqkit command
        cmd = f"seqkit subseq --bed {bed_file} {self.config.genome_fa} > {output_file}"

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, " 序列提取|Sequence extraction")

        if success:
            # 统计提取的序列数|Count extracted sequences
            try:
                result = subprocess.run(
                    ["grep", "-c", "^>", output_file],
                    capture_output=True,
                    text=True,
                    check=False
                )
                seq_count = int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0
                self.logger.info(f" 成功提取{seq_count}条序列|Successfully extracted {seq_count} sequences")
            except Exception as e:
                self.logger.warning(f" 无法统计序列数量|Cannot count sequences: {e}")

            self.logger.info(f" 序列文件已保存|Sequence file saved: {output_file}")

        return success if success else False

# ===== END FILE =====
