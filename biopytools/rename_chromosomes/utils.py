"""
Chromosome Rename Utilities
染色体重命名工具类
"""

import logging
import sys
from datetime import datetime
import os


class ChromosomeRenameLogger:
    """染色体重命名日志管理器 | Chromosome Rename Logger Manager"""

    def __init__(self, output_dir: str):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            output_dir: 输出目录 | Output directory
        """
        self.output_dir = output_dir

        # 创建日志目录 | Create log directory
        os.makedirs(output_dir, exist_ok=True)

        # 设置日志文件名 | Set log file name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(output_dir, f"rename_chromosomes_{timestamp}.log")

        # 配置日志处理器 | Configure log handlers
        self._setup_logging(log_file)

    def _setup_logging(self, log_file: str):
        """
        设置日志 | Setup logging

        Args:
            log_file: 日志文件路径 | Log file path
        """
        # 创建logger | Create logger
        self.logger = logging.getLogger("ChromosomeRename")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 创建formatter | Create formatter
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout处理器 - 只输出INFO及以上 | stdout handler - INFO and above only
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr处理器 - 输出WARNING及以上 | stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 文件处理器 - 输出所有级别 | File handler - all levels
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # 添加处理器 | Add handlers
        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class FastaRenamer:
    """FASTA序列重命名器 | FASTA Sequence Renamer"""

    def __init__(self, logger):
        """
        初始化重命名器 | Initialize renamer

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def rename_sequences(
        self,
        input_file: str,
        output_file: str,
        chromosome_number: int
    ) -> bool:
        """
        重命名FASTA序列 | Rename FASTA sequences

        Args:
            input_file: 输入FASTA文件 | Input FASTA file
            output_file: 输出FASTA文件 | Output FASTA file
            chromosome_number: 染色体数量 | Number of chromosomes

        Returns:
            是否成功 | Whether successful
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始重命名FASTA序列")
            self.logger.info("=" * 60)
            self.logger.info(f"输入文件: {input_file}")
            self.logger.info(f"输出文件: {output_file}")
            self.logger.info(f"染色体数量: {chromosome_number}")

            # 使用awk进行流式处理 | Use awk for stream processing
            awk_script = f'''
BEGIN {{
    chr_count = 0
    scaf_count = 0
}}
/^>/ {{
    chr_count++
    if (chr_count <= {chromosome_number}) {{
        printf ">Chr%02d\\n", chr_count
    }} else {{
        scaf_count++
        printf ">HiC_scaffold_%02d\\n", scaf_count
    }}
    next
}}
{{ print }}
'''

            # 执行awk命令 | Execute awk command
            import subprocess
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                result = subprocess.run(
                    ['awk', awk_script],
                    stdin=infile,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True
                )

                if result.returncode != 0:
                    self.logger.error(f"AWK处理失败: {result.stderr}")
                    return False

            # 统计序列数量 | Count sequences
            total_seqs = self._count_sequences(output_file)
            scaffold_total = total_seqs - chromosome_number

            self.logger.info("=" * 60)
            self.logger.info("重命名完成")
            self.logger.info("=" * 60)
            self.logger.info(f"总序列数: {total_seqs}")
            self.logger.info(f"染色体: {chromosome_number}")
            self.logger.info(f"Scaffolds: {scaffold_total}")

            return True

        except Exception as e:
            self.logger.error(f"重命名失败: {str(e)}", exc_info=True)
            return False

    def _count_sequences(self, fasta_file: str) -> int:
        """
        统计FASTA序列数量 | Count FASTA sequences

        Args:
            fasta_file: FASTA文件路径 | FASTA file path

        Returns:
            序列数量 | Number of sequences
        """
        count = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count

    def preview_output(self, output_file: str, lines: int = 10):
        """
        预览输出文件 | Preview output file

        Args:
            output_file: 输出文件路径 | Output file path
            lines: 显示行数 | Number of lines to show
        """
        self.logger.info("=" * 60)
        self.logger.info(f"输出文件预览 (前{lines}行)")
        self.logger.info("=" * 60)

        with open(output_file, 'r') as f:
            for i, line in enumerate(f):
                if i >= lines:
                    break
                self.logger.info(line.rstrip())
