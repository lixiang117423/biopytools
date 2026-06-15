"""
MAFFT多序列比对模块|MAFFT Multiple Sequence Alignment Module
"""

import os
import subprocess
from pathlib import Path
from typing import Optional


class MAFFTAligner:
    """MAFFT多序列比对器|MAFFT Multiple Sequence Aligner"""

    def __init__(self, mafft_path='mafft', threads=4, logger=None):
        """
        初始化MAFFT比对器|Initialize MAFFT aligner

        Parameters
        ----------
        mafft_path : str
            MAFFT可执行文件路径|MAFFT executable path
        threads : int
            线程数|Number of threads
        logger : logging.Logger
            日志器|Logger object
        """
        self.mafft_path = mafft_path
        self.threads = threads
        self.logger = logger

    def check_mafft(self) -> bool:
        """
        检查MAFFT是否可用|Check if MAFFT is available

        Returns
        -------
        bool
            MAFFT是否可用|Whether MAFFT is available
        """
        try:
            result = subprocess.run(
                [self.mafft_path, '--help'],
                capture_output=True,
                text=True,
                timeout=10
            )
            # MAFFT --help 返回非0退出码，检查输出是否包含MAFFT标识
            return 'MAFFT' in result.stdout or 'MAFFT' in result.stderr
        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            return False

    def run_alignment(
        self,
        input_file: str | Path,
        output_file: str | Path,
        mafft_params: str = '--auto'
    ) -> bool:
        """
        运行MAFFT比对|Run MAFFT alignment

        Parameters
        ----------
        input_file : str or Path
            输入FASTA文件|Input FASTA file
        output_file : str or Path
            输出比对文件|Output alignment file
        mafft_params : str
            MAFFT参数|MAFFT parameters (default: --auto)

        Returns
        -------
        bool
            比对是否成功|Whether alignment was successful
        """
        input_file = Path(input_file)
        output_file = Path(output_file)

        # 检查输入文件|Check input file
        if not input_file.exists():
            if self.logger:
                self.logger.error(f"输入文件不存在|Input file does not exist: {input_file}")
            return False

        # 创建输出目录|Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # 构建MAFFT命令|Build MAFFT command
        cmd = f"{self.mafft_path} {mafft_params} --thread {self.threads} "
        cmd += f"{input_file.resolve()} > {output_file.resolve()}"

        if self.logger:
            self.logger.info(f"运行MAFFT比对|Running MAFFT alignment")
            self.logger.info(f"命令|Command: {self.mafft_path} {mafft_params} --thread {self.threads}")

        try:
            # 运行MAFFT|Run MAFFT
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            if result.returncode == 0:
                # 检查输出文件|Check output file
                if output_file.exists() and output_file.stat().st_size > 0:
                    if self.logger:
                        self.logger.info(f"比对完成|Alignment completed: {output_file}")
                    return True
                else:
                    if self.logger:
                        self.logger.error(f"输出文件为空|Output file is empty: {output_file}")
                    return False
            else:
                if self.logger:
                    self.logger.error(f"MAFFT执行失败|MAFFT execution failed")
                    if result.stderr:
                        self.logger.error(f"错误信息|Error: {result.stderr}")
                return False

        except subprocess.TimeoutExpired:
            if self.logger:
                self.logger.error(f"MAFFT执行超时|MAFFT execution timeout")
            return False
        except Exception as e:
            if self.logger:
                self.logger.error(f"MAFFT执行异常|MAFFT execution error: {e}")
            return False

    def detect_sequence_type(self, fasta_file: str | Path) -> Optional[str]:
        """
        检测序列类型|Detect sequence type

        Parameters
        ----------
        fasta_file : str or Path
            FASTA文件路径|FASTA file path

        Returns
        -------
        str or None
            序列类型 ('protein' 或 'nucleotide')|Sequence type
        """
        fasta_file = Path(fasta_file)
        if not fasta_file.exists():
            return None

        try:
            with open(fasta_file, 'r') as f:
                sequence = ''
                for line in f:
                    if not line.startswith('>'):
                        sequence += line.strip()
                    if len(sequence) >= 1000:  # 只检测前1000个碱基/氨基酸
                        break

            if not sequence:
                return None

            # 统计ATGCUN比例|Count ATGCUN ratio
            atgc_count = sum(1 for c in sequence.upper() if c in 'ATGCUN')
            total_count = sum(1 for c in sequence.upper() if c in 'ATGCUN-BDEFHIJKLMPQRSVWXYZ*')

            if total_count == 0:
                return None

            ratio = atgc_count / total_count

            # 如果ATGCUN比例>90%，判断为核苷酸|If ATGCUN ratio > 90%, it's nucleotide
            if ratio > 0.9:
                return 'nucleotide'
            else:
                return 'protein'

        except Exception as e:
            if self.logger:
                self.logger.warning(f"序列类型检测失败|Sequence type detection failed: {e}")
            return None
