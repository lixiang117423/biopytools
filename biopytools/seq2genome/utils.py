"""
序列到基因组比对工具函数模块|Sequence to Genome Alignment Utility Functions Module
"""

import os
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Tuple
from Bio import SeqIO


class Pep2GenomeLogger:
    """蛋白质到基因组比对日志管理器|Protein to Genome Alignment Logger Manager"""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = output_dir / f"pep2genome_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"pep2genome_processing_{timestamp}")

        # 设置日志级别|Set log level
        self.logger.setLevel(logging.DEBUG)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        self.logger = logger

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间（秒），None表示无限制|Timeout in seconds, None means no limit

        Returns:
            bool: 执行成功返回True，失败返回False|True if successful, False otherwise
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        if timeout:
            self.logger.info(f"超时设置|Timeout: {timeout}秒|seconds ({timeout/3600:.1f}小时|hours)")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            self.logger.info(f"{description} 完成|completed")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error(f"{description} 超时|timed out after {timeout}秒|seconds")
            return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class FileValidator:
    """文件验证器|File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在|Check if file exists

        Args:
            file_path: 文件路径|File path
            description: 文件描述|File description

        Returns:
            bool: 文件存在返回True，否则返回False|True if exists, False otherwise
        """
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False


def detect_sequence_type(fasta_file: str, sample_size: int = 10000, check_sequences: int = 5) -> Tuple[str, float]:
    """
    自动检测FASTA文件中的序列类型（DNA或蛋白质）|Auto-detect sequence type in FASTA file (DNA or protein)

    检测逻辑|Detection logic:
    - DNA序列: 只包含ATCGN字符，且ATCG含量高
    - 蛋白质序列: 包含20种氨基酸字符

    Args:
        fasta_file: FASTA文件路径|FASTA file path
        sample_size: 每条序列采样的字符数|Sample characters per sequence
        check_sequences: 检查的序列数量|Number of sequences to check

    Returns:
        Tuple[str, float]: (序列类型|Sequence type, 置信度|Confidence)
                          - 'dna' 或 'protein' 或 'unknown'
                          - 置信度 0.0-1.0
    """
    logger = logging.getLogger(__name__)

    try:
        # 读取序列|Read sequences
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq).upper())
            if len(sequences) >= check_sequences:
                break

        if not sequences:
            logger.error("FASTA文件中没有序列|No sequences found in FASTA file")
            return "unknown", 0.0

        # 统计字符|Count characters
        total_chars = 0
        dna_chars = 0
        protein_chars = 0
        ambiguous_chars = 0

        for seq in sequences[:check_sequences]:
            sample_seq = seq[:sample_size] if len(seq) > sample_size else seq
            total_chars += len(sample_seq)

            for char in sample_seq:
                if char in "ATCG":
                    dna_chars += 1
                elif char in "N":
                    ambiguous_chars += 1
                elif char in "DEFGHIKLMNPQRSTVWY":
                    protein_chars += 1
                else:
                    # 其他字符（如BZJOUX*等 ambiguous amino acids）
                    protein_chars += 1

        # 计算比例|Calculate ratios
        if total_chars == 0:
            return "unknown", 0.0

        dna_ratio = (dna_chars + ambiguous_chars) / total_chars
        protein_ratio = protein_chars / total_chars

        logger.debug(f"序列统计|Sequence statistics:")
        logger.debug(f"  总字符数|Total chars: {total_chars}")
        logger.debug(f"  DNA字符|DNA chars (ATCGN): {dna_chars + ambiguous_chars} ({dna_ratio*100:.1f}%)")
        logger.debug(f"  蛋白质字符|Protein chars: {protein_chars} ({protein_ratio*100:.1f}%)")

        # 判断序列类型|Determine sequence type
        # 如果DNA字符比例 > 90%，判定为DNA
        # If DNA char ratio > 90%, classify as DNA
        if dna_ratio > 0.90:
            confidence = dna_ratio
            logger.info(f"检测到DNA序列|Detected DNA sequence (置信度|confidence: {confidence:.2f})")
            return "dna", confidence

        # 如果蛋白质字符比例 > 30%，判定为蛋白质
        # If protein char ratio > 30%, classify as protein
        # (蛋白质序列通常含有更多样化的氨基酸)
        elif protein_ratio > 0.30:
            confidence = protein_ratio
            logger.info(f"检测到蛋白质序列|Detected protein sequence (置信度|confidence: {confidence:.2f})")
            return "protein", confidence

        else:
            logger.warning(f"无法确定序列类型|Cannot determine sequence type (DNA: {dna_ratio:.2f}, Protein: {protein_ratio:.2f})")
            return "unknown", 0.0

    except Exception as e:
        logger.error(f"序列类型检测失败|Sequence type detection failed: {e}")
        return "unknown", 0.0
