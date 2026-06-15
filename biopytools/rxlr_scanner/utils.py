"""
RxLR扫描工具函数模块|RxLR Scanner Utility Functions Module
"""

import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterator, Optional


class RxLRLogger:
    """RxLR扫描日志管理器|RxLR Scanner Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, verbose: bool = False):
        self.log_file = log_file
        self.verbose = verbose
        self.logger = None

        self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器|Setup logger"""
        self.logger = logging.getLogger("rxlr_scanner")
        self.logger.setLevel(logging.DEBUG)

        # 清除现有处理器|Clear existing handlers
        self.logger.handlers.clear()

        # 日志格式|Log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 控制台处理器|Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO if not self.verbose else logging.DEBUG)
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # 文件处理器|File handler
        if self.log_file:
            try:
                log_dir = os.path.dirname(os.path.abspath(self.log_file))
                Path(log_dir).mkdir(parents=True, exist_ok=True)

                file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
                file_handler.setLevel(logging.DEBUG)
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)
            except Exception as e:
                self.logger.warning(f"无法创建日志文件|Cannot create log file: {e}")

    def get_logger(self):
        """获取日志记录器|Get logger"""
        return self.logger


class FastaParser:
    """FASTA文件解析器|FASTA File Parser"""

    @staticmethod
    def parse_fasta(file_path: str) -> Iterator[Tuple[str, str]]:
        """
        解析FASTA文件|Parse FASTA file

        Args:
            file_path: FASTA文件路径|FASTA file path

        Yields:
            (序列ID, 序列)|(sequence ID, sequence)
        """
        current_id = None
        current_seq = []

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()

                    if not line:
                        continue

                    if line.startswith('>'):
                        # 保存前一条序列|Save previous sequence
                        if current_id is not None:
                            yield current_id, ''.join(current_seq)

                        # 开始新序列|Start new sequence
                        current_id = line[1:].split()[0]  # 取第一个词作为ID
                        current_seq = []
                    else:
                        # 累积序列|Accumulate sequence
                        current_seq.append(line)

                # 保存最后一条序列|Save last sequence
                if current_id is not None:
                    yield current_id, ''.join(current_seq)

        except Exception as e:
            raise IOError(f"读取FASTA文件失败|Failed to read FASTA file: {e}")

    @staticmethod
    def validate_sequence(seq_id: str, sequence: str, min_length: int) -> Tuple[bool, str]:
        """
        验证序列|Validate sequence

        Args:
            seq_id: 序列ID|Sequence ID
            sequence: 序列|Sequence
            min_length: 最小长度|Minimum length

        Returns:
            (是否有效, 错误信息)|(is_valid, error_message)
        """
        if not sequence:
            return False, f"序列为空|Empty sequence: {seq_id}"

        # 检查是否为标准氨基酸|Check for valid amino acids
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
        sequence_upper = sequence.upper()

        invalid_chars = set(sequence_upper) - valid_aa
        if invalid_chars and '*' not in invalid_chars:
            # 允许终止符*，但不允许其他非法字符
            return False, f"序列包含非法氨基酸|Sequence contains invalid amino acids: {seq_id} - {invalid_chars}"

        # 检查长度|Check length
        if len(sequence) < min_length:
            return False, f"序列长度小于{min_length}|Sequence length < {min_length}: {seq_id}"

        return True, ""


class SequenceExtractor:
    """序列片段提取器|Sequence Fragment Extractor"""

    @staticmethod
    def extract_window(sequence: str, start: int, end: int) -> Optional[str]:
        """
        提取序列窗口|Extract sequence window

        Args:
            sequence: 完整序列|Full sequence
            start: 起始位置（0-based）|Start position (0-based)
            end: 结束位置（不包含）|End position (exclusive)

        Returns:
            提取的序列片段|Extracted sequence fragment, or None if out of range
        """
        seq_len = len(sequence)

        if start >= seq_len:
            return None

        # 调整end不超过序列长度|Adjust end not to exceed sequence length
        actual_end = min(end, seq_len)

        return sequence[start:actual_end]

    @staticmethod
    def search_motif(sequence: str, motif: str) -> List[int]:
        """
        搜索基序位置|Search motif positions

        Args:
            sequence: 序列|Sequence
            motif: 基序|Motif

        Returns:
            匹配位置列表（0-based）|List of match positions (0-based)
        """
        positions = []
        seq_len = len(sequence)
        motif_len = len(motif)

        for i in range(seq_len - motif_len + 1):
            if sequence[i:i+motif_len] == motif:
                positions.append(i)

        return positions
