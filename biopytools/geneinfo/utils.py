"""
GFF3工具函数模块|GFF3 Utility Functions Module
"""

import logging
import sys
import time
from pathlib import Path
from typing import Dict, Any


class GFFLogger:
    """GFF3处理日志管理器|GFF3 Processing Logger Manager"""

    def __init__(self, output_file: str, verbose: bool = False, quiet: bool = False):
        self.output_dir = Path(output_file).parent
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f"gff_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"gff_processing_{timestamp}")

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

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

    def step(self, message: str):
        """记录步骤|Log step"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)

class AttributeParser:
    """属性解析器|Attribute Parser"""

    def __init__(self, logger):
        self.logger = logger

    def parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        解析GFF3第九列的属性字符串|Parse GFF3 column 9 attributes string

        Args:
            attr_string: 属性字符串，如 "ID=gene1;Name=GeneA"|Attribute string, e.g., "ID=gene1;Name=GeneA"

        Returns:
            Dict[str, str]: 属性字典，如 {'ID': 'gene1', 'Name': 'GeneA'}|Attribute dictionary, e.g., {'ID': 'gene1', 'Name': 'GeneA'}
        """
        attributes = {}

        if not attr_string or attr_string.strip() == '.':
            return attributes

        try:
            for part in attr_string.strip().split(';'):
                if '=' in part:
                    key, value = part.split('=', 1)
                    attributes[key.strip()] = value.strip()
        except Exception as e:
            self.logger.warning(f"解析属性字符串时出错|Error parsing attributes: {attr_string}, {e}")

        return attributes

class GFFValidator:
    """GFF3文件验证器|GFF3 File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def validate_gff_line(self, line: str) -> bool:
        """
        验证GFF3行格式|Validate GFF3 line format

        Args:
            line: GFF3行|GFF3 line

        Returns:
            bool: 是否有效|Whether valid
        """
        if line.startswith('#'):
            return False

        parts = line.strip().split('\t')
        if len(parts) != 9:
            return False

        # 检查必需字段|Check required fields
        try:
            # 检查起始和结束位置是否为数字|Check if start and end positions are numeric
            int(parts[3])  # start
            int(parts[4])  # end

            # 检查链方向|Check strand
            if parts[6] not in ['+', '-', '.', '?']:
                return False

        except (ValueError, IndexError):
            return False

        return True

    def check_file_format(self, file_path: str) -> bool:
        """
        检查文件格式|Check file format

        Args:
            file_path: 文件路径|File path

        Returns:
            bool: 是否为有效的GFF3文件|Whether valid GFF3 file
        """
        try:
            with open(file_path, 'r') as f:
                # 检查前几行|Check first few lines
                for i, line in enumerate(f):
                    if i > 100:  # 只检查前100行|Only check first 100 lines
                        break

                    if line.startswith('##gff-version'):
                        return True

                    if not line.startswith('#') and line.strip():
                        # 检查数据行格式|Check data line format
                        if self.validate_gff_line(line):
                            return True
                        else:
                            self.logger.warning(f"可能的格式错误|Possible format error at line {i+1}: {line.strip()}")
                            return False

            return True

        except Exception as e:
            self.logger.error(f"文件格式检查失败|File format check failed: {e}")
            return False
