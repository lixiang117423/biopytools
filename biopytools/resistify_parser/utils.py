"""
Resistify Parser工具函数模块|Resistify Parser Utility Functions Module
"""

import pandas as pd
import logging
import sys
from pathlib import Path


class ResistifyParserLogger:
    """Resistify Parser日志管理器|Resistify Parser Logger Manager"""

    def __init__(self, log_file: str = "resistify_parser_analysis.log"):
        self.log_file = log_file
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def validate_input_directory(input_dir: str) -> bool:
    """
    验证输入目录是否存在|Validate input directory exists

    Args:
        input_dir: 输入目录路径|Input directory path

    Returns:
        bool: 目录是否存在|Whether directory exists
    """
    return Path(input_dir).is_dir()


def check_required_files(input_dir: str, required_files: list) -> tuple:
    """
    检查必需文件是否存在|Check required files exist

    Args:
        input_dir: 输入目录路径|Input directory path
        required_files: 必需文件列表|Required files list

    Returns:
        tuple: (all_exist, missing_files)|Whether all exist, list of missing files
    """
    input_path = Path(input_dir)
    missing = []

    for file in required_files:
        file_path = input_path / file
        if not file_path.exists():
            missing.append(file)

    return (len(missing) == 0, missing)


def format_classification(classification: str) -> str:
    """
    格式化分类名称|Format classification name

    Args:
        classification: 分类字符串|Classification string

    Returns:
        str: 格式化后的分类|Formatted classification
    """
    return classification.strip()


def create_output_filename(prefix: str, suffix: str, output_dir: str = '.', extension: str = 'csv') -> Path:
    """
    创建输出文件名|Create output filename

    Args:
        prefix: 输出前缀|Output prefix
        suffix: 输出后缀|Output suffix
        output_dir: 输出目录|Output directory
        extension: 文件扩展名|File extension

    Returns:
        Path: 输出文件路径|Output file path
    """
    filename = f"{prefix}_{suffix}.{extension}"
    return Path(output_dir) / filename


def print_summary_statistics(df: pd.DataFrame, logger=None):
    """
    打印汇总统计信息|Print summary statistics

    Args:
        df: 数据框|DataFrame
        logger: 日志对象|Logger object (optional)
    """
    msg = [
        "="*60,
        "统计摘要|Statistics Summary",
        "="*60,
        f"总NLR基因数|Total NLR genes: {len(df)}",
    ]

    # 分类统计|Classification statistics
    if 'Classification' in df.columns:
        class_counts = df['Classification'].value_counts()
        msg.append("\n分类统计|Classification statistics:")
        for cls, count in class_counts.items():
            msg.append(f"  {cls}: {count}")

    # Domain统计|Domain statistics
    domain_cols = ['has_TIR', 'has_NB_ARC', 'has_CC', 'has_LRR', 'has_RPW8']
    for col in domain_cols:
        if col in df.columns:
            count = df[col].sum()
            msg.append(f"\n{col}: {count}条|{count} sequences")

    msg.append("="*60)

    summary = "\n".join(msg)
    if logger:
        logger.info(summary)
    else:
        print(summary)
