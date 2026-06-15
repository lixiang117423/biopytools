"""
PlotSR工具函数模块|PlotSR Utility Functions Module
"""

import logging
import sys
import os
import glob
import re
import subprocess
from typing import List, Tuple


class PlotSRLogger:
    """PlotSR日志管理器|PlotSR Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 设置日志级别|Set log level
        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        # 配置logging|Configure logging
        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies() -> Tuple[bool, List[str]]:
    """
    检查依赖工具是否可用|Check if required tools are available

    Returns:
        tuple: (是否全部可用|All available, 缺失工具列表|Missing tools list)
    """
    required_tools = ['minimap2', 'samtools', 'syri', 'plotsr']
    missing_tools = []

    for tool in required_tools:
        try:
            result = subprocess.run(
                ['which', tool],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode != 0:
                missing_tools.append(tool)
        except Exception:
            missing_tools.append(tool)

    return len(missing_tools) == 0, missing_tools


def discover_genomes_in_folder(folder: str) -> List[str]:
    """
    自动发现文件夹内的基因组文件|Auto-discover genome files in folder

    Args:
        folder: 文件夹路径|Folder path

    Returns:
        list: 排序后的基因组文件列表|Sorted list of genome files
    """
    # 查找FASTA文件|Find FASTA files
    extensions = ['*.fa', '*.fa.gz', '*.fasta', '*.fasta.gz', '*.fna', '*.fna.gz']
    fasta_files = []

    for ext in extensions:
        pattern = os.path.join(folder, ext)
        fasta_files.extend(glob.glob(pattern))

    if not fasta_files:
        raise ValueError(f"文件夹内未找到基因组文件|No genome files found in: {folder}")

    # 智能排序|Smart sort
    return smart_sort(fasta_files)


def smart_sort(files: List[str]) -> List[str]:
    """
    智能排序文件列表（支持数字序号）|Smart sort file list (support number sequencing)

    Args:
        files: 文件列表|List of files

    Returns:
        list: 排序后的文件列表|Sorted list of files
    """
    def extract_number(filepath):
        """从文件名提取数字|Extract number from filename"""
        basename = os.path.basename(filepath)
        # 查找所有数字|Find all numbers
        matches = re.findall(r'(\d+)', basename)
        if matches:
            # 使用第一个数字作为排序键|Use first number as sort key
            return int(matches[0])
        return 0

    # 尝试按数字排序|Try to sort by number
    try:
        return sorted(files, key=extract_number)
    except Exception:
        # 失败则按字母序|Fallback to alphabetical sort
        return sorted(files)


def extract_chromosome_lengths(fasta_file: str, output_file: str):
    """
    提取染色体长度|Extract chromosome lengths

    Args:
        fasta_file: FASTA文件路径|FASTA file path
        output_file: 输出文件路径|Output file path
    """
    try:
        # 首先确保索引存在|First ensure index exists
        fai_file = f"{fasta_file}.fai"
        if not os.path.exists(fai_file):
            cmd = ['samtools', 'faidx', fasta_file]
            subprocess.run(cmd, check=True, capture_output=True)

        # 读取索引文件|Read index file
        lines = []
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    lines.append(f"{fields[0]}\t{fields[1]}\n")

        # 写入输出文件|Write output file
        with open(output_file, 'w') as f:
            f.writelines(lines)

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"samtools faidx失败|samtools faidx failed: {e.stderr}")
    except IOError as e:
        raise RuntimeError(f"文件读写失败|File I/O failed: {e}")


def run_command(cmd: List[str], logger: logging.Logger, check: bool = True) -> bool:
    """
    运行命令|Run command

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code

    Returns:
        bool: 是否成功|Whether successful
    """
    try:
        logger.info(f"运行命令|Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            check=check,
            capture_output=True,
            text=True
        )
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"命令执行异常|Command execution error: {e}")
        return False


def format_number(num: int) -> str:
    """
    格式化数字|Format number

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
