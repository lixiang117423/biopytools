"""
GEC校正工具函数模块|GEC Correction Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path


class GECLogger:
    """GEC校正日志管理器|GEC Correction Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 删除旧日志文件|Delete old log file
        if self.log_file and os.path.exists(self.log_file):
            os.remove(self.log_file)

        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir or Path.cwd()

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


def format_number(num: float) -> str:
    """格式化数字|Format number"""
    if isinstance(num, (int, float)):
        if num >= 1_000_000:
            return f"{num / 1_000_000:.2f}M"
        elif num >= 1_000:
            return f"{num / 1_000:.2f}K"
        else:
            return f"{num:.2f}"
    return str(num)


def format_scientific(num: float) -> str:
    """格式化科学计数法|Format scientific notation"""
    return f"{num:.2e}"


def validate_chromosome_format(pfile: str, chrom_col: str = 'CHR',
                               logger = None) -> tuple:
    """
    验证并检查染色体编号格式|Validate and check chromosome format

    Args:
        pfile: P值文件路径|P-value file path
        chrom_col: 染色体列名|Chromosome column name
        logger: 日志器|Logger

    Returns:
        (is_valid, format_type, examples): (是否有效, 格式类型, 示例)
    """
    import gzip

    chrom_formats = {
        'numeric': [],      # 纯数字: 1, 2, 22
        'chr_prefix': [],   # chr前缀: chr1, chr2, chrX
        'Chr_prefix': [],   # Chr前缀: Chr1, Chr2, ChrX
        'other': []         # 其他格式
    }

    try:
        open_func = gzip.open if pfile.endswith('.gz') else open

        with open_func(pfile, 'rt') as f:
            # 读取表头|Read header
            header = f.readline().strip()

            if not header:
                return False, 'empty_file', []

            # 查找染色体列索引|Find chromosome column index
            columns = header.split('\t')
            try:
                chrom_idx = columns.index(chrom_col)
            except ValueError:
                # 尝试不区分大小写匹配|Try case-insensitive match
                chrom_idx = None
                for i, col in enumerate(columns):
                    if col.upper() == chrom_col.upper():
                        chrom_idx = i
                        break

                if chrom_idx is None:
                    return False, 'column_not_found', []

            # 读取前100行数据检查格式|Read first 100 lines to check format
            sample_chroms = set()
            for line_num, line in enumerate(f, 2):
                if line_num > 102:  # 表头+100行数据|header + 100 data lines
                    break

                line = line.strip()
                if not line:
                    continue

                fields = line.split('\t')
                if len(fields) > chrom_idx:
                    chrom = fields[chrom_idx].strip()
                    if chrom:
                        sample_chroms.add(chrom)

                        # 分类染色体格式|Classify chromosome format
                        if chrom.isdigit():
                            chrom_formats['numeric'].append(chrom)
                        elif chrom.startswith('chr') or chrom.startswith('Chr'):
                            if chrom.startswith('chr'):
                                chrom_formats['chr_prefix'].append(chrom)
                            else:
                                chrom_formats['Chr_prefix'].append(chrom)
                        else:
                            chrom_formats['other'].append(chrom)

        # 判断主要格式类型|Determine main format type
        format_counts = {
            k: len(set(v)) for k, v in chrom_formats.items()
        }

        main_format = max(format_counts, key=format_counts.get)

        if format_counts[main_format] == 0:
            return False, 'no_data', []

        examples = list(set(chrom_formats[main_format]))[:5]

        if logger:
            logger.info(f"染色体编号格式检测|Chromosome format detection: {main_format}")
            logger.info(f"样本染色体|Sample chromosomes: {', '.join(examples)}")

        return True, main_format, examples

    except Exception as e:
        if logger:
            logger.error(f"染色体格式检测失败|Chromosome format detection failed: {str(e)}")
        return False, 'error', []
