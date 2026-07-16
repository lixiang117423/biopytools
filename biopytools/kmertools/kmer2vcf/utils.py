"""
工具函数模块|Utility Functions Module
"""

import gzip
import logging
import sys
import subprocess
# zstandard 延迟导入: 仅处理 .zst 文件时需要,避免硬依赖导致整个模块无法导入
# zstandard is lazily imported only when handling .zst files (avoids hard dependency blocking module import)


def open_input(filepath, logger=None):
    """根据后缀自动选择打开方式，支持 .gz / .zst / 普通文本|Open file with auto-detection of compression format

    Args:
        filepath: 文件路径|File path
        logger: 日志器|Logger

    Returns:
        文件对象|File object in text mode
    """
    if filepath.endswith('.gz'):
        if logger:
            logger.info("检测到gzip压缩文件|Detected gzip compressed file")
        return gzip.open(filepath, 'rt')
    elif filepath.endswith('.zst'):
        if logger:
            logger.info("检测到zstd压缩文件|Detected zstd compressed file")
        try:
            import zstandard as zstd
        except ImportError:
            raise ImportError(
                "处理 .zst 文件需要 zstandard 库|zstandard library required to handle .zst files. "
                "请安装|Please install: conda install -c conda-forge zstandard  或|or  pip install zstandard"
            )
        dctx = zstd.ZstdDecompressor()
        with open(filepath, 'rb') as f_compressed:
            decompressed = dctx.decompress(f_compressed.read(), max_output_size=5 * 1024 * 1024 * 1024)
        import io
        return io.TextIOWrapper(io.BytesIO(decompressed), encoding='utf-8')
    else:
        return open(filepath, 'r')


class Kmer2VcfLogger:
    """Kmer转VCF日志管理器|Kmer to VCF Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

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
    else:
        return str(num)


def print_progress(current: int, total: int, prefix: str = "Processed"):
    """
    打印进度|Print progress

    Args:
        current: 当前进度|Current progress
        total: 总数|Total count
        prefix: 前缀文本|Prefix text
    """
    if total == 0:
        return

    percent = (current / total) * 100
    print(f"\r{prefix}: {format_number(current)}/{format_number(total)} ({percent:.2f}%)", end='', flush=True)


def print_interval_progress(current: int, interval: int = 1000000):
    """
    按间隔打印进度|Print progress at intervals

    Args:
        current: 当前进度|Current progress
        interval: 间隔|Interval
    """
    if current % interval == 0:
        print(f"  Processed {format_number(current)} lines...", end='\r')


def sort_file_by_kmer(input_file, output_file, temp_dir, sort_memory='10G', logger=None):
    """
    按kmer ID排序文件|Sort file by kmer ID

    Args:
        input_file: 输入文件路径|Input file path
        output_file: 输出文件路径|Output file path
        temp_dir: 临时文件目录|Temporary directory for sort
        sort_memory: 排序最大内存使用|Maximum memory for sorting (default: 10G)
        logger: 日志器|Logger object

    Returns:
        bool: 是否成功|Success status
    """
    if logger:
        logger.info("=" * 100)
        logger.info(" Step 2: 按kmer ID排序文件|Step 2: Sorting file by kmer ID")
        logger.info("=" * 100)

    try:
        # 构建sort命令|Build sort command
        cmd = [
            'sort',
            '-T', temp_dir,      # 临时文件目录|Temporary directory
            '-S', sort_memory,   # 最大内存|Maximum memory
            '-k1,1',             # 按第一列排序|Sort by first column
            input_file,
            '-o', output_file
        ]

        if logger:
            logger.info(f"排序命令|Sort command: {' '.join(cmd)}")
            logger.info(f"内存限制|Memory limit: {sort_memory}")

        # 执行排序|Execute sorting
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if logger:
            logger.info(f"排序完成|Sorting completed")
            logger.info(f"排序后文件|Sorted file: {output_file}")

        return True

    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"排序失败|Sorting failed: {e.stderr}")
        return False
    except Exception as e:
        if logger:
            logger.error(f"排序异常|Sorting error: {e}")
        return False


def compress_vcf_with_bgzip(input_file, output_file, logger=None):
    """
    使用bgzip压缩VCF文件|Compress VCF file with bgzip

    Args:
        input_file: 输入VCF文件|Input VCF file
        output_file: 输出压缩文件|Output compressed file (.vcf.gz)
        logger: 日志器|Logger object

    Returns:
        bool: 是否成功|Success status
    """
    try:
        cmd = ['bgzip', '-c', input_file]

        if logger:
            logger.info(f"压缩VCF文件|Compressing VCF file: {input_file}")

        with open(output_file, 'wb') as f_out:
            result = subprocess.run(
                cmd,
                stdout=f_out,
                stderr=subprocess.PIPE,
                check=True
            )

        if logger:
            logger.info(f"压缩完成|Compression completed: {output_file}")

        return True

    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"压缩失败|Compression failed: {e.stderr}")
        return False
    except Exception as e:
        if logger:
            logger.error(f"压缩异常|Compression error: {e}")
        return False
