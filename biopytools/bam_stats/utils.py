"""
BAM统计分析工具函数模块|BAM Statistics Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import List


def setup_logger(
    log_dir: str, log_name: str = 'pipeline.log', log_level: str = 'INFO'
) -> logging.Logger:
    """设置日志管理器|Setup logger manager"""
    log_path = Path(log_dir) / '99_logs'
    log_path.mkdir(parents=True, exist_ok=True)
    log_file = log_path / log_name

    logger = logging.getLogger('bam_stats')
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, threads: int = 24):
        self.logger = logger
        self.threads = threads

    def run(
        self,
        cmd: str,
        description: str = "",
        use_threads: bool = True,
    ) -> tuple:
        """
        执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description
            use_threads: 是否注入线程参数|Whether to inject thread params

        Returns:
            (success, stdout_or_stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        if use_threads and self.threads > 1:
            if 'samtools' in cmd and '-@' not in cmd and '--threads' not in cmd:
                for subcmd in ['view', 'sort', 'index', 'depth', 'flagstat']:
                    if f'samtools {subcmd}' in cmd:
                        cmd = cmd.replace(
                            f'samtools {subcmd}',
                            f'samtools {subcmd} -@ {self.threads}',
                        )
                        break

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600,
            )
            return True, result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败|Command execution failed: {description}"
            )
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr[:500]}")
            return False, e.stderr
        except subprocess.TimeoutExpired:
            self.logger.error(
                f"命令超时|Command timed out: {description}"
            )
            return False, 'Timeout'


def get_sample_name(bam_file: str) -> str:
    """从BAM文件路径提取样品名称|Extract sample name from BAM file path"""
    return Path(bam_file).stem.replace('.sorted', '')


def check_dependencies(logger: logging.Logger) -> bool:
    """检查依赖软件是否已安装|Check if required software is installed"""
    missing_libs = []
    for lib in ['pandas', 'openpyxl', 'tqdm']:
        try:
            __import__(lib)
        except ImportError:
            missing_libs.append(lib)

    if missing_libs:
        logger.error(
            f"缺少Python库|Missing Python libraries: "
            f"{' '.join(missing_libs)}"
        )
        return False

    try:
        result = subprocess.run(
            ['samtools', '--version'],
            capture_output=True, text=True, timeout=10,
        )
        version = result.stdout.strip().split('\n')[0]
        logger.info(f"samtools可用|samtools available: {version}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("未找到samtools|samtools not found")
        return False

    return True
