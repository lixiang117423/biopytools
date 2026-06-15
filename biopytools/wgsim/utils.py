"""
Wgsim工具类|Wgsim Utilities
"""

import logging
import sys
import os
import subprocess
import glob
from typing import List


class WgsimLogger:
    """Wgsim日志管理器|Wgsim Logger Manager"""

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        log_file = os.path.join(output_dir, "wgsim.log")
        self._setup_logging(log_file)

    def _setup_logging(self, log_file: str):
        self.logger = logging.getLogger("Wgsim")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        return self.logger


def get_input_files(input_path: str, extensions: list) -> List[str]:
    """获取输入基因组文件列表|Get input genome file list"""
    if os.path.isfile(input_path):
        return [input_path]

    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(input_path, f"*{ext}")))
    return sorted(set(files))


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    if num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def run_wgsim(
    logger: logging.Logger,
    wgsim_path: str,
    genome_file: str,
    tmp1: str,
    tmp2: str,
    num_reads: int,
    read_length: int,
    seed: int = 0,
    error_rate: float = 0.020,
    mutation_rate: float = 0.001,
    outer_distance: int = 500,
    inner_distance: int = 0,
) -> bool:
    """执行wgsim模拟|Execute wgsim simulation"""
    cmd = [wgsim_path, genome_file, tmp1, tmp2]
    cmd.extend(['-N', str(num_reads)])
    cmd.extend(['-1', str(read_length)])
    cmd.extend(['-2', str(read_length)])
    cmd.extend(['-s', str(seed)])
    cmd.extend(['-e', str(error_rate)])
    cmd.extend(['-r', str(mutation_rate)])
    cmd.extend(['-d', str(outer_distance)])
    cmd.extend(['-D', str(inner_distance)])

    logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        if result.stderr:
            logger.debug(f"wgsim stderr|wgsim stderr: {result.stderr.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"wgsim执行失败|wgsim execution failed: {e.stderr.strip()}")
        return False
    except FileNotFoundError:
        logger.error(f"wgsim不存在|wgsim not found: {wgsim_path}")
        return False


def fix_quality_lines(logger: logging.Logger, file_path: str) -> bool:
    """
    修复wgsim输出的质量行，替换为高Phred值字符

    wgsim根据-e参数生成恒定低质量值（如Q18='2'），下游fastp等工具
    会因默认质量阈值过高而丢弃所有reads。此函数将质量行替换为'I'（Q40），
    同时保留原始序列长度对齐。

    Args:
        logger: 日志对象|Logger object
        file_path: FASTQ文件路径|FASTQ file path

    Returns:
        是否成功|Whether succeeded
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for i, line in enumerate(lines):
            if i % 4 == 3:
                # 质量行：替换为与上一行序列等长的'I'字符（Q40）
                seq_line = new_lines[-1] if new_lines else ''
                seq_len = len(seq_line.rstrip('\n\r'))
                new_lines.append('I' * seq_len + '\n')
            else:
                new_lines.append(line)

        with open(file_path, 'w') as f:
            f.writelines(new_lines)

        logger.debug(f"质量行已修复|Quality lines fixed: {file_path}")
        return True
    except Exception as e:
        logger.error(f"质量行修复失败|Quality line fix failed: {e}")
        return False


def compress_file(logger: logging.Logger, file_path: str) -> bool:
    """使用gzip压缩文件|Compress file with gzip"""
    try:
        cmd = ['gzip', '-f', file_path]
        logger.info(f"命令|Command: {' '.join(cmd)}")
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"压缩失败|Compression failed: {e.stderr.strip()}")
        return False
