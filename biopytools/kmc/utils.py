"""
KMC工具函数模块|KMC Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import shutil
import re
from pathlib import Path
from typing import List, Dict, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command_string(command: str, args: str) -> str:
    """
    构建conda run命令字符串（用于需要shell特性的命令）|Build conda run command string (for commands needing shell features)

    Args:
        command: 命令名称|Command name
        args: 命令参数字符串|Command arguments string

    Returns:
        完整命令字符串|Complete command string
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return f"conda run -n {conda_env} --no-capture-output {command} {args}"
    else:
        return f"{command} {args}"


class KMCLogger:
    """KMC日志管理器|KMC Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "kmc_analysis.log", log_level: str = "INFO"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        level = getattr(logging, log_level.upper(), logging.INFO)

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=level,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令（自动检测conda环境）|Execute command (auto-detect conda environment)"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        # 提取命令名称用于检测conda环境|Extract command name for conda environment detection
        cmd_parts = cmd.strip().split()
        if cmd_parts:
            cmd_name = os.path.basename(cmd_parts[0])

            # 自动检测conda环境|Auto-detect conda environment
            conda_env = get_conda_env(cmd_name)

            if conda_env:
                # 使用conda run包装命令|Use conda run to wrap command
                full_cmd = f"conda run -n {conda_env} --no-capture-output {cmd}"
                self.logger.debug(f"检测到conda环境|Detected conda environment: {conda_env}")
            else:
                # 直接执行命令|Execute command directly
                full_cmd = cmd
                self.logger.debug(f"未检测到conda环境，直接执行|No conda environment detected, executing directly")
        else:
            full_cmd = cmd

        try:
            result = subprocess.run(
                full_cmd,
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

    def run_parallel(self, commands: List[str], description: str = "", max_threads: int = 4) -> Dict[str, bool]:
        """并行执行命令|Execute commands in parallel"""
        if description:
            self.logger.info(f"并行执行|Parallel execution: {description}")

        results = {}

        # 使用简化的并行执行|Use simplified parallel execution
        # 在实际生产中可以使用进程池|In production, can use process pool
        for i, cmd in enumerate(commands):
            self.logger.info(f"执行命令 {i+1}/{len(commands)}|Executing command {i+1}/{len(commands)}")

            success = self.run(cmd, f"{description} - 任务|Task {i+1}")
            results[cmd] = success

        return results


def format_number(num: int) -> str:
    """格式化数字为大单位显示|Format number to large unit display"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def reverse_complement(seq: str) -> str:
    """计算反向互补序列|Calculate reverse complement sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])


def canonical_kmer(seq: str) -> str:
    """获取canonical k-mer（按字典序较小的序列）|Get canonical k-mer (lexicographically smaller sequence)"""
    rc = reverse_complement(seq)
    return min(seq, rc)


def validate_fasta_or_fastq(file_path: str) -> bool:
    """验证FASTA或FASTQ文件格式|Validate FASTA or FASTQ file format"""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()

            if first_line.startswith('>'):
                return True  # FASTA
            elif first_line.startswith('@'):
                return True  # FASTQ
            else:
                return False
    except Exception:
        return False


def get_file_extension(file_path: str) -> str:
    """获取文件扩展名|Get file extension"""
    _, ext = os.path.splitext(file_path)
    return ext.lower()


def is_gzipped(file_path: str) -> bool:
    """检查文件是否为gzip压缩格式|Check if file is gzipped"""
    return file_path.endswith('.gz')


def count_sequences(file_path: str) -> int:
    """统计序列数量|Count sequences"""
    try:
        if is_gzipped(file_path):
            cmd = f"zcat '{file_path}' | wc -l"
        else:
            cmd = f"wc -l '{file_path}'"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = int(result.stdout.strip().split()[0])

        # FASTQ: 4行/序列, FASTA: 2行/序列|FASTQ: 4 lines/seq, FASTA: 2 lines/seq
        if file_path.endswith('.fq') or file_path.endswith('.fq.gz') or \
           file_path.endswith('.fastq') or file_path.endswith('.fastq.gz'):
            return lines // 4
        else:
            return lines // 2

    except Exception as e:
        return 0
