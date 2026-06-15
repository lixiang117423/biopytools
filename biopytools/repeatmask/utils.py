"""
重复序列屏蔽工具函数模块|Repeat Masking Utility Functions Module
"""

import logging
import sys
import subprocess
import os
import shutil
import re
from pathlib import Path
from typing import Optional


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


class RepeatMaskLogger:
    """重复序列屏蔽日志管理器|Repeat Masking Logger Manager"""

    def __init__(self, output_path: Path, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_path: 输出目录路径|Output directory path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = output_path / "repeatmask.log"
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, output_path: Path):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器对象|Logger object
            output_path: 输出目录路径|Output directory path
        """
        self.logger = logger
        self.output_path = output_path

    def run(self, command: str, description: str,
            timeout: Optional[int] = None, check: bool = True) -> bool:
        """
        执行命令（自动检测conda环境）|Execute command (auto-detect conda environment)

        Args:
            command: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds
            check: 是否检查返回码|Whether to check return code

        Returns:
            bool: 命令是否成功执行|Whether command executed successfully
        """
        self.logger.info(f"{description}|{description}")

        # 提取命令名称用于检测conda环境|Extract command name for conda environment detection
        cmd_parts = command.strip().split()
        if cmd_parts:
            cmd_name = os.path.basename(cmd_parts[0])

            # 自动检测conda环境|Auto-detect conda environment
            conda_env = get_conda_env(cmd_name)

            if conda_env:
                # 使用conda run包装命令|Use conda run to wrap command
                full_command = f"conda run -n {conda_env} {command}"
                self.logger.debug(f"检测到conda环境|Detected conda environment: {conda_env}")
            else:
                # 直接执行命令|Execute command directly
                full_command = command
                self.logger.debug(f"未检测到conda环境，直接执行|No conda environment detected, executing directly")
        else:
            full_command = command

        self.logger.debug(f"执行命令|Executing command: {full_command}")

        try:
            result = subprocess.run(
                full_command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=check
            )

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

            return True

        except subprocess.TimeoutExpired as e:
            self.logger.error(f"命令执行超时|Command execution timeout: {description}")
            return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            if e.stderr:
                self.logger.error(f"错误输出|Stderr: {e.stderr[:500]}")
            return False

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {e}")
            return False

    def run_with_log(self, command: str, log_file: Path,
                     description: str, timeout: Optional[int] = None) -> bool:
        """
        执行命令并输出到日志文件|Execute command and output to log file

        Args:
            command: 要执行的命令|Command to execute
            log_file: 日志文件路径|Log file path
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            bool: 命令是否成功执行|Whether command executed successfully
        """
        self.logger.info(f"{description}|{description}")
        self.logger.debug(f"执行命令|Executing command: {command}")

        try:
            with open(log_file, 'w') as f:
                process = subprocess.Popen(
                    command,
                    shell=True,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True
                )

                process.wait(timeout=timeout)

                if process.returncode == 0:
                    self.logger.info(f"{description}完成|{description} completed")
                    return True
                else:
                    self.logger.warning(f"{description}返回非零退出码|{description} returned non-zero exit code: {process.returncode}")
                    return True  # 某些工具可能返回非零但实际成功

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command execution timeout: {description}")
            process.kill()
            return False

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {e}")
            return False


def get_genome_stats(genome_file: str, logger) -> dict:
    """
    获取基因组统计信息|Get genome statistics

    Args:
        genome_file: 基因组FASTA文件路径|Genome FASTA file path
        logger: 日志器对象|Logger object

    Returns:
        dict: 基因组统计信息|Genome statistics
    """
    logger.info(f"分析基因组统计信息|Analyzing genome statistics: {genome_file}")

    stats = {
        'sequences': 0,
        'total_length': 0,
        'n50': 0,
        'gc_content': 0.0
    }

    try:
        sequences = []
        total_length = 0
        gc_count = 0

        from Bio import SeqIO
        for record in SeqIO.parse(genome_file, "fasta"):
            sequences.append(len(record.seq))
            total_length += len(record.seq)
            gc_count += record.seq.count('G') + record.seq.count('C')

        if sequences:
            stats['sequences'] = len(sequences)
            stats['total_length'] = total_length
            stats['gc_content'] = (gc_count / total_length * 100) if total_length > 0 else 0

            # 计算N50|Calculate N50
            sequences.sort(reverse=True)
            half_length = total_length / 2
            cumulative = 0
            for seq_len in sequences:
                cumulative += seq_len
                if cumulative >= half_length:
                    stats['n50'] = seq_len
                    break

        logger.info(f"基因组统计|Genome stats: 序列数|Sequences={stats['sequences']}, "
                   f"总长度|Total length={stats['total_length']:,} bp, "
                   f"N50={stats['n50']:,} bp, "
                   f"GC含量|GC content={stats['gc_content']:.2f}%")

    except Exception as e:
        logger.error(f"获取基因组统计信息失败|Failed to get genome statistics: {e}")

    return stats
