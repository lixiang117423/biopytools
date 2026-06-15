"""
基因组组装质量评估工具函数模块|Genome Assembly Quality Control Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Any, Optional, List, Callable
from functools import wraps


def get_conda_env_from_path(conda_env_path: str) -> str:
    """
    从conda环境路径提取环境名称|Extract conda environment name from path

    Args:
        conda_env_path: conda环境路径|conda environment path (e.g., ~/miniforge3/envs/BUSCO_v.6.0.0)

    Returns:
        环境名称|Environment name (e.g., BUSCO_v.6.0.0)
    """
    expanded_path = os.path.expanduser(conda_env_path)
    match = re.search(r'/envs/([^/]+)', expanded_path)
    if match:
        return match.group(1)
    raise ValueError(f"无法从路径提取conda环境名|Cannot extract conda env name from: {conda_env_path}")


def build_conda_command(env_name: str, command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to execute software in conda environment

    Args:
        env_name: conda环境名称|conda environment name (e.g., BUSCO_v.6.0.0)
        command: 命令名称|Command name (e.g., busco)
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表 (适用于subprocess.run(shell=False))|Complete command list (for subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('BUSCO_v.6.0.0', 'busco', ['--version'])
        ['conda', 'run', '-n', 'BUSCO_v.6.0.0', 'busco', '--version']
    """
    full_cmd = ['conda', 'run', '-n', env_name, command] + args
    return full_cmd


class AssemblyQCLogger:
    """基因组组装质量评估日志管理器|Genome Assembly Quality Control Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "assembly_qc.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 创建logger|Create logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        注意|Note:
            - 始终传入列表，使用shell=False（更安全）|Always pass a list, use shell=False (safer)
            - 禁止使用shell=True传入列表，仅执行第一个元素|Never use shell=True with list, only first element executes
            - 若必须使用shell执行字符串，请先: cmd_str = " ".join(cmd)|
              If shell string is needed: cmd_str = " ".join(cmd)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                capture_output=True,
                text=True,
                check=True,
                cwd=str(self.working_dir)
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:500]}")
            return False


def calculate_genome_stats(genome_file: str, logger: logging.Logger) -> Dict[str, Any]:
    """计算基因组统计信息|Calculate genome statistics"""
    logger.info("计算基因组统计信息|Calculating genome statistics")

    stats = {
        'total_size': 0,
        'total_size_mb': 0.0,
        'contig_count': 0,
        'scaffold_count': 0,
        'contig_n50': 0,
        'scaffold_n50': 0,
        'contig_n50_mb': 0.0,
        'scaffold_n50_mb': 0.0,
        'gc_content': 0.0,
        'sequences': []
    }

    try:
        total_length = 0
        lengths = []
        gc_count = 0

        with open(genome_file, 'r') as f:
            current_seq = []
            current_id = ""

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存上一个序列|Save previous sequence
                    if current_seq:
                        seq = ''.join(current_seq)
                        seq_len = len(seq)
                        lengths.append(seq_len)
                        total_length += seq_len
                        gc_count += seq.count('G') + seq.count('C')

                        stats['sequences'].append({
                            'id': current_id,
                            'length': seq_len
                        })

                        if 'contig' in current_id.lower() or 'ctg' in current_id.lower():
                            stats['contig_count'] += 1
                        elif 'scaffold' in current_id.lower() or 'scf' in current_id.lower():
                            stats['scaffold_count'] += 1
                        else:
                            # 未明确标记，默认为contig|Not explicitly marked, default to contig
                            stats['contig_count'] += 1

                    # 开始新序列|Start new sequence
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            # 保存最后一个序列|Save last sequence
            if current_seq:
                seq = ''.join(current_seq)
                seq_len = len(seq)
                lengths.append(seq_len)
                total_length += seq_len
                gc_count += seq.count('G') + seq.count('C')

                stats['sequences'].append({
                    'id': current_id,
                    'length': seq_len
                })

                if 'contig' in current_id.lower() or 'ctg' in current_id.lower():
                    stats['contig_count'] += 1
                elif 'scaffold' in current_id.lower() or 'scf' in current_id.lower():
                    stats['scaffold_count'] += 1
                else:
                    stats['contig_count'] += 1

        stats['total_size'] = total_length
        stats['total_size_mb'] = total_length / 1_000_000
        stats['gc_content'] = (gc_count / total_length * 100) if total_length > 0 else 0

        # 计算N50|Calculate N50
        lengths_sorted = sorted(lengths, reverse=True)
        stats['contig_n50'] = _calculate_n50(lengths_sorted, total_length)
        stats['contig_n50_mb'] = stats['contig_n50'] / 1_000_000

        # 如果有scaffold，计算scaffold N50|Calculate scaffold N50 if available
        if stats['scaffold_count'] > 0:
            scaffold_lengths = [s['length'] for s in stats['sequences']
                              if 'scaffold' in s['id'].lower() or 'scf' in s['id'].lower()]
            if scaffold_lengths:
                scaffold_lengths_sorted = sorted(scaffold_lengths, reverse=True)
                total_scaffold = sum(scaffold_lengths)
                stats['scaffold_n50'] = _calculate_n50(scaffold_lengths_sorted, total_scaffold)
                stats['scaffold_n50_mb'] = stats['scaffold_n50'] / 1_000_000

        logger.info(f"基因组统计完成|Genome statistics completed: {stats['total_size_mb']:.2f} Mb, "
                   f"N50: {stats['contig_n50_mb']:.2f} Mb")

        return stats

    except Exception as e:
        logger.error(f"计算基因组统计失败|Calculate genome statistics failed: {e}")
        return stats


def _calculate_n50(sorted_lengths: list, total_length: int) -> int:
    """计算N50|Calculate N50"""
    target = total_length / 2
    cumulative = 0

    for length in sorted_lengths:
        cumulative += length
        if cumulative >= target:
            return length

    return 0


def parse_busco_summary(summary_file: str) -> Optional[Dict[str, Any]]:
    """解析BUSCO摘要文件|Parse BUSCO summary file"""
    if not os.path.exists(summary_file):
        return None

    try:
        with open(summary_file, 'r') as f:
            for line in f:
                if line.startswith('\tComplete BUSCOs'):
                    parts = line.split()
                    if len(parts) >= 5:
                        return {
                            'complete': float(parts[1].rstrip('%')),
                            'single': float(parts[3].split('(')[0].strip()),
                            'duplicated': float(parts[3].split('(')[1].rstrip('%)')),
                        }
                elif line.startswith('\tFragmented'):
                    parts = line.split()
                    if len(parts) >= 2:
                        return {
                            'fragmented': float(parts[1].rstrip('%'))
                        }
                elif line.startswith('\tMissing'):
                    parts = line.split()
                    if len(parts) >= 2:
                        return {
                            'missing': float(parts[1].rstrip('%'))
                        }

        return None

    except Exception:
        return None


def parse_lai_output(lai_file: str) -> Optional[float]:
    """解析LAI输出文件|Parse LAI output file"""
    if not os.path.exists(lai_file):
        return None

    try:
        with open(lai_file, 'r') as f:
            for line in f:
                if line.startswith('LAI'):
                    parts = line.split()
                    if len(parts) >= 2:
                        return float(parts[1])

        return None

    except Exception:
        return None


def retry_on_io_error(max_retries: int = 3, retry_delay: int = 5,
                      output_files_to_clean: Optional[List[str]] = None):
    """
    重试装饰器：处理NFS I/O错误|Retry decorator for handling NFS I/O errors

    Args:
        max_retries: 最大重试次数|Max retry count
        retry_delay: 重试间隔（秒）|Retry delay in seconds
        output_files_to_clean: 失败时需要清理的输出文件列表|Output files to clean on failure

    Returns:
        装饰器函数|Decorator function
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # self对象通常是第一个参数|self object is usually first arg
            self_obj = args[0] if args else None
            logger = getattr(self_obj, 'logger', None) if self_obj else None

            for attempt in range(1, max_retries + 1):
                try:
                    return func(*args, **kwargs)

                except Exception as e:
                    error_str = str(e)
                    is_io_error = (
                        "Input/output error" in error_str or
                        "error reading" in error_str or
                        "error closing" in error_str
                    )

                    if is_io_error and attempt < max_retries:
                        if logger:
                            logger.warning(
                                f"{func.__name__}遇到I/O错误，{retry_delay}秒后重试 ({attempt}/{max_retries})|"
                                f"{func.__name__} I/O error, retrying in {retry_delay}s ({attempt}/{max_retries})"
                            )
                            logger.warning(f"错误信息|Error: {error_str.strip()}")

                        # 清理部分写入的文件|Clean partially written files
                        if output_files_to_clean and self_obj:
                            for file_path_pattern in output_files_to_clean:
                                # 从kwargs中获取文件路径|Get file path from kwargs
                                if file_path_pattern in kwargs:
                                    file_to_clean = Path(kwargs[file_path_pattern])
                                    if file_to_clean.exists():
                                        try:
                                            file_to_clean.unlink()
                                        except Exception:
                                            pass

                        time.sleep(retry_delay)
                        continue

                    elif is_io_error:
                        if logger:
                            logger.error(
                                f"{func.__name__}失败，已达最大重试次数|"
                                f"{func.__name__} failed after {max_retries} retries"
                            )
                            logger.error(f"最后错误|Final error: {error_str}")
                        raise

                    else:
                        # 非I/O错误，直接抛出|Non-I/O error, raise immediately
                        raise

            return None

        return wrapper

    return decorator
