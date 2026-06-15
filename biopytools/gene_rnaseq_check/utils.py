"""候选基因RNA-seq转录验证工具函数模块|Candidate Gene RNA-seq Validation Utility Functions Module"""

import os
import re
import subprocess
import sys
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple


# ============================================================
# conda 环境检测与命令构建|Conda environment detection and command building
# ============================================================

def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在conda环境|Detect conda environment for a command"""
    import shutil

    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_base = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base, 'envs')
        if os.path.isdir(envs_dir):
            cmd_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', cmd_name)
                if os.path.exists(env_bin):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建conda run命令|Build conda run command

    Args:
        command: 命令完整路径|Full command path
        args: 参数列表|Argument list

    Returns:
        命令列表（适用于 subprocess.run(shell=False)）|Command list for subprocess.run(shell=False)
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


# ============================================================
# 日志管理器|Logger Manager
# ============================================================

class GeneRnaseqCheckLogger:
    """候选基因RNA-seq转录验证日志管理器|Gene RNA-seq Check Logger Manager"""

    def __init__(self, log_dir: Path, verbose: bool = False, quiet: bool = False):
        self.log_dir = log_dir
        self.log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = log_dir / f"gene_rnaseq_check_{timestamp}.log"

        self.logger = logging.getLogger(f"gene_rnaseq_check_{timestamp}")
        self.logger.propagate = False

        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
        )

        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger实例|Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤分隔|Log step separator"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


# ============================================================
# 命令执行器|Command Runner
# ============================================================

class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def run(self, cmd: List[str], description: str = "", cwd: str = None,
            timeout: int = 86400) -> Tuple[bool, str, str]:
        """执行命令|Execute command

        Args:
            cmd: 命令列表|Command list (from build_conda_command)
            description: 命令描述|Command description
            cwd: 工作目录|Working directory
            timeout: 超时秒数|Timeout in seconds

        Returns:
            (success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        if timeout and timeout < 3600:
            pass
        elif timeout:
            self.logger.info(
                f"超时设置|Timeout: {timeout}秒|seconds ({timeout / 3600:.1f}小时|hours)"
            )

        try:
            result = subprocess.run(
                cmd, shell=False, capture_output=True, text=True,
                check=True, cwd=cwd, timeout=timeout,
            )
            self.logger.info(f"{description} 完成|completed")
            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误输出|Error output: {e.stderr[:500]}")
            return False, e.stdout, e.stderr

        except subprocess.TimeoutExpired:
            self.logger.error(
                f"{description} 超时|timed out after {timeout}秒|seconds "
                f"({timeout / 3600:.1f}小时|hours)"
            )
            return False, '', 'timeout'

        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, '', f'Command not found: {cmd[0]}'


# ============================================================
# 辅助函数|Utility Functions
# ============================================================

def format_number(num: int) -> str:
    """格式化大数字|Format large numbers"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def generate_software_versions(output_dir: str, tools: Dict[str, str],
                               params: Dict, start_time=None):
    """生成software_versions.yml|Generate software_versions.yml"""
    try:
        import yaml
        from datetime import datetime
    except ImportError:
        return

    if start_time is None:
        start_time = datetime.now()

    versions = {}
    for name, path in tools.items():
        try:
            cmd = build_conda_command(path, ['--version'])
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            version = r.stdout.strip().split('\n')[0] if r.returncode == 0 else 'unknown'
        except Exception:
            version = 'unknown'
        versions[name] = {'version': version, 'path': path}

    end_time = datetime.now()
    info = {
        'pipeline': {
            'name': 'biopytools gene_rnaseq_check',
            'version': '1.0.0',
        },
        'tools': versions,
        'parameters': params,
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': int((end_time - start_time).total_seconds()),
        },
    }

    info_dir = Path(output_dir) / '00_pipeline_info'
    info_dir.mkdir(parents=True, exist_ok=True)
    info_file = info_dir / 'software_versions.yml'
    with open(info_file, 'w', encoding='utf-8') as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True)


def get_chrom_lengths(bam_file: str) -> Dict[str, int]:
    """从BAM header获取染色体长度|Get chromosome lengths from BAM header"""
    lengths = {}
    try:
        import pysam
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for ref_info in bam.header.references:
                ref_name = bam.header.references[bam.header.references.index(ref_info)]
                lengths[ref_name] = bam.header.get_reference_length(ref_name)
    except Exception:
        pass
    return lengths
