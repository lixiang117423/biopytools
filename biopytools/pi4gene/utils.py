"""
Pi4Gene工具函数模块|Pi4Gene Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


class Pi4GeneLogger:
    """Pi4Gene日志管理器|Pi4Gene Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "pi4gene_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / '99_logs' / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        logger = logging.getLogger("Pi4GeneAnalyzer")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, return environment name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    if os.path.exists(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件
    Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        cmd_name = Path(command).name
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', cmd_name] + args
    else:
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "", timeout: int = 3600, stdout_file: str = None) -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            stdout_fh = open(stdout_file, 'w') if stdout_file else None
            try:
                result = subprocess.run(
                    cmd,
                    shell=False,
                    stdout=stdout_fh,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    timeout=timeout,
                    cwd=self.working_dir
                )

                self.logger.info(f"命令执行成功|Command executed successfully: {description}")

                if not stdout_file and result.stdout:
                    self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

                return True
            finally:
                if stdout_fh:
                    stdout_fh.close()

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:5000]}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command execution timeout: {description}")
            return False


def parse_id_file(id_file: str) -> Dict[str, List[str]]:
    """
    解析分组ID文件，返回分组到序列ID列表的映射
    Parse group ID file, return mapping of group to sequence ID list

    Args:
        id_file: 分组ID文件路径|Group ID file path (column 1: group, column 2: seq_id)

    Returns:
        {分组名: [序列ID列表]}|{group_name: [seq_id_list]}
    """
    group_to_ids: Dict[str, List[str]] = {}

    with open(id_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if '\t' in line:
                parts = line.split('\t')
            elif ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()

            if len(parts) >= 2:
                group_name = parts[0].strip()
                seq_id = parts[1].strip()
                if group_name not in group_to_ids:
                    group_to_ids[group_name] = []
                group_to_ids[group_name].append(seq_id)

    return group_to_ids


def get_software_version(tool_path: str, logger) -> str:
    """自动检测软件版本|Auto-detect software version"""
    try:
        cmd = build_conda_command(tool_path, ['--version'])
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            return version
        else:
            return "unknown"
    except Exception as e:
        logger.warning(f"版本检测失败|Version detection failed: {e}")
        return "unknown"


def format_number(num: int) -> str:
    """格式化大数字|Format large number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return f"{num:,}"
