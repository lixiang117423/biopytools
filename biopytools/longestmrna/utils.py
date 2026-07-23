"""
最长转录本提取工具函数模块|Longest mRNA Extraction Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import List, Optional, Tuple


class LongestMRNALogger:
    """最长转录本提取日志管理器|Longest mRNA Extraction Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger("longestmrna")
        logger.setLevel(level)
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
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda env, return env name

    Args:
        command: 命令名称或路径|Command name or path
    Returns:
        conda环境名称或None|Conda env name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

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


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command for conda env software

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list
    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        self.logger = logger

    def run(self, cmd: list, description: str = "", check: bool = True) -> subprocess.CompletedProcess:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description
            check: 是否检查返回码|Whether to check return code
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                encoding='utf-8',
                check=check
            )

            if result.returncode == 0:
                self.logger.info(f"命令执行成功|Command executed successfully: {description}")
                if result.stdout.strip():
                    self.logger.debug(f"标准输出|Stdout: {result.stdout.strip()}")
            else:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                if result.stderr.strip():
                    self.logger.error(f"错误信息|Error: {result.stderr.strip()}")
                if check:
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)

            return result

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            raise


class TempFileManager:
    """临时文件管理器|Temporary File Manager

    可选 base_dir 把临时文件重定向到 output/tmp(避免超算系统 /tmp 爆满)|
    Optional base_dir redirects temp files to output/tmp (avoids filling the
    supercomputer's system /tmp). base_dir=None 时回退 tempfile 默认目录(向后兼容)|
    When base_dir is None, falls back to tempfile's default dir (backward compatible).
    """

    def __init__(self, logger, base_dir: Optional[str] = None):
        self.logger = logger
        # None 时回退系统默认(向后兼容)|None -> system default (backward compatible)
        self.base_dir = base_dir
        self.temp_files = []
        if base_dir:
            # 存在则 no-op,不存在则创建(支持重复 manager)|no-op if exists, create otherwise
            os.makedirs(base_dir, exist_ok=True)

    def create_temp_file(self, mode='w+', delete=False, suffix='', encoding='utf-8'):
        """创建临时文件|Create temporary file

        dir=self.base_dir(None 时 tempfile 用系统默认)|
        dir=self.base_dir (None -> tempfile uses system default)
        """
        temp_file = tempfile.NamedTemporaryFile(
            mode=mode,
            delete=delete,
            suffix=suffix,
            encoding=encoding,
            dir=self.base_dir
        )
        self.temp_files.append(temp_file.name)
        return temp_file

    def cleanup(self):
        """清理临时文件|Cleanup temporary files"""
        for temp_file in self.temp_files:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    self.logger.debug(f"清理临时文件|Cleaned up temp file: {temp_file}")
            except Exception as e:
                self.logger.warning(f"清理临时文件失败|Failed to cleanup temp file {temp_file}: {e}")
        self.temp_files.clear()  # 清空列表,避免复用时重复删除|Clear to avoid redundant deletes on reuse
