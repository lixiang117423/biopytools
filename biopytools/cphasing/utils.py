"""
CPhasing工具函数模块|CPhasing Utility Functions Module

提供日志管理、命令执行、conda环境包装等辅助功能
Provides logging, command execution, conda environment wrapping and other utilities
"""

import os
import re
import logging
import sys
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional

from ..common.paths import get_tool_path, expand_path


# CPhasing软件目录（包含bin下的二进制工具）
# CPhasing software directory (contains binary tools in bin/)
CPHASING_DIR = get_tool_path(
    'cphasing',
    '~/software/CPhasing/CPhasing-main',
    'CPHASING_DIR'
)


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, returns environment name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|Conda environment name or None
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
    构建conda run命令来运行conda环境中的软件
    Build conda run command to execute software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Full command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def get_cphasing_env() -> dict:
    """
    获取CPhasing运行所需的环境变量
    Get environment variables needed for CPhasing execution

    清理所有可能指向旧版CPhasing的环境变量和路径，
    确保conda run加载正确的版本
    Clean all env vars/paths pointing to old CPhasing to ensure correct version
    """
    env = os.environ.copy()

    # 移除CPHASING_DIR环境变量，防止指向旧版CPhasing
    # Remove CPHASING_DIR to prevent pointing to old CPhasing
    env.pop('CPHASING_DIR', None)

    # 清理PYTHONPATH中的旧CPhasing路径
    # Clean old CPhasing paths from PYTHONPATH
    pythonpath = env.get('PYTHONPATH', '')
    if pythonpath:
        cleaned = ':'.join(
            p for p in pythonpath.split(':')
            if p and 'cphasing' not in p.lower() and 'CPhasing' not in p
        )
        env['PYTHONPATH'] = cleaned if cleaned else ''

    # 清理PATH中的旧CPhasing bin路径
    # Clean old CPhasing bin paths from PATH
    path = env.get('PATH', '')
    cleaned_path = ':'.join(
        p for p in path.split(':')
        if p and not (
            '/cphasing/CPhasing' in p or
            '/Cphasing/CPhasing' in p
        )
    )
    env['PATH'] = cleaned_path

    # 添加新版CPhasing的bin目录到PATH（包含cphasing-rs, allhic等二进制工具）
    # 不使用os.path.isdir检查，NFS缓存可能导致误判
    # Add new CPhasing bin to PATH (contains cphasing-rs, allhic, etc.)
    # Skip os.path.isdir check due to NFS cache issues
    cphasing_bin = os.path.join(expand_path('~/software/CPhasing/CPhasing-main'), 'bin')
    env['PATH'] = f"{cphasing_bin}:{env['PATH']}"

    return env


class CPhasingLogger:
    """
    CPhasing日志管理器|CPhasing Logger Manager

    遵循超算日志分离规范：
    - INFO -> stdout -> .out 文件
    - WARNING+ -> stderr -> .err 文件
    - 全部 -> 本地文件|All -> Local file
    """

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志系统|Setup logging system"""
        self.logger = logging.getLogger("CPhasing")
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
        self.logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        if self.log_file:
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    # 即便退出码为0，stderr里出现这些模式也认为是失败
    # |Even with exit code 0, these stderr patterns indicate failure
    # 用于防御上游工具（如CPhasing）吞异常的bug
    # |Defends against upstream tools (e.g. CPhasing) that swallow exceptions
    SILENT_FAILURE_PATTERNS = (
        'Traceback (most recent call last)',
        'AssertionError',
        'AttributeError',
        'RuntimeError',
        'ValueError: ',
        'KeyError: ',
        'IndexError: ',
        'FileNotFoundError',
        'PermissionError',
    )

    def __init__(self, logger: logging.Logger, working_dir: Optional[str] = None):
        self.logger = logger
        self.working_dir = working_dir

    def _detect_silent_failure(self, stderr: str) -> Optional[str]:
        """扫描stderr检测静默失败|Scan stderr for silent failure patterns"""
        if not stderr:
            return None
        for pattern in self.SILENT_FAILURE_PATTERNS:
            if pattern in stderr:
                return pattern
        return None

    def run_command(
        self,
        cmd: List[str],
        description: str = "",
        cwd: Optional[str] = None,
        timeout: Optional[int] = None,
        extra_env: Optional[dict] = None,
    ) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description
            cwd: 工作目录|Working directory
            timeout: 超时时间(秒)|Timeout in seconds
            extra_env: 额外环境变量|Extra environment variables
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        env = None
        if extra_env:
            env = os.environ.copy()
            env.update(extra_env)

        work_dir = cwd or self.working_dir

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=work_dir,
                timeout=timeout,
                env=env,
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed (rc={result.returncode}): {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Error output: {result.stderr[:500]}")
                return False, result.stdout, result.stderr

            # 兜底：即便退出码为0，stderr里检测异常模式
            # |Safety net: detect exception patterns in stderr even with exit code 0
            silent_pattern = self._detect_silent_failure(result.stderr or "")
            if silent_pattern:
                self.logger.error(
                    f"命令静默失败|Silent failure detected: {description} "
                    f"(rc=0 但stderr含|but stderr contains '{silent_pattern}')"
                )
                # 打印完整stderr方便排查|Print full stderr for debugging
                if result.stderr:
                    self.logger.error(f"完整stderr|Full stderr:\n{result.stderr[-2000:]}")
                return False, result.stdout, result.stderr

            if result.stdout:
                self.logger.debug(f"输出|Output:\n{result.stdout[:200]}")

            return True, result.stdout, result.stderr

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令超时|Command timeout: {description}")
            return False, "", "Timeout expired"

        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, "", f"Command not found: {cmd[0]}"

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {str(e)}")
            return False, "", str(e)
