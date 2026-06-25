"""
CPhasing工具函数模块|CPhasing Utility Functions Module

提供日志管理、命令执行等辅助功能
Provides logging, command execution and other utilities

 重要|IMPORTANT:
    CPhasing 是用 pixi 安装的，不是普通 conda env。用户必须先 source activate_cphasing
    才能使用 biopytools cphasing。本模块不再使用 conda run 包装。
    |CPhasing is installed via pixi, not regular conda env. Users must source
    activate_cphasing before using biopytools cphasing. This module no longer
    wraps with conda run.
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


def check_cphasing_available() -> Tuple[bool, str]:
    """
    检查 cphasing 是否在当前 shell 的 PATH 中可用
    |Check if cphasing is available in the current shell's PATH

    CPhasing 必须通过 `source ~/software/CPhasing_v0.3.0/bin/activate_cphasing`
    激活后才能使用（激活脚本通过 pixi 设置 PATH/PYTHONPATH/LD_LIBRARY_PATH）。
    |CPhasing must be activated via `source .../activate_cphasing` before use
    (the activation script uses pixi to set PATH/PYTHONPATH/LD_LIBRARY_PATH).

    Returns:
        (available, cphasing_path or error_message)
    """
    cphasing_path = shutil.which('cphasing')
    if cphasing_path:
        return True, cphasing_path
    return False, (
        "cphasing 不在 PATH 中|cphasing not in PATH.\n"
        "请先激活 CPhasing 环境|Please activate CPhasing env first:\n"
        "  source ~/software/CPhasing_v0.3.0/bin/activate_cphasing\n"
        "（路径根据你的实际安装位置调整|Adjust path to your actual install）"
    )


def check_cphasing_rs_available() -> Tuple[bool, str]:
    """
    检查 cphasing-rs（Rust 二进制）是否可用
    |Check if cphasing-rs (Rust binary) is available

    CPhasing 内部 18 处调用 cphasing-rs，缺失会导致静默失败。
    |CPhasing internally calls cphasing-rs in 18 places; missing it causes
    silent failures.
    """
    path = shutil.which('cphasing-rs')
    if path:
        return True, path
    return False, (
        "cphasing-rs 不在 PATH 中|cphasing-rs not in PATH.\n"
        "通常意味着 activate_cphasing 没正确执行。"
        "|Usually means activate_cphasing was not sourced correctly."
    )


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
