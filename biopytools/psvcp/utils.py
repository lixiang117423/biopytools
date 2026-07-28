"""PSVCP 工具函数模块|PSVCP Utility Functions Module

提供 conda 命令包装、三 handler 日志、命令执行器。
统一走 env psvcp_v.1.0.1:nucmer/assemblytics/bedtools/samtools/python3/Rscript 都在
该 env;vendored helper 内部调用 bedtools / 用 pandas,故 helper 也走 conda run。
|Provides conda wrapping, 3-handler logging, command runners. All tools live in the
psvcp_v.1.0.1 env; vendored helpers internally call bedtools / use pandas, so they too
run via `conda run -n psvcp_v.1.0.1`.
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# --------------------------------------------------------------------------- #
# conda 包装(§13)|conda wrapping
# --------------------------------------------------------------------------- #
def get_conda_env(command: str) -> Optional[str]:
    """检测命令所属 conda 环境,返回环境名|Detect conda env of a command, return env name

    传完整路径(禁止 basename,§13.6.1)。|Pass FULL path (never basename).
    """
    cmd_path = shutil.which(command) or command
    match = re.search(r'/envs/([^/]+)/bin/', cmd_path) or re.search(r'/envs/([^/]+)', cmd_path)
    if match:
        return match.group(1)

    # 兜底:搜索所有 conda 环境|fallback: search all conda envs
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        envs_dir = os.path.join(os.path.dirname(os.path.dirname(conda_exe)), 'envs')
        base_name = os.path.basename(command)
        if os.path.isdir(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, 'bin', base_name)):
                    return env_name
    return None


def build_conda_command(command: str, args) -> List[str]:
    """构建 conda run 命令列表(配合 subprocess.run(shell=False))|Build conda run command list

    必须传完整路径;自动加 --no-capture-output(§13.2.0,避免大数据 OOM)。
    |Must pass full path; auto-adds --no-capture-output (avoids OOM on large data).
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + list(args)
    return [command] + list(args)


# --------------------------------------------------------------------------- #
# 日志(三 handler 分离,§2.3)|logging (3-handler separation)
# --------------------------------------------------------------------------- #
class PSVCPLogger:
    """PSVCP 日志管理器|PSVCP Logger Manager

    stdout(INFO)→.out, stderr(WARNING+)→.err, file(DEBUG+)→psvcp.log
    """

    def __init__(self, log_file: Optional[str] = None, verbose: bool = False):
        self.log_file = log_file
        self.setup_logging(verbose)

    def setup_logging(self, verbose: bool):
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        logger = logging.getLogger("PSVCP")
        logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        logger.handlers.clear()
        logger.propagate = False

        # stdout - INFO|stdout - INFO
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

        # stderr - WARNING+|stderr - WARNING+
        eh = logging.StreamHandler(sys.stderr)
        eh.setLevel(logging.WARNING)
        eh.setFormatter(formatter)
        logger.addHandler(eh)

        # file - 所有级别|file - all levels
        if self.log_file:
            fh = logging.FileHandler(self.log_file, encoding='utf-8')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

        self.logger = logger

    def get_logger(self):
        return self.logger


# --------------------------------------------------------------------------- #
# 命令执行器|command runners
# --------------------------------------------------------------------------- #
def run_command(cmd: List[str], logger: logging.Logger,
                description: str = "", cwd: Optional[str] = None,
                capture: bool = True) -> Tuple[str, str]:
    """执行列表命令(shell=False,§13)|Execute a list command (shell=False)

    Args:
        cmd: 命令列表(由 build_conda_command 构建或裸工具路径)|command list
        logger: 日志器|logger
        description: 步骤描述(先记描述再记完整命令,§2.2.1)|step description
        cwd: 工作目录|working directory
        capture: 是否捕获 stdout(发射器脚本需捕获)|capture stdout?
    Returns:
        (stdout, stderr)
    """
    if description:
        logger.info(f"执行|Executing: {description}")
    logger.info(f"命令|Command: {' '.join(cmd)}")

    result = subprocess.run(
        cmd, shell=False, cwd=cwd,
        capture_output=capture, text=True
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"命令执行失败|Command failed (exit={result.returncode}): "
            f"{(result.stderr or '').strip() or description}"
        )
    return result.stdout or "", result.stderr or ""


def run_shell_command(cmd_str: str, logger: logging.Logger,
                      description: str = "", cwd: Optional[str] = None) -> bool:
    """执行 shell 字符串命令(awk/sort/grep 管道)|Execute a shell-string command

    仅用于 vendored 发射器脚本产出的 awk 管道(含 sort 等),无法用列表表达。
    命令由编排层用绝对路径构建,文件名机器友好无空格(§12.1)。
    |Only for awk pipelines emitted by vendored scripts. Built from absolute paths;
    filenames are machine-friendly (no spaces).
    """
    if description:
        logger.info(f"执行|Executing: {description}")
    logger.info(f"命令|Command: {cmd_str}")
    result = subprocess.run(cmd_str, shell=True, cwd=cwd,
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"shell 命令失败|shell command failed (exit={result.returncode}): "
            f"{result.stderr.strip() or description}"
        )
    return True


def format_number(num: int) -> str:
    """格式化大数字(M 单位,§5.3)|Format large numbers (M unit)"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    if num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def check_assemblytics_numpy(python_path: str, logger: logging.Logger) -> bool:
    """预检 env python 能否导入 numpy(assemblytics 依赖)|Pre-check numpy importable

    env 若被 GraalPy 顶掉(BLOCKER 1),assemblytics 会因 numpy C 扩展失败。
    明确报错引导用户修 env,而非黑箱失败(符合 [[graceful-degradation-preference]])。
    |If the env python is hijacked by GraalPy, assemblytics fails on numpy C-extensions.
    Surface a clear error instead of a black-box failure.
    """
    try:
        cmd = build_conda_command(python_path, ['-c', 'import numpy; print(numpy.__version__)'])
        result = subprocess.run(cmd, shell=False, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            logger.error(
                "env python 无法导入 numpy|env python cannot import numpy.\n"
                "  可能 conda 环境 psvcp_v.1.0.1 的 python 被 GraalPy 顶掉|"
                "the env python may be hijacked by GraalPy.\n"
                "  修复|Fix: conda remove -n psvcp_v.1.0.1 graalpy && "
                "conda install -n psvcp_v.1.0.1 'python=3.10' numpy\n"
                f"  stderr: {(result.stderr or '').strip()}"
            )
            return False
        return True
    except Exception as e:
        logger.error(f"numpy 预检异常|numpy pre-check error: {e}")
        return False
