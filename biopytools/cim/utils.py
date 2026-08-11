"""
CIM分析工具函数模块|CIM Analysis Utility Functions Module
"""

import logging
import os
import re
import shutil
import signal
import subprocess
import sys
from typing import List, Optional, Tuple


class CIMLogger:
    """CIM分析日志管理器|CIM Analysis Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers,
            force=True
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger instance
        """
        self.logger = logger
        # 当前子进程(供中断/超时时整组杀掉) | current child proc (for group-kill on interrupt)
        self._child_proc = None
        # 注册SIGTERM处理: 作业调度取消/kill时连带杀子进程组, 避免孤儿conda/R |
        # SIGTERM handler: scheduler-cancel/kill also kills the child group (no orphans)
        try:
            signal.signal(signal.SIGTERM, self._on_terminate_signal)
        except (ValueError, OSError):
            pass  # 非主线程无法注册signal(测试环境) | cannot register outside main thread

    def _on_terminate_signal(self, signum, frame):
        """SIGTERM回调: 杀子进程组后退出|SIGTERM callback: kill child group then exit"""
        self._kill_child_group()
        sys.exit(128 + signum)

    def _kill_child_group(self):
        """杀掉当前子进程所在的整个进程组(sh/conda/Rscript全部)|kill the child's whole
        process group (sh + conda + Rscript, all descendants).

        子进程以start_new_session=True启动, 自成独立进程组; 先SIGTERM优雅退出,
        5秒后SIGKILL兜底。|child runs in its own session via start_new_session=True;
        SIGTERM first, SIGKILL after a 5s grace.
        """
        proc = self._child_proc
        if proc is None:
            return
        try:
            pgid = os.getpgid(proc.pid)
        except ProcessLookupError:
            return  # 进程已退出 | already gone
        self.logger.info(f"终止子进程组|Terminating child process group (pgid={pgid})")
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=5)  # 给优雅退出的时间 | grace period
            return
        except subprocess.TimeoutExpired:
            pass
        try:
            os.killpg(pgid, signal.SIGKILL)  # 兜底强杀 | force kill
        except ProcessLookupError:
            pass

    def run(self, cmd: str, description: str = "",
            timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        子进程及其所有后代(shell/conda/Rscript)以独立进程组启动, 在用户中断(Ctrl-C)、
        超时、或被SIGTERM(作业调度取消)时整组杀掉, 不留孤儿进程。
        |Child + all descendants (sh/conda/Rscript) run in their own process group and
        are killed as a group on Ctrl-C / timeout / SIGTERM — no orphaned processes.

        Args:
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            tuple: (是否成功|success, 标准输出|stdout, 标准错误|stderr)
        """
        if description:
            self.logger.info(f"执行|Running: {description}")
        self.logger.info(f"命令|Command: {cmd}")

        try:
            # start_new_session=True: 子进程树进入独立进程组, 便于整组杀 |
            # child tree in its own session/group so we can kill the whole group
            proc = subprocess.Popen(
                cmd, shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, start_new_session=True,
            )
        except Exception as e:
            self.logger.error(f"命令启动失败|Command launch error: {e}")
            return False, "", str(e)

        self._child_proc = proc
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command timed out ({timeout}s): {cmd}")
            self._kill_child_group()
            try:
                stdout, stderr = proc.communicate(timeout=10)
            except Exception:
                stdout, stderr = "", ""
            self._child_proc = None
            return False, stdout or "", stderr or "Timeout"
        except KeyboardInterrupt:
            # 用户Ctrl-C: 杀掉子进程组(含conda/Rscript), 再向上抛出给main.py处理 |
            # user Ctrl-C: kill child group (incl. conda/Rscript), then propagate up
            self.logger.error("用户中断, 终止子进程组|User interrupted, killing child process group")
            self._kill_child_group()
            try:
                proc.communicate(timeout=10)
            except Exception:
                pass
            self._child_proc = None
            raise
        self._child_proc = None

        if proc.returncode != 0:
            self.logger.error(f"命令执行失败|Command failed (exit {proc.returncode})")
            if stderr:
                # 只记录最后20行错误信息|Only log last 20 lines of stderr
                err_lines = stderr.strip().split('\n')
                for line in err_lines[-20:]:
                    self.logger.error(f"  {line}")
            return False, stdout, stderr

        return True, stdout, stderr

    def run_conda(self, env_name: str, cmd: str, description: str = "",
                  timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """
        在conda环境中执行命令|Execute command in conda environment

        Args:
            env_name: conda环境名|Conda environment name
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            tuple: (是否成功|success, 标准输出|stdout, 标准错误|stderr)
        """
        if '/' in env_name or '~' in env_name:
            full_cmd = f"conda run --no-capture-output -p {env_name} {cmd}"
        else:
            full_cmd = f"conda run --no-capture-output -n {env_name} {cmd}"
        return self.run(full_cmd, description, timeout)


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中,返回环境名称|Detect if command is in a conda environment

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
    """构建conda run命令来运行conda环境中的软件|Build conda run command for conda env software

    Args:
        command: 命令名称或完整路径(建议传完整路径以正确检测env)|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list

    Note:
        必须使用--no-capture-output避免conda缓冲输出导致大数据内存溢出
        |Must use --no-capture-output to avoid conda buffering causing OOM on large data
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args
