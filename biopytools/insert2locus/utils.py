"""insert2locus工具函数|insert2locus utilities"""

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """从完整路径提取conda环境名|Extract conda env name from full path

    必须传完整路径(禁basename提取,会丢/envs/段)|
    Always pass the full path (never basename, which loses the /envs/ segment)
    """
    match = re.search(r"/envs/([^/]+)/", command)
    if match:
        return match.group(1)
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r"/envs/([^/]+)/", cmd_path)
        if match:
            return match.group(1)
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """conda run包装(必须含--no-capture-output)|Wrap with conda run (must include --no-capture-output)"""
    conda_env = get_conda_env(command)
    if conda_env:
        return ["conda", "run", "-n", conda_env, "--no-capture-output", command] + args
    return [command] + args


def format_number(num: int) -> str:
    """大数字M单位|Format big numbers with M unit"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


class Insert2locusLogger:
    """三handler日志管理器(stdout/stderr/file)|Three-handler logger manager"""

    def __init__(self, log_dir: Path, log_name: str, log_level: str = "INFO"):
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = str(log_dir / log_name)  # 必须在setup前赋值|Must set before setup
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """配置命名logger(不污染root)|Configure named logger (no root pollution)"""
        logger = logging.getLogger("insert2locus")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False
        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S")
        # stdout=INFO(超算.out)|stderr=WARNING+(超算.err)|file=DEBUG(本地完整)
        stdout_h = logging.StreamHandler(sys.stdout)
        stdout_h.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        stdout_h.setFormatter(formatter)
        logger.addHandler(stdout_h)
        stderr_h = logging.StreamHandler(sys.stderr)
        stderr_h.setLevel(logging.WARNING)
        stderr_h.setFormatter(formatter)
        logger.addHandler(stderr_h)
        file_h = logging.FileHandler(self.log_file)
        file_h.setLevel(logging.DEBUG)
        file_h.setFormatter(formatter)
        logger.addHandler(file_h)
        self._logger = logger

    def get_logger(self) -> logging.Logger:
        """获取logger|Get logger"""
        return self._logger


class CommandRunner:
    """命令执行器(列表shell=False)|Command runner (list form, shell=False)"""

    def __init__(self, logger: logging.Logger, working_dir: Path):
        self.logger = logger
        self.working_dir = Path(working_dir)

    def run(self, cmd: List[str], description: str = "",
            stdout_file: Optional[Path] = None) -> bool:
        """执行单命令,只关心成败|Run single command, success only"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}"
                         + (f" > {stdout_file}" if stdout_file else ""))
        try:
            if stdout_file is not None:
                Path(stdout_file).parent.mkdir(parents=True, exist_ok=True)
                with open(stdout_file, "wb") as fh:
                    result = subprocess.run(cmd, shell=False,
                                            cwd=str(self.working_dir),
                                            stdout=fh, stderr=subprocess.PIPE)
            else:
                result = subprocess.run(cmd, shell=False,
                                        cwd=str(self.working_dir),
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.PIPE, text=True)
            if result.returncode != 0 and result.stderr:
                err = result.stderr
                if isinstance(err, bytes):
                    err = err.decode(errors="replace")
                self.logger.error(f"命令stderr|Command stderr: {err[:2000]}")
            return result.returncode == 0
        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return False

    def run_capture(self, cmd: List[str], description: str = "") -> Optional[str]:
        """执行单命令并捕获stdout|Run single command and capture stdout"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, cwd=str(self.working_dir),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True)
            if result.returncode != 0:
                self.logger.error(
                    f"命令失败|Command failed (rc={result.returncode}): "
                    f"{result.stderr[:2000]}")
                return None
            return result.stdout
        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return None

    def run_capture_stderr(self, cmd: List[str],
                           description: str = "") -> Optional[str]:
        """执行并捕获stdout+stderr(如spades版本探测)|Capture stdout+stderr combined"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, cwd=str(self.working_dir),
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True)
            if result.returncode != 0:
                return None
            return result.stdout
        except FileNotFoundError:
            return None

    def run_pipeline(self, cmds: List[List[str]], description: str = "") -> bool:
        """Popen直连管道(禁conda run入管道)|Popen-chained pipeline (no conda run in pipe)"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(
            f"命令|Command: {' | '.join(' '.join(c) for c in cmds)}")
        procs = []
        prev = None
        try:
            for i, cmd in enumerate(cmds):
                stdout = subprocess.PIPE if i < len(cmds) - 1 else subprocess.DEVNULL
                procs.append(subprocess.Popen(
                    cmd, shell=False, stdin=prev, stdout=stdout,
                    stderr=subprocess.PIPE, cwd=str(self.working_dir)))
                if prev:
                    prev.close()
                prev = procs[-1].stdout
            ok = True
            for p in procs:
                err = p.stderr.read() if p.stderr else b""
                rc = p.wait()
                if rc != 0:
                    ok = False
                    err_text = err.decode(errors="replace") if isinstance(err, bytes) else err
                    if err_text.strip():
                        self.logger.error(
                            f"管道成员失败(rc={rc})|Pipeline member failed: "
                            f"{err_text[:2000]}")
            return ok
        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return False

    def run_pipeline_capture(self, cmds: List[List[str]],
                             description: str = "") -> Optional[str]:
        """Popen管道并捕获末端stdout(如bwa mem|samtools view收read名)|
        Popen pipeline capturing final stdout

        末端必须stdout=PIPE否则输出漏到终端且捕获到None;
        中间成员stderr=DEVNULL防大量stderr撑爆管道缓冲死锁(bwa mem进度)|
        Final member must use PIPE or output leaks and None is captured;
        intermediate stderr=DEVNULL avoids pipe-buffer deadlock (bwa progress)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(
            f"命令|Command: {' | '.join(' '.join(c) for c in cmds)}")
        procs = []
        prev = None
        try:
            for i, cmd in enumerate(cmds):
                is_last = i == len(cmds) - 1
                procs.append(subprocess.Popen(
                    cmd, shell=False, stdin=prev, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE if is_last else subprocess.DEVNULL,
                    text=True, cwd=str(self.working_dir)))
                if prev:
                    prev.close()
                prev = procs[-1].stdout
            output, err = procs[-1].communicate()
            ok = procs[-1].returncode == 0
            for p in procs[:-1]:
                if p.wait() != 0:
                    ok = False
            if not ok:
                err_text = (err or "").strip()
                self.logger.error(
                    f"管道执行失败|Pipeline failed"
                    + (f": {err_text[:2000]}" if err_text else ""))
                return None
            return output
        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return None
