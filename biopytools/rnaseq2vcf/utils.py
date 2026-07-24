"""rnaseq2vcf 工具函数|utils (logger, conda wrapping, command runner, checkpoints, sample discovery)"""

import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import threading
from typing import List, Optional, Tuple

LOG_FORMAT = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
LOG_DATEFMT = '%Y-%m-%d %H:%M:%S'


class Rnaseq2vcfLogger:
    """模块日志管理器(三 handler: stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_dir: str, log_file: Optional[str] = None,
                 log_level: str = "INFO", verbose: bool = False, quiet: bool = False):
        os.makedirs(log_dir, exist_ok=True)
        if not log_file:
            log_file = os.path.join(log_dir, "pipeline.log")
        self.log_file = log_file
        self.verbose = verbose
        self.quiet = quiet

        logger = logging.getLogger("rnaseq2vcf")
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG if verbose else getattr(logging, log_level.upper(), logging.INFO))

        formatter = logging.Formatter(LOG_FORMAT, LOG_DATEFMT)

        stdout_h = logging.StreamHandler(sys.stdout)
        stdout_h.setLevel(logging.DEBUG if verbose else logging.INFO)
        stdout_h.setFormatter(formatter)
        logger.addHandler(stdout_h)

        stderr_h = logging.StreamHandler(sys.stderr)
        stderr_h.setLevel(logging.WARNING)
        stderr_h.setFormatter(formatter)
        logger.addHandler(stderr_h)

        file_h = logging.FileHandler(log_file)
        file_h.setLevel(logging.DEBUG)
        file_h.setFormatter(formatter)
        logger.addHandler(file_h)

        self._logger = logger

    def get_logger(self):
        return self._logger

    def step(self, message: str):
        self._logger.info(f"==== {message} ====")


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所属 conda 环境|Detect conda env of a command (by path)"""
    match = re.search(r'/envs/([^/]+)', command)
    if match:
        return match.group(1)
    resolved = shutil.which(command)
    if resolved:
        m2 = re.search(r'/envs/([^/]+)', resolved)
        if m2:
            return m2.group(1)
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建 conda run 命令(带 --no-capture-output)|Build conda run cmd with --no-capture-output"""
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


class CommandRunner:
    """命令执行器|Command runner (list→shell=False; str→shell=True; long→run_with_progress)"""

    def __init__(self, logger: logging.Logger, output_dir: str, dry_run: bool = False):
        self.logger = logger
        self.output_dir = output_dir
        self.dry_run = dry_run

    def run(self, command, description: str = "") -> bool:
        """执行命令(list 或 str)|Execute command (list or str)"""
        cmd_str = command if isinstance(command, str) else ' '.join(command)
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {cmd_str}")
        if self.dry_run:
            return True
        try:
            if isinstance(command, str):
                subprocess.run(command, shell=True, check=True, cwd=self.output_dir)
            else:
                subprocess.run(command, shell=False, check=True, cwd=self.output_dir)
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令失败|Command failed (rc={e.returncode}): {cmd_str}")
            return False
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {e}")
            return False

    def run_with_progress(self, command: str, description: str = "", timeout: int = 86400) -> bool:
        """流式执行长命令(Popen+读线程)|Stream long command via Popen + reader thread"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {command}")
        if self.dry_run:
            return True
        try:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, bufsize=1,
                                       cwd=self.output_dir)
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {e}")
            return False

        stop = threading.Event()

        def _drain():
            for line in iter(process.stdout.readline, ''):
                if stop.is_set():
                    break
                if line.strip():
                    self.logger.info(line.rstrip())
            process.stdout.close()

        reader = threading.Thread(target=_drain, daemon=True)
        reader.start()
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            self.logger.error(f"命令超时|Command timed out after {timeout}s: {command}")
            return False
        finally:
            stop.set()
            reader.join(timeout=5)
        if process.returncode != 0:
            self.logger.error(f"命令失败|Command failed (rc={process.returncode}): {command}")
            return False
        return True


class CheckpointManager:
    """断点续传管理|Checkpoint manager (file-existence based)"""

    def __init__(self, checkpoint_dir: str, logger: logging.Logger):
        self.checkpoint_dir = checkpoint_dir
        self.logger = logger
        os.makedirs(checkpoint_dir, exist_ok=True)

    def _path(self, step: str) -> str:
        return os.path.join(self.checkpoint_dir, f"{step}.done")

    def exists(self, step: str) -> bool:
        return os.path.exists(self._path(step))

    def create(self, step: str):
        try:
            open(self._path(step), 'w').close()
        except Exception as e:
            self.logger.warning(f"写断点失败|Failed to write checkpoint {step}: {e}")

    def remove(self, step: str):
        if self.exists(step):
            os.remove(self._path(step))

    def list_completed(self) -> List[str]:
        return [os.path.basename(p)[:-5] for p in glob.glob(os.path.join(self.checkpoint_dir, "*.done"))]


class FileManager:
    """文件工具|File helpers"""

    @staticmethod
    def ensure_directory(directory: str):
        os.makedirs(directory, exist_ok=True)

    @staticmethod
    def find_files(directory: str, pattern: str) -> List[str]:
        return sorted(glob.glob(os.path.join(directory, pattern)))


class SystemChecker:
    """系统/工具检查|System & tool checks"""

    @staticmethod
    def check_command_exists(command: str, logger: logging.Logger) -> bool:
        if not (shutil.which(command) or os.path.exists(command)):
            logger.error(f"工具缺失|Tool not found: {command}")
            return False
        return True

    @staticmethod
    def check_disk_space(path: str, required_gb: int, logger: logging.Logger) -> bool:
        try:
            usage = shutil.disk_usage(path)
            free_gb = usage.free / (1024 ** 3)
            if free_gb < required_gb:
                logger.warning(f"磁盘空间不足|Low disk: {free_gb:.1f}GB < {required_gb}GB ({path})")
            return free_gb >= required_gb
        except Exception:
            return True


def discover_samples(input_dir: str, read1_pattern: str, read2_pattern: str) -> List[Tuple[str, str, str]]:
    """发现成对样本(R1 glob + R2 确定性替换推导,避免 S1 误匹配 S10)|
    Discover paired samples (R1 glob + deterministic R2, avoid S1 matching S10)"""
    samples: List[Tuple[str, str, str]] = []
    r1_files = sorted(glob.glob(os.path.join(input_dir, f'*{read1_pattern}')))
    for r1 in r1_files:
        base = os.path.basename(r1)
        sample = base[:-len(read1_pattern)] if base.endswith(read1_pattern) else os.path.splitext(base)[0]
        r2 = os.path.join(input_dir, sample + read2_pattern)
        if not os.path.exists(r2):
            raise FileNotFoundError(f"样本|R2 not found for sample {sample}: {r2}")
        samples.append((sample, r1, r2))
    return samples
