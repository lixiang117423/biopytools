"""
TGS-GapCloser工具函数模块|TGS-GapCloser Utility Functions Module
"""

import logging
import re
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import List, Optional


class TGSGapCloserLogger:
    """TGS-GapCloser日志管理器(named logger,三handler,§2.3)|Logger (named, 3 handlers)"""

    LOGGER_NAME = "tgsgapcloser"

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
        self.log_file = log_file
        self.log_level = log_level
        self._setup_logging()

    def _setup_logging(self):
        """设置日志系统(named logger,propagate=False,§2.3.3)|Setup (named, no root propagation)"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        logger = logging.getLogger(self.LOGGER_NAME)
        logger.handlers.clear()
        logger.setLevel(logging.DEBUG)
        # §2.3.3:不传播到 root,避免与 biopytools 其他模块串扰/重复输出|don't propagate to root logger
        logger.propagate = False

        # 删除旧日志|Delete old log
        if self.log_file and Path(self.log_file).exists():
            Path(self.log_file).unlink()

        # stdout handler - INFO级别|stdout - INFO
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr - WARNING+
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file handler - DEBUG级别(全量)|file - DEBUG (all)
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所属 conda 环境(按路径)|Detect conda env of a command (by path)"""
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
    """构建 conda run 命令(带 --no-capture-output,§13.2.0)|Build conda run cmd with --no-capture-output"""
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


class CommandRunner:
    """命令执行器(支持 dry_run + 流式)|Command runner (dry_run + streaming)"""

    def __init__(self, logger, output_dir, dry_run: bool = False):
        """初始化命令执行器|Initialize command runner"""
        self.logger = logger
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run

    def run_command(self, command, check=True, timeout=None, cwd=None):
        """
        执行命令(捕获输出)|Execute command (capture output)

        短命令用此方法;长命令(大输出)用 run_with_progress 避免 capture 缓冲 OOM|
        Short cmds here; long cmds → run_with_progress to avoid capture-buffer OOM.
        stderr 降为 DEBUG(工具正常进度常走 stderr,避免误报 WARNING 污染 .err,§2.3)|
        stderr → DEBUG (tools emit progress to stderr; avoid WARNING noise in .err).
        """
        cmd_str = command if isinstance(command, str) else ' '.join(str(c) for c in command)
        self.logger.info(f"命令|Command: {cmd_str}")

        if self.dry_run:
            self.logger.info("dry-run 模式,跳过执行|dry-run, skipping execution")
            return subprocess.CompletedProcess(command, 0)

        try:
            result = subprocess.run(
                command,
                check=check,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=cwd
            )

            if result.stdout:
                for line in result.stdout.splitlines():
                    self.logger.debug(f"stdout| {line}")

            if result.stderr:
                for line in result.stderr.splitlines():
                    self.logger.debug(f"stderr| {line}")

            self.logger.info(f"命令执行完成|Command completed (rc={result.returncode})")
            return result

        except subprocess.TimeoutExpired as e:
            self.logger.error(f"命令执行超时|Command execution timeout: {e}")
            raise

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command failed (rc={e.returncode})")
            if e.stdout:
                for line in e.stdout.splitlines():
                    self.logger.debug(f"stdout| {line}")
            if e.stderr:
                for line in e.stderr.splitlines():
                    self.logger.error(f"stderr| {line}")
            raise

    def run_with_progress(self, command, description: str = "", timeout: int = 86400, cwd=None) -> bool:
        """
        流式执行长命令(Popen+读线程)|Stream long command via Popen + reader thread

        避免 capture_output 缓冲全部输出导致 OOM(§13.2.0);stdout+stderr 合并流式记录到 INFO|
        Avoids capture_output buffering (§13.2.0); merged stream logged to INFO.
        command 为 list→shell=False;为 str→shell=True(仅信任内部拼接的命令)|
        list→shell=False; str→shell=True (only for internally-built commands).
        """
        cmd_str = command if isinstance(command, str) else ' '.join(str(c) for c in command)
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {cmd_str}")

        if self.dry_run:
            self.logger.info("dry-run 模式,跳过执行|dry-run, skipping execution")
            return True

        try:
            process = subprocess.Popen(
                command,
                shell=isinstance(command, str),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=cwd or str(self.output_dir)
            )
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
            self.logger.error(f"命令超时|Command timed out after {timeout}s: {cmd_str}")
            return False
        finally:
            stop.set()
            reader.join(timeout=5)

        if process.returncode != 0:
            self.logger.error(f"命令失败|Command failed (rc={process.returncode}): {cmd_str}")
            return False
        return True

    def build_tgsgapcloser_command(self, config):
        """
        构建TGS-GapCloser命令|Build TGS-GapCloser command

        Args:
            config: TGSGapCloserConfig对象|TGSGapCloserConfig object

        Returns:
            list: 命令列表|Command list
        """
        cmd = [config.tgsgapcloser_path]

        # 必需参数|Required parameters (使用长参数格式|Use long option format)
        cmd.extend(['--scaff', config.scaff_file])
        cmd.extend(['--reads', config.reads_file])
        cmd.extend(['--output', config.output_prefix])

        # TGS类型|TGS type
        cmd.extend(['--tgstype', config.tgstype])

        # 纠错模式|Error correction mode
        if config.mode == 'none':
            cmd.append('--ne')  # no error correction
        elif config.mode == 'racon':
            if config.racon_path:
                cmd.extend(['--racon', config.racon_path])
        elif config.mode == 'pilon':
            if config.ngs_file:
                cmd.extend(['--ngs', config.ngs_file])
            if config.pilon_path:
                cmd.extend(['--pilon', config.pilon_path])
            if config.samtools_path:
                cmd.extend(['--samtools', config.samtools_path])
            if config.java_path:
                cmd.extend(['--java', config.java_path])

        # 过滤参数|Filter parameters
        if config.min_idy is not None:
            cmd.extend(['--min_idy', str(config.min_idy)])
        if config.min_match is not None:
            cmd.extend(['--min_match', str(config.min_match)])

        # 线程数|Threads
        cmd.extend(['--thread', str(config.threads)])

        # 分块数量|Chunk count
        cmd.extend(['--chunk', str(config.chunk)])

        # Gap大小差异检查|Gap size difference check
        if config.g_check:
            cmd.append('--g_check')

        # 最小/最大reads数量|Min/Max read count
        cmd.extend(['--min_nread', str(config.min_nread)])
        cmd.extend(['--max_nread', str(config.max_nread)])

        # 最大候选数|Max candidates
        cmd.extend(['--max_candidate', str(config.max_candidate)])

        # Racon参数|Racon parameters
        if config.mode == 'racon':
            cmd.extend(['--r_round', str(config.racon_round)])

        # Pilon参数|Pilon parameters
        if config.mode == 'pilon':
            cmd.extend(['--pilon_mem', config.pilon_mem])
            cmd.extend(['--p_round', str(config.pilon_round)])

        # 自定义minimap2参数|Custom minimap2 parameters
        if config.minmap_arg:
            cmd.extend(['--minmap_arg', config.minmap_arg])

        return cmd
