"""
EDTA工具函数模块|EDTA Utility Functions Module
"""

import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

from ..common.conda_runner import build_conda_command


class EDTALogger:
    """EDTA日志管理器(named logger,stdout/stderr 分离 §2.3)|EDTA Logger Manager

    stdout(<=INFO) + stderr(>=WARNING) + 单日志文件;named logger 不污染 root
    |named logger (no root pollution), stdout/stderr separation for HPC .out/.err
    """

    def __init__(self, output_dir: Path, log_name: str = "edta_analysis.log"):
        output_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = output_dir / log_name
        self.logger = self._setup_logging()

    def _setup_logging(self) -> logging.Logger:
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')

        logger = logging.getLogger('biopytools.edta')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        logger.info("EDTA植物基因组TE注释工具启动|EDTA Plant Genome TE Annotation Tool Started")
        logger.info(f"日志文件|Log file: {self.log_file}")
        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner

    命令首元素必须是**完整工具路径**(§13.2.3,禁止 basename 提取);
    经 common 层 build_conda_command 包装(自动 conda run + --no-capture-output)
    |cmd[0] must be a full tool path; wrapped via the common conda runner
    """

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.process = None
        self.start_time = None

    def run(self, cmd: List[str], description: str = "", check: bool = True,
            show_progress: bool = True) -> subprocess.CompletedProcess:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        wrapped_cmd = build_conda_command(cmd[0], list(cmd[1:])) if cmd else cmd

        cmd_str = " ".join(wrapped_cmd) if isinstance(wrapped_cmd, list) else wrapped_cmd
        self.logger.info(f"命令|Command: {cmd_str}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        self.start_time = time.time()

        try:
            # 实时输出(合并 stderr)|Real-time output (stderr merged)
            self.process = subprocess.Popen(
                wrapped_cmd,
                shell=False if isinstance(wrapped_cmd, list) else True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=self.working_dir,
                universal_newlines=True
            )

            output_lines = []
            if show_progress:
                self.logger.info("开始监控进度|Starting progress monitoring...")

            while True:
                output = self.process.stdout.readline()
                if output == '' and self.process.poll() is not None:
                    break
                if output:
                    line = output.strip()
                    output_lines.append(line)
                    if show_progress:
                        self.logger.info(f"EDTA: {line}")

            self.process.wait()

            if check and self.process.returncode != 0:
                raise subprocess.CalledProcessError(
                    self.process.returncode, wrapped_cmd, '\n'.join(output_lines)
                )

            result = subprocess.CompletedProcess(
                wrapped_cmd, self.process.returncode,
                stdout='\n'.join(output_lines), stderr=''
            )

            elapsed_time = time.time() - self.start_time
            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            self.logger.info(f"耗时|Time elapsed: {elapsed_time:.2f} seconds")

            return result

        except subprocess.CalledProcessError as e:
            elapsed_time = time.time() - self.start_time if self.start_time else 0
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"失败前耗时|Time before failure: {elapsed_time:.2f} seconds")
            self.logger.error(f"错误代码|Error code: {e.returncode}")

            if hasattr(e, 'stdout') and e.stdout:
                # 只记尾部,避免整段 EDTA 输出刷爆日志|Tail only, avoid log flooding
                self.logger.error(f"输出尾部|Output tail: {e.stdout[-2000:]}")

            if check:
                raise
            return e

    def terminate_process(self):
        """终止当前进程|Terminate current process"""
        if self.process and self.process.poll() is None:
            self.logger.warning("终止EDTA进程|Terminating EDTA process...")
            self.process.terminate()
            try:
                self.process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self.logger.warning("强制杀死EDTA进程|Force killing EDTA process...")
                self.process.kill()


def check_edta_dependencies(edta_pl: str, logger) -> str:
    """检查 EDTA 可执行性并解析版本(经 conda 包装跑 --help)
    |Check EDTA runs and parse its version (via conda wrapper)

    Args:
        edta_pl: EDTA.pl 完整路径(已解析)|resolved full EDTA.pl path
        logger: 日志器|logger

    Returns:
        EDTA 版本字符串(解析失败为 'unknown')|version string ('unknown' if unparsed)

    EDTA.pl 无参数/坏参数时退出码非 0 亦属正常(--help 非官方选项);
    仅 FileNotFoundError/超时视为不可用。
    |Non-zero exit is acceptable (EDTA.pl has no real --help); only
    FileNotFoundError/timeout mean unavailable.
    """
    logger.info("检查EDTA依赖软件|Checking EDTA dependencies")

    cmd = build_conda_command(edta_pl, ['--help'])
    logger.info(f"命令|Command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, shell=False, capture_output=True,
                                text=True, timeout=120)
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"EDTA检查失败|EDTA check failed: {e}")
        raise RuntimeError(
            f"EDTA 无法运行|EDTA cannot run: {edta_pl}") from e

    text = (result.stdout or '') + (result.stderr or '')
    m = re.search(r'Version\s*:?\s*v?([\d.]+)', text)
    version = m.group(1) if m else 'unknown'
    tail = text.strip().splitlines()
    preview = tail[0] if tail else ''
    logger.info(f"EDTA可用|EDTA available (version={version}; {preview[:60]})")
    return version
