"""
互作转录组工具函数模块|Dual RNA-seq Utility Functions Module
"""

import os
import re
import shlex
import shutil
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str, preferred: Optional[str] = None) -> Optional[str]:
    """检测命令所在的conda环境名称|Detect conda env name where the command resides"""
    conda_exe = os.environ.get('CONDA_EXE')
    envs_dir = None
    if conda_exe:
        envs_dir = os.path.join(os.path.dirname(os.path.dirname(conda_exe)), 'envs')

    if preferred and envs_dir and os.path.exists(os.path.join(envs_dir, preferred, 'bin', command)):
        return preferred

    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    if envs_dir and os.path.isdir(envs_dir):
        for env_name in os.listdir(envs_dir):
            if os.path.exists(os.path.join(envs_dir, env_name, 'bin', command)):
                return env_name

    return None


# dual_rnaseq 涉及的工具,用于在整条shell命令(含管道)中检测conda环境|
# Tools used by dual_rnaseq, for detecting the conda env within a whole shell command
DUAL_RNASEQ_TOOLS = ['hisat2', 'hisat2-build', 'samtools', 'stringtie']


class DualRNASeqLogger:
    """互作转录组日志管理器|Dual RNA-seq Logger Manager"""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = output_dir / f"dual_rnaseq_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"dual_rnaseq_processing_{timestamp}")

        # 设置日志级别|Set log level
        self.logger.setLevel(logging.DEBUG)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger


class CommandRunner:
    """命令执行器(自动conda包装)|Command runner (auto conda-wrap)"""

    def __init__(self, logger):
        self.logger = logger

    def _conda_wrap(self, cmd: str) -> str:
        """
        把整条shell命令(含管道)包进 conda run -n ENV bash -c '...'|Wrap a whole shell
        command (pipes allowed) in a single conda run activation.

        检测命令中的工具,取第一个能解析出 conda 环境的,整条命令(含 hisat2|samtools 管道)
        在该环境下运行,符合 §13.2.1(避免 conda run | conda run 双重包装)。找不到环境则
        原样执行(依赖 PATH)。|Detects the first tool with a resolvable conda env and runs
        the whole command under it (§13.2.1 compliant). Falls back to PATH if no env.
        """
        for tool in DUAL_RNASEQ_TOOLS:
            if re.search(rf'(^|[\s|;]){re.escape(tool)}\b', cmd):
                env = get_conda_env(tool)
                if env:
                    self.logger.info(f"使用conda环境|Using conda env: {env} (for {tool})")
                    return f"conda run -n {env} --no-capture-output bash -c {shlex.quote(cmd)}"
                return cmd
        return cmd

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令(可为含管道的shell字符串)|Command (shell string, pipes allowed)
            description: 命令描述|Command description
            timeout: 超时时间（秒），None表示无限制|Timeout in seconds, None means no limit

        Returns:
            bool: 执行成功返回True，失败返回False|True if successful, False otherwise
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        # 自动conda包装(含管道)|auto conda-wrap (pipes included)
        full_cmd = self._conda_wrap(cmd)
        if full_cmd != cmd:
            self.logger.info(f"原始命令|Original: {cmd}")

        self.logger.info(f"命令|Command: {full_cmd}")

        if timeout:
            self.logger.info(f"超时设置|Timeout: {timeout}秒|seconds ({timeout/3600:.1f}小时|hours)")

        try:
            result = subprocess.run(
                full_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            self.logger.info(f"{description} 完成|completed")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error(f"{description} 超时|timed out after {timeout}秒|seconds")
            return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class FileValidator:
    """文件验证器|File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在|Check if file exists

        Args:
            file_path: 文件路径|File path
            description: 文件描述|File description

        Returns:
            bool: 文件存在返回True，否则返回False|True if exists, False otherwise
        """
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False
