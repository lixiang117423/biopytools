"""
RNA-seq分析工具函数模块|RNA-seq Analysis Utility Functions Module
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
    """
    检测命令所在的conda环境名称|Detect conda env name where the command resides

    策略|Strategy:
        1. 优先 preferred 环境(若该环境含此命令)|Prefer user-specified env if it has the command
        2. 从 which 解析出的路径提取 /envs/<name>/|Extract /envs/<name>/ from which path
        3. 搜索所有 conda 环境的 bin 目录|Search all conda envs' bin directories
    """
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


def build_conda_command(command: str, args: List[str], preferred_env: Optional[str] = None) -> List[str]:
    """构建conda run命令(单工具,非管道)|Build conda run command (single tool, not a pipeline)"""
    conda_env = get_conda_env(command, preferred=preferred_env)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


# rnaseq 涉及的工具,用于在整条shell命令(含管道)中检测conda环境|
# Tools used by rnaseq, for detecting the conda env within a whole shell command (incl. pipes)
RNASEQ_TOOLS = ['hisat2', 'hisat2-build', 'samtools', 'stringtie']


class RNASeqLogger:
    """RNA-seq分析日志管理器|RNA-seq Analysis Logger Manager"""

    def __init__(self, output_dir: Path, verbose: bool = False, quiet: bool = False, log_name: str = "rnaseq_analysis.log"):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = output_dir / f"rnaseq_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"rnaseq_processing_{timestamp}")

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

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

    def step(self, message: str):
        """记录步骤|Log step"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


class CommandRunner:
    """命令执行器(自动conda包装+dry_run支持)|Command runner (auto conda-wrap + dry_run)"""

    def __init__(self, logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run

    def _conda_wrap(self, cmd: str) -> str:
        """
        把整条shell命令(支持管道)包进 conda run -n ENV bash -c '...'|Wrap a whole shell
        command (pipes allowed) in a single conda run activation.

        检测命令中出现的 rnaseq 工具,取第一个能解析出 conda 环境的,把整条命令(含管道)
        在该环境下执行。这样 hisat2|samtools 管道两端都在同一环境下运行,符合 §13.2.1
        (避免 conda run | conda run 双重包装)。找不到环境则原样执行(依赖 PATH)。|
        Detects the first rnaseq tool in the command with a resolvable conda env and runs
        the whole command (pipe included) under that env. hisat2|samtools thus share one
        env activation, complying with §13.2.1 (no conda run | conda run double-wrap).
        Falls back to running as-is (PATH) when no env is found.
        """
        for tool in RNASEQ_TOOLS:
            if re.search(rf'(^|[\s|;]){re.escape(tool)}\b', cmd):
                env = get_conda_env(tool)
                if env:
                    wrapped = f"conda run -n {env} --no-capture-output bash -c {shlex.quote(cmd)}"
                    self.logger.info(f"使用conda环境|Using conda env: {env} (for {tool})")
                    return wrapped
                # 找到工具但无conda环境,说明在PATH,原样执行|tool on PATH, run as-is
                return cmd
        return cmd

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令(可为含管道的shell字符串)|Command (shell string, pipes allowed)
            description: 命令描述|Command description
            timeout: 超时时间（秒），None表示无限制|Timeout in seconds, None means no limit
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        # dry_run: 只记录命令不执行|dry_run: log command without executing
        if self.dry_run:
            self.logger.info(f"[DRY RUN] 命令|Command: {cmd}")
            self.logger.info(f"[DRY RUN] 跳过执行|Skipping execution")
            return True

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
            self.logger.error(f"{description} 超时|timed out after {timeout}秒|seconds ({timeout/3600:.1f}小时|hours)")
            self.logger.error(f"跳过该步骤继续处理|Skipping this step and continuing...")
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
        """检查文件是否存在|Check if file exists"""
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False
