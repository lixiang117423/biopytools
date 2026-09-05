"""
BAM统计分析工具函数模块|BAM Statistics Analysis Utility Functions Module
"""

import logging
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Union

from ..common.conda_runner import build_conda_command, check_tools


def setup_logger(
    log_dir: str, log_name: str = 'pipeline.log', log_level: str = 'INFO'
) -> logging.Logger:
    """设置日志管理器|Setup logger manager"""
    log_path = Path(log_dir) / '99_logs'
    log_path.mkdir(parents=True, exist_ok=True)
    log_file = log_path / log_name

    logger = logging.getLogger('bam_stats')
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


class CommandRunner:
    """命令执行器|Command Runner

    支持列表命令(conda环境自动包装, shell=False)与管道命令字符串
    (方案B: 管道内直调域环境二进制, 禁止 conda run, shell=True)
    |Supports list commands (auto conda wrap, shell=False) and pipeline
    strings (solution B: direct domain-env binary, no conda run, shell=True)
    """

    THREAD_SUBCMDS = ('view', 'sort', 'index', 'depth', 'flagstat')

    def __init__(self, logger: logging.Logger, threads: int = 24):
        self.logger = logger
        self.threads = threads

    def run(
        self,
        cmd: Union[List[str], str],
        description: str = "",
        use_threads: bool = True,
    ) -> tuple:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表或管道命令字符串|Command list or pipeline string
            description: 步骤描述|Step description
            use_threads: 是否注入线程参数|Whether to inject thread params

        Returns:
            (success, stdout_or_stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        if isinstance(cmd, (list, tuple)):
            if use_threads and self.threads > 1:
                cmd = self._inject_threads_list(list(cmd))
            # 列表命令: conda环境自动包装|List command: auto conda wrap
            full_cmd = build_conda_command(cmd[0], cmd[1:])
            use_shell = False
            display = ' '.join(full_cmd)
        else:
            if use_threads and self.threads > 1:
                cmd = self._inject_threads_string(cmd)
            # 管道字符串: 方案B(§13.2.2), 工具为域环境二进制直调, 不 conda run
            # |Pipeline string: solution B, direct domain binary, no conda run
            full_cmd = cmd
            use_shell = True
            display = cmd

        self.logger.info(f"命令|Command: {display}")

        try:
            result = subprocess.run(
                full_cmd,
                shell=use_shell,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600,
            )
            return True, result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败|Command execution failed: {description}"
            )
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr[:500]}")
            return False, e.stderr
        except subprocess.TimeoutExpired:
            self.logger.error(
                f"命令超时|Command timed out: {description}"
            )
            return False, 'Timeout'

    def _inject_threads_list(self, cmd: List[str]) -> List[str]:
        """列表命令注入线程参数|Inject thread params into list command"""
        if any(a in ('-@', '--threads') for a in cmd):
            return cmd
        # 在已知子命令后插入 -@ N|Insert -@ N after known subcommand
        for i, token in enumerate(cmd):
            if token in self.THREAD_SUBCMDS:
                cmd[i + 1:i + 1] = ['-@', str(self.threads)]
                break
        return cmd

    def _inject_threads_string(self, cmd: str) -> str:
        """管道字符串注入线程参数|Inject thread params into pipeline string"""
        if '-@' in cmd or '--threads' in cmd:
            return cmd
        # 匹配 (路径)samtools 子命令, 在子命令后插入 -@ N
        # |Match (path)samtools subcommand, insert -@ N after subcommand
        pattern = re.compile(
            r'((?:^|\s)[^\s|]*samtools)\s+(view|sort|index|depth|flagstat)(?=\s)'
        )
        m = pattern.search(cmd)
        if m:
            return cmd[:m.end()] + f' -@ {self.threads}' + cmd[m.end():]
        return cmd


def get_sample_name(bam_file: str) -> str:
    """从BAM文件路径提取样品名称|Extract sample name from BAM file path"""
    return Path(bam_file).stem.replace('_sorted', '').replace('.sorted', '')


def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查依赖软件是否已安装|Check if required software is installed"""
    missing_libs = []
    for lib in ['pandas', 'openpyxl', 'tqdm']:
        try:
            __import__(lib)
        except ImportError:
            missing_libs.append(lib)

    if missing_libs:
        logger.error(
            f"缺少Python库|Missing Python libraries: "
            f"{' '.join(missing_libs)}"
        )
        return False

    missing = check_tools(
        [(config.samtools_path, "samtools", ["--version"])], logger
    )
    if missing:
        logger.error("未找到samtools|samtools not found")
        return False

    return True
