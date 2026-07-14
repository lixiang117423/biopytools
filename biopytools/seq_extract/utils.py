"""
序列提取工具函数|Sequence extraction utility functions
"""

import logging
import re
import os
import shutil
import subprocess
import sys
from typing import List, Optional, Tuple


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中|Detect if command is in a conda environment

    Args:
        command: 命令名称或路径|Command name or path

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
    """
    构建conda run命令|Build conda run command

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Full command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class SeqExtractLogger:
    """序列提取日志管理器|SeqExtract Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
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
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


class SeqExtractRunner:
    """序列提取执行器|Sequence extraction runner"""

    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def run(self) -> bool:
        """
        执行序列提取|Execute sequence extraction

        Returns:
            bool: 成功/失败|Success/Failure
        """
        cmd = self._build_command()
        desc = self._get_description()

        self.logger.info(f"查询类型|Query type: {self.config.query_type}")
        self.logger.info(f"执行|Executing: {desc}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, check=True)
            self.logger.info(f"输出文件|Output file: {os.path.abspath(self.config.effective_output)}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed (exit code: {e.returncode})")
            return False
        except FileNotFoundError as e:
            self.logger.error(f"seqkit未找到|seqkit not found: {e}")
            return False

    def _build_command(self) -> List[str]:
        """构建seqkit命令|Build seqkit command"""
        if self.config.query_type == "bed_file":
            # seqkit subseq --bed <bed> <sequence> -o <output>
            args = ["subseq", "--bed", self.config.input_query,
                    self.config.sequence_file, "-o", self.config.effective_output]
        elif self.config.query_type == "id_file":
            # seqkit grep -f <id_file> <sequence> -o <output>
            args = ["grep", "-f", self.config.input_query,
                    self.config.sequence_file, "-o", self.config.effective_output]
        else:
            # seqkit grep -p <id> <sequence> -o <output>
            args = ["grep", "-p", self.config.input_query,
                    self.config.sequence_file, "-o", self.config.effective_output]

        return build_conda_command(self.config.seqkit_path, args)

    def _get_description(self) -> str:
        """获取步骤描述|Get step description"""
        descriptions = {
            "single_id": f"单个ID提取序列|Single ID sequence extraction: {self.config.input_query}",
            "id_file": f"ID文件批量提取序列|ID file batch sequence extraction: {self.config.input_query}",
            "bed_file": f"BED区间提取序列|BED region sequence extraction: {self.config.input_query}",
        }
        return descriptions.get(self.config.query_type, "序列提取|Sequence extraction")
