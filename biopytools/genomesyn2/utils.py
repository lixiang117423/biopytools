"""
GenomeSyn2工具函数模块|GenomeSyn2 Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import shutil
import re
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'perl' or '/path/to/perl')

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        # 例如: /miniforge3/envs/genomesyn2/bin/perl
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    # 尝试找到conda基础目录|Try to find conda base directory
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        # 需要获取envs目录|Need to get envs directory
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令|Search for command in all environments
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    # 找到了|Found it
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    # 检查是否在conda环境中|Check if in conda environment
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用 conda run|Use conda run
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 直接调用|Direct call
        full_cmd = [command] + args

    return full_cmd


def build_conda_command_string(command: str, args: str = "") -> str:
    """
    构建conda run命令字符串（用于需要shell特性的命令）|Build conda run command string (for commands needing shell features)

    Args:
        command: 命令名称|Command name
        args: 命令参数字符串|Command arguments string

    Returns:
        完整命令字符串|Complete command string
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用 conda run|Use conda run
        if args:
            full_cmd = f"conda run -n {conda_env} --no-capture-output {command} {args}"
        else:
            full_cmd = f"conda run -n {conda_env} --no-capture-output {command}"
    else:
        # 直接调用|Direct call
        if args:
            full_cmd = f"{command} {args}"
        else:
            full_cmd = command

    return full_cmd


class GenomeSyn2Logger:
    """GenomeSyn2日志管理器|GenomeSyn2 Logger Manager"""

    def __init__(self, output_dir: Optional[Path] = None, log_name: str = "genomesyn2.log"):
        self.output_dir = output_dir
        self.log_name = log_name

        if output_dir:
            self.log_file = output_dir / log_name
        else:
            self.log_file = None

        self.setup_logging()

    def setup_logging(self, log_level: str = "INFO"):
        """设置日志|Setup logging"""
        # 删除已存在的日志文件|Delete existing log file
        if self.log_file and self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class GenomeSyn2CommandRunner:
    """GenomeSyn2命令执行器|GenomeSyn2 Command Runner"""

    def __init__(self, logger, working_dir: Optional[Path] = None):
        self.logger = logger
        self.working_dir = working_dir or Path.cwd()

    def run_perl_script(self, script_path: str, args: list, description: str = "") -> bool:
        """执行Perl脚本|Execute Perl script"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 构建命令，使用build_conda_command包装|Build command, wrap with build_conda_command
        cmd = build_conda_command(script_path, args)
        cmd_str = " ".join(cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                shell=False,  # 使用列表形式时必须使用shell=False|Must use shell=False with list
                cwd=self.working_dir
            )

            # 记录输出|Log output
            if result.stdout:
                # 只记录重要信息|Log only important messages
                for line in result.stdout.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        self.logger.debug(f"STDOUT: {line}")

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error message: {result.stderr}")
                return False

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True

        except FileNotFoundError as e:
            self.logger.error(f"文件未找到|File not found: {e}")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {e}")
            return False

    def run_perl_script_via_shell(self, perl_path: str, script_path: str,
                                  args: list, description: str = "") -> bool:
        """通过shell执行Perl脚本|Execute Perl script via shell"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 构建命令字符串，使用build_conda_command_string包装|Build command string, wrap with build_conda_command_string
        args_str = " ".join(args)
        cmd_str = build_conda_command_string(perl_path, f"{script_path} {args_str}")
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd_str,
                shell=True,
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir
            )

            # 记录输出|Log output
            if result.stdout:
                for line in result.stdout.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        self.logger.debug(f"STDOUT: {line}")

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error message: {result.stderr}")
                return False

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {e}")
            return False


def format_genome_size(size_bp: int) -> str:
    """格式化基因组大小|Format genome size"""
    if size_bp >= 1_000_000_000:
        return f"{size_bp / 1_000_000_000:.2f}Gb"
    elif size_bp >= 1_000_000:
        return f"{size_bp / 1_000_000:.2f}Mb"
    elif size_bp >= 1_000:
        return f"{size_bp / 1_000:.2f}Kb"
    else:
        return f"{size_bp}bp"


def check_perl_module(perl_path: str, module_name: str) -> bool:
    """检查Perl模块是否可用|Check if Perl module is available"""
    try:
        args = ['-M' + module_name, '-e', 'print "OK\n"']
        cmd = build_conda_command(perl_path, args)
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, shell=False)
        return result.returncode == 0 and "OK" in result.stdout
    except Exception:
        return False
