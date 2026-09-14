"""
GenomeSyn2工具函数模块|GenomeSyn2 Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
# build_conda_command_string 依赖 conda_run_prefix 构造shell字符串
from ..common.conda_runner import build_conda_command, conda_run_prefix


def build_conda_command_string(command: str, args: str = "") -> str:
    """
    构建conda run命令字符串（用于需要shell特性的命令）|Build conda run command string (for commands needing shell features)

    Args:
        command: 命令名称|Command name
        args: 命令参数字符串|Command arguments string

    Returns:
        完整命令字符串|Complete command string
    """
    prefix = conda_run_prefix(command)

    if prefix:
        # 使用 conda run 前缀(同源conda绝对路径+run -p,§13)|Use the conda run prefix
        if args:
            full_cmd = f"{prefix} {command} {args}"
        else:
            full_cmd = f"{prefix} {command}"
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
