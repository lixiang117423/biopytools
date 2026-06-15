"""
YaHS工具函数模块|YaHS Utility Functions Module

提供日志管理、命令执行、工具检测等辅助功能
Provides logging, command execution, tool detection and other utilities
"""

import logging
import sys
import subprocess
import shutil
import os
import re
from typing import Optional, List, Tuple
from pathlib import Path


class YaHSLogger:
    """
    YaHS日志管理器|YaHS Logger Manager

    遵循超算日志分离规范：
    - INFO → stdout → .out 文件
    - WARNING+ → stderr → .err 文件
    - 全部 → 本地文件|All → Local file
    """

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志系统|Setup logging system"""
        # 创建logger|Create logger
        self.logger = logging.getLogger("YaHS")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        # → 本地完整日志|→ Local complete log
        if self.log_file:
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        """
        获取日志器|Get logger

        Returns:
            配置好的日志器|Configured logger
        """
        return self.logger


class CommandRunner:
    """
    命令执行器|Command Runner

    执行外部命令，提供日志记录和错误处理
    Execute external commands with logging and error handling
    """

    def __init__(self, logger: logging.Logger, working_dir: Optional[str] = None):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger instance
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = working_dir

    def run_command(
        self,
        cmd: List[str],
        description: str = "",
        check: bool = True,
        capture_output: bool = True,
        timeout: Optional[int] = None,
        use_conda_wrapper: bool = True
    ) -> Tuple[bool, str, str, int]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（不使用shell=True）|Command list (no shell=True)
            description: 步骤描述|Step description
            check: 是否检查返回码|Whether to check return code
            capture_output: 是否捕获输出|Whether to capture output
            timeout: 超时时间(秒)|Timeout in seconds
            use_conda_wrapper: 是否使用conda run包装|Whether to use conda run wrapper

        Returns:
            (成功状态, 标准输出, 标准错误, 返回码)|(Success, stdout, stderr, returncode)
        """
        # 使用conda run包装命令（如果需要）|Use conda run wrapper if needed
        if use_conda_wrapper and len(cmd) > 0:
            cmd = build_conda_command(cmd[0], cmd[1:])

        if description:
            self.logger.info(f"执行|Executing: {description}")
            # 显示完整命令（用于论文写作和调试）|Show full command (for paper writing and debugging)
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
        else:
            # 即使没有描述也显示命令|Show command even without description
            self.logger.info(f"执行|Executing: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                capture_output=capture_output,
                text=True,
                check=False,
                cwd=self.working_dir,
                timeout=timeout
            )

            success = result.returncode == 0

            if not success and check:
                self.logger.error(f"命令执行失败|Command failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Error output: {result.stderr[:500]}")
                return False, result.stdout, result.stderr, result.returncode

            if result.stdout:
                self.logger.debug(f"输出|Output:\n{result.stdout[:200]}")

            return True, result.stdout, result.stderr, result.returncode

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令超时|Command timeout: {description}")
            return False, "", "Timeout expired", -1

        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, "", "Command not found", -1

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {str(e)}")
            return False, "", str(e), -1

    def run_bwa_pipe_samtools(
        self,
        bwa_cmd: List[str],
        samtools_cmd: List[str],
        output_file: str,
        description: str = ""
    ) -> bool:
        """
        运行BWA | samtools管道|Run BWA | samtools pipeline

        自动使用conda run包装命令|Automatically wrap commands with conda run

        Args:
            bwa_cmd: BWA命令列表|BWA command list
            samtools_cmd: samtools命令列表|samtools command list
            output_file: 输出文件|Output file
            description: 步骤描述|Step description

        Returns:
            是否成功|Whether successful
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 使用conda run包装命令|Wrap commands with conda run
        if len(bwa_cmd) > 0:
            bwa_cmd = build_conda_command(bwa_cmd[0], bwa_cmd[1:])
        if len(samtools_cmd) > 0:
            samtools_cmd = build_conda_command(samtools_cmd[0], samtools_cmd[1:])

        # 显示完整管道命令（用于论文写作和调试）|Show full pipeline command
        self.logger.info(f"命令|Command: {' '.join(bwa_cmd)} | {' '.join(samtools_cmd)} > {output_file}")

        try:
            # 启动进程|Start processes
            proc1 = subprocess.Popen(
                bwa_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            proc2 = subprocess.Popen(
                samtools_cmd,
                stdin=proc1.stdout,
                stdout=open(output_file, 'wb'),
                stderr=subprocess.PIPE
            )

            # 关闭管道|Close pipe
            proc1.stdout.close()

            # 等待完成|Wait for completion
            _, err1 = proc1.communicate()
            _, err2 = proc2.communicate()

            if proc2.returncode != 0:
                self.logger.error(f"管道执行失败|Pipeline failed: {description}")
                if err2:
                    self.logger.error(f"错误|Error: {err2.decode()[:500]}")
                return False

            self.logger.info(f"管道执行成功|Pipeline successful: {description}")
            return True

        except Exception as e:
            self.logger.error(f"管道执行异常|Pipeline error: {str(e)}")
            return False


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, return environment name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains envs
        # 例如: /miniforge3/envs/BUSCO_v.6.0.0/bin/busco
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令|Search command in all environments
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件
    Build conda run command to execute software in conda environment

    ⚠️ 重要|IMPORTANT:
        - 直接使用完整路径调用，避免conda run的开销和问题
        - Call full path directly to avoid conda run overhead and issues
        - 必须传递完整路径（以便设置LD_LIBRARY_PATH）|Must pass full path (for LD_LIBRARY_PATH setup)

    Args:
        command: 命令完整路径|Command full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Full command list

    Examples:
        >>> # ✅ 正确：conda环境，直接调用完整路径
        >>> build_conda_command('/miniforge3/envs/BUSCO_v.6.0.0/bin/busco', ['--version'])
        ['/miniforge3/envs/BUSCO_v.6.0.0/bin/busco', '--version']

        >>> # ✅ 正确：非conda环境，使用完整路径
        >>> build_conda_command('/usr/bin/tool', ['--help'])
        ['/usr/bin/tool', '--help']
    """
    # 直接使用完整路径调用，不使用conda run包装
    # 这样避免了：
    # 1. conda run的内存开销
    # 2. conda run的参数解析问题（-n冲突等）
    # 3. --no-capture-output的问题
    # Use full path directly without conda run wrapper to avoid:
    # 1. conda run memory overhead
    # 2. conda run argument parsing issues (e.g., -n conflicts)
    # 3. --no-capture-output issues
    full_cmd = [command] + args

    return full_cmd


def format_number(num: int) -> str:
    """
    格式化数字为大单位显示|Format number to large unit display

    Args:
        num: 数字|Number

    Returns:
        格式化后的字符串|Formatted string

    Examples:
        >>> format_number(10000000)
        '10.00M'
        >>> format_number(1500)
        '1500'
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def calculate_runtime(start_time: float) -> str:
    """
    计算运行时间|Calculate runtime

    Args:
        start_time: 起始时间戳（秒）|Start timestamp (seconds)

    Returns:
        格式化的运行时间|Formatted runtime
    """
    import time
    end_time = time.time()
    duration = int(end_time - start_time)

    hours = duration // 3600
    minutes = (duration % 3600) // 60
    seconds = duration % 60

    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def check_file_exists(file_path: str, min_size: int = 0) -> bool:
    """
    检查文件是否存在且有效|Check if file exists and is valid

    Args:
        file_path: 文件路径|File path
        min_size: 最小文件大小（字节）|Minimum file size (bytes)

    Returns:
        文件是否有效|Whether file is valid
    """
    path = Path(file_path)
    if path.is_file() and path.stat().st_size > min_size:
        return True
    return False


def generate_software_versions_yml(
    config,
    output_file: str,
    logger: logging.Logger
) -> None:
    """
    生成software_versions.yml文件|Generate software_versions.yml file

    Args:
        config: 配置对象|Configuration object
        output_file: 输出文件路径|Output file path
        logger: 日志器|Logger
    """
    import yaml
    from datetime import datetime

    try:
        info = config.get_software_info()
        info['execution'] = {
            'timestamp': datetime.now().isoformat()
        }

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, sort_keys=False)

        logger.info(f"版本信息已保存|Version info saved: {output_file}")

    except Exception as e:
        logger.warning(f"生成版本信息文件失败|Failed to generate version info: {e}")
