"""
rMVP工具函数|rMVP Utility Functions
"""

import os
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional, List


class RMVPLogger:
    """rMVP日志类|rMVP Logger Class"""

    def __init__(self, log_file: Path, log_level: str = "INFO"):
        """
        初始化日志|Initialize logger

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR)
        """
        import sys

        self.log_file = Path(log_file)
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)

        # 确保日志目录存在|Ensure log directory exists
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 创建logger|Create logger
        self.logger = logging.getLogger("RMVP")
        self.logger.setLevel(logging.DEBUG)

        # 清除已有的处理器|Clear existing handlers
        self.logger.handlers = []
        self.logger.propagate = False

        # 格式化器|Formatter
        # 格式: YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 消息中文|Message English
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
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
        file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger对象|Get logger object"""
        return self.logger


def detect_r_environment(r_env_name: Optional[str] = None,
                        r_path: Optional[str] = None,
                        logger: Optional[logging.Logger] = None) -> Tuple[bool, str, Optional[str]]:
    """
    检测R环境|Detect R environment

    Args:
        r_env_name: conda环境名称或路径|conda env name or path
        r_path: R可执行文件路径|R executable path
        logger: 日志对象|Logger object

    Returns:
        (是否成功, conda环境名称或R可执行文件路径, 类型标识)| (success, env_name_or_path, type)
        类型标识|type: "conda" 或 "direct"|"conda" or "direct"
    """
    if logger:
        logger.info(" 检测R环境|Detecting R environment")

    # 处理conda环境|Handle conda environment
    if r_env_name:
        # 检查是否是路径|Check if it's a path
        if '/' in r_env_name or '\\' in r_env_name:
            # 是路径，提取环境名称|It's a path, extract env name
            env_path = Path(r_env_name)
            env_name = env_path.name  # 提取最后一部分作为环境名称|Extract last part as env name
            if logger:
                logger.info(f"   检测到环境路径|Detected env path: {r_env_name}")
                logger.info(f"   提取环境名称|Extracted env name: {env_name}")
            return True, env_name, "conda"
        else:
            # 是环境名称，直接使用|It's an env name, use directly
            if logger:
                logger.info(f"   使用conda环境名称|Using conda env name: {r_env_name}")
            return True, r_env_name, "conda"

    # 2. 如果指定了R路径，直接使用|2. Use specified R path
    if r_path:
        r_path_expanded = Path(r_path).expanduser()
        if r_path_expanded.exists():
            if logger:
                logger.info(f"   使用指定的R路径|Using specified R path: {r_path_expanded}")
            return True, str(r_path_expanded), "direct"
        else:
            if logger:
                logger.warning(f"   指定的R路径不存在|Specified R path does not exist: {r_path}")

    # 3. 在系统PATH中查找R|3. Find R in system PATH
    if logger:
        logger.info("   在系统PATH中查找R|Searching for R in system PATH")
    try:
        which_result = subprocess.run(
            ["which", "R"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if which_result.returncode == 0:
            r_executable = which_result.stdout.strip()
            if logger:
                logger.info(f"   找到R|Found R: {r_executable}")
            return True, r_executable, "direct"
    except Exception as e:
        if logger:
            logger.warning(f"   系统PATH查找失败|System PATH search failed: {e}")

    # 未找到R|R not found
    if logger:
        logger.error(" R未找到|R not found")
        logger.info("   请安装R或指定正确的路径|Please install R or specify the correct path")
        logger.info("   conda安装|Install via conda: conda install -c conda-forge r-base")
    return False, "", None


def check_r_mvp(r_env_name_or_path: str, r_env_type: Optional[str] = None, logger: Optional[logging.Logger] = None) -> bool:
    """
    检查rMVP包是否已安装|Check if rMVP package is installed

    Args:
        r_env_name_or_path: conda环境名称或R可执行文件路径|conda env name or R executable path
        r_env_type: R环境类型标识|"conda" or "direct"
        logger: 日志对象|Logger object

    Returns:
        是否安装了rMVP|Whether rMVP is installed
    """
    if logger:
        logger.info(" 检查rMVP包|Checking rMVP package")

    try:
        # 测试rMVP是否可用|Test if rMVP is available
        if r_env_type == "conda":
            # 使用conda run -n <env_name>|Use conda run -n <env_name>
            # 注意：<env_name> 应该是环境名称（如 rMVP），而不是路径|Note: <env_name> should be env name (e.g., rMVP) not path
            cmd = f"conda run -n {r_env_name_or_path} R --slave -e \"suppressPackageStartupMessages(library(rMVP))\""
        else:
            # 直接调用R|Direct R call
            cmd = f'{r_env_name_or_path} R --slave -e "suppressPackageStartupMessages(library(rMVP))"'

        if logger:
            logger.info(f"   命令|Command: {cmd}")

        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            if logger:
                logger.info("   rMVP包已安装|rMVP package is installed")
            return True
        else:
            if logger:
                logger.error("   rMVP包未安装|rMVP package is not installed")
                logger.info("   安装命令|Install command: conda run -p <env> R -e \"install.packages('rMVP')\"")
                if result.stderr:
                    logger.debug(f"   错误详情|Error details: {result.stderr}")
            return False

    except Exception as e:
        if logger:
            logger.error(f"   rMVP检查失败|rMVP check failed: {e}")
        return False


def prepare_vcf_for_rmvp(vcf_file: str, output_dir: str,
                         logger: Optional[logging.Logger] = None) -> str:
    """
    为rMVP准备VCF文件（解压gzip/bgzip为未压缩格式）|Prepare VCF file for rMVP (decompress gzip/bgzip)

    rMVP的C++ VCF parser不支持任何gzip压缩格式（包括bgzip和普通gzip），必须使用未压缩的VCF文件。
    The C++ VCF parser in rMVP does not support any gzip format (bgzip or plain gzip),
    uncompressed VCF is required.

    Args:
        vcf_file: 原始VCF文件路径|Original VCF file path
        output_dir: 输出目录|Output directory
        logger: 日志对象|Logger object

    Returns:
        可用于rMVP的VCF文件路径|VCF file path usable by rMVP
    """
    vcf_path = Path(vcf_file)

    if not vcf_path.name.lower().endswith('.gz'):
        return vcf_file

    if logger:
        logger.info("   检测到压缩VCF文件，正在解压|Compressed VCF detected, decompressing")

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # .vcf.gz -> .vcf, .gz -> 去掉.gz后缀
    if vcf_path.name.lower().endswith('.vcf.gz'):
        uncompressed_name = vcf_path.name[:-3]  # 去掉 .gz
    else:
        uncompressed_name = vcf_path.name[:-3]  # 去掉 .gz

    uncompressed_vcf = out_path / uncompressed_name

    if uncompressed_vcf.exists():
        if logger:
            logger.info(f"   解压文件已存在，跳过|Decompressed file exists, skipping: {uncompressed_vcf}")
        return str(uncompressed_vcf)

    try:
        result = subprocess.run(
            f"zcat -f '{vcf_file}' > '{uncompressed_vcf}'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=7200
        )
        if result.returncode == 0 and uncompressed_vcf.exists():
            if logger:
                orig_size = vcf_path.stat().st_size / 1024 / 1024
                new_size = uncompressed_vcf.stat().st_size / 1024 / 1024
                logger.info(f"   解压完成|Decompression done: {uncompressed_vcf} ({orig_size:.1f}MB -> {new_size:.1f}MB)")
            return str(uncompressed_vcf)
        else:
            if logger:
                logger.error(f"   VCF解压失败|VCF decompression failed (returncode: {result.returncode})")
    except Exception as e:
        if logger:
            logger.error(f"   VCF解压失败|VCF decompression failed: {e}")

    return vcf_file


def validate_vcf_file(vcf_file: Path, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证VCF文件|Validate VCF file

    Args:
        vcf_file: VCF文件路径|VCF file path
        logger: 日志对象|Logger object

    Returns:
        是否有效|Whether valid
    """
    if not vcf_file.exists():
        if logger:
            logger.error(f"VCF文件不存在|VCF file does not exist: {vcf_file}")
        return False

    # 检查VCF头部|Check VCF header
    try:
        if vcf_file.suffix.lower() == '.gz':
            cmd = f"zcat -f '{vcf_file}' 2>/dev/null | head -20 | grep '##fileformat=VCF'"
        else:
            cmd = f"head -20 '{vcf_file}' | grep '##fileformat=VCF'"

        if logger:
            logger.info(f"   命令|Command: {cmd}")

        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=120
        )

        if result.returncode == 0 and "VCF" in result.stdout:
            if logger:
                logger.info("   VCF文件格式正确|VCF file format is valid")
            return True
        else:
            if logger:
                logger.error("VCF文件格式错误|VCF file format is invalid")
            return False

    except Exception as e:
        if logger:
            logger.error(f"VCF文件验证失败|VCF file validation failed: {e}")
        return False


def validate_phenotype_file(pheno_file: Path, logger: Optional[logging.Logger] = None) -> Tuple[bool, int]:
    """
    验证表型文件并返回列数|Validate phenotype file and return column count

    Args:
        pheno_file: 表型文件路径|Phenotype file path
        logger: 日志对象|Logger object

    Returns:
        (是否有效, 列数)|(Whether valid, column count)
    """
    if not pheno_file.exists():
        if logger:
            logger.error(f"表型文件不存在|Phenotype file does not exist: {pheno_file}")
        return False, 0

    try:
        with open(pheno_file, 'r') as f:
            header = f.readline().strip()
            columns = header.split('\t')
            n_cols = len(columns)

            if n_cols < 2:
                if logger:
                    logger.error(f"表型文件至少需要2列（样本ID + 至少1个表型）|Phenotype file needs at least 2 columns")
                return False, 0

            if logger:
                logger.info(f"   表型文件验证通过|Phenotype file is valid")
                logger.info(f"   检测到{n_cols - 1}个表型|Detected {n_cols - 1} traits: {columns[1:]}")

            return True, n_cols

    except Exception as e:
        if logger:
            logger.error(f"表型文件验证失败|Phenotype file validation failed: {e}")
        return False, 0


class CommandRunner:
    """命令执行器（列表命令，shell=False）|Command Runner (list cmd, shell=False)"""

    def __init__(self, logger, working_dir):
        """
        初始化|Initialize

        Args:
            logger: 日志对象|Logger object
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = str(working_dir)

    def run(self, cmd: List[str], description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        Returns:
            是否成功|Whether successful
        """
        if description:
            self.logger.info(f"   {description}")

        # 记录完整命令到INFO级别|Log complete command at INFO level
        self.logger.info(f"   命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir
            )

            if result.returncode != 0:
                self.logger.error(f"   命令执行失败|Command failed: {description}")
                self.logger.error(f"   错误代码|Error code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"   错误信息|Error message: {result.stderr.strip()}")
                return False

            if result.stdout:
                self.logger.debug(f"   标准输出|Stdout: {result.stdout.strip()}")

            self.logger.info(f"   命令执行成功|Command succeeded: {description}")
            return True

        except Exception as e:
            self.logger.error(f"   命令执行异常|Command error: {description}")
            self.logger.error(f"   异常信息|Exception: {e}")
            return False


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    import shutil
    import re

    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)
    # 命令本身是绝对路径时，直接从路径解析|If command is an absolute path, parse from it
    if '/' in command:
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list

    重要|IMPORTANT:
        必须使用--no-capture-output避免conda缓冲输出导致内存问题
        Must use --no-capture-output to avoid conda buffering output causing memory issues
    """
    # 从路径中提取命令名称|Extract command name from path
    command_name = os.path.basename(command)

    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用，只使用命令名称|Use conda run with command name only
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command_name] + args
    else:
        # 非conda环境，使用完整路径或命令名称|Non-conda environment, use full path or command name
        full_cmd = [command] + args

    return full_cmd

