"""
TeloComp工具函数模块|TeloComp Utility Functions Module
"""

import logging
import sys
import subprocess
import os


class TeloCompLogger:
    """TeloComp日志管理器|TeloComp Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
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
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CondaCommandRunner:
    """Conda环境命令执行器|Conda Environment Command Runner"""

    def __init__(self, logger, conda_env):
        """初始化命令执行器|Initialize command runner"""
        self.logger = logger
        self.conda_env = conda_env
        self.conda_python = f"{conda_env}/bin/python"

    def run_telocomp_command(self, command_name, args, check=True):
        """
        在conda环境中运行TeloComp命令|Run TeloComp command in conda environment

        Args:
            command_name: 命令名称|Command name (e.g., 'telocomp_Filter_1')
            args: 参数列表|Argument list
            check: 是否检查返回码|Whether to check return code

        Returns:
            subprocess.CompletedProcess: 命令执行结果|Command execution result
        """
        # 构建完整命令|Build full command
        full_command = [self.conda_python, f"{os.environ.get('TELOCOMP_BIN', '')}/{command_name}"]
        full_command.extend(args)

        self.logger.info(f"执行命令|Executing command: {' '.join(command_name.split())} {' '.join(args)}")

        try:
            result = subprocess.run(
                full_command,
                check=check,
                capture_output=True,
                text=True,
                env=self._get_conda_env()
            )

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            if result.stderr and result.returncode != 0:
                self.logger.warning(f"标准错误|Stderr: {result.stderr}")

            return result

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            raise

    def _get_conda_env(self):
        """获取conda环境变量|Get conda environment variables"""
        env = os.environ.copy()

        # 添加TeloComp bin目录到PATH|Add TeloComp bin directory to PATH
        telocomp_bin = os.environ.get('TELOCOMP_BIN', '')
        genomesyn_bin = os.environ.get('GENOMESYN_BIN', '')

        path = env.get('PATH', '')
        if telocomp_bin:
            path = f"{telocomp_bin}:{path}"
        if genomesyn_bin:
            path = f"{genomesyn_bin}:{path}"

        env['PATH'] = path
        env['PYTHONPATH'] = f"{self.conda_env}/lib/python3.12/site-packages"

        return env


def check_genome_index(genome_file):
    """
    检查基因组索引文件是否存在|Check if genome index file exists

    Args:
        genome_file: 基因组FASTA文件路径|Genome FASTA file path

    Returns:
        bool: 如果索引存在返回True，否则返回False
               Returns True if index exists, False otherwise
    """
    fai_file = f"{genome_file}.fai"
    return os.path.exists(fai_file)


def create_genome_index(genome_file, logger, conda_env):
    """
    创建基因组索引文件|Create genome index file

    Args:
        genome_file: 基因组FASTA文件路径|Genome FASTA file path
        logger: 日志器|Logger
        conda_env: Conda环境路径|Conda environment path

    Returns:
        bool: 成功返回True，失败返回False
               Returns True if successful, False otherwise
    """
    logger.info(f"创建基因组索引|Creating genome index: {genome_file}")

    samtools = f"{conda_env}/bin/samtools"
    cmd = [samtools, "faidx", genome_file]

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        logger.info("基因组索引创建成功|Genome index created successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"基因组索引创建失败|Failed to create genome index: {e.stderr}")
        return False
