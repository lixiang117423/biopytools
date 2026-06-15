"""
微观共线性分析工具函数|Microsynteny Analysis Utilities
"""

import logging
import os
import sys
import subprocess
import shutil
import re
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
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
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, command] + args
    else:
        return [command] + args


class MicrosyntenyLogger:
    """微观共线性分析日志管理器|Microsynteny Logger Manager"""

    def __init__(self, log_file: Optional[Path] = None, log_level: str = "INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR/CRITICAL)
        """
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        # 格式: YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 消息中文|Message English
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

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
        """获取日志器|Get logger

        Returns:
            logging.Logger: 日志器对象|Logger object
        """
        return self.logger


def run_command(
    command: List[str],
    logger: logging.Logger,
    cwd: Optional[str] = None,
    env: Optional[dict] = None
) -> bool:
    """执行命令（自动检测conda环境）|Execute command (auto-detect conda environment)

    Args:
        command: 命令列表|Command list
        logger: 日志器|Logger
        cwd: 工作目录|Working directory
        env: 环境变量|Environment variables

    Returns:
        bool: 是否成功|Whether successful
    """
    # 自动包装conda环境的命令|Auto-wrap conda environment commands
    if command:
        cmd_name = os.path.basename(command[0])
        wrapped_cmd = build_conda_command(cmd_name, command[1:])
    else:
        wrapped_cmd = command

    logger.info(f"执行命令|Executing command: {' '.join(wrapped_cmd)}")

    # 确保环境变量包含必要路径|Ensure env includes necessary paths
    if env is None:
        env = os.environ.copy()

    try:
        result = subprocess.run(
            wrapped_cmd,
            cwd=cwd,
            env=env,
            check=True,
            capture_output=True,
            text=True
        )

        if result.stdout:
            logger.debug(f"标准输出|Stdout: {result.stdout}")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"执行异常|Execution error: {str(e)}")
        return False


def run_jcvi_command(
    jcvi_python: str,
    module: str,
    action: str,
    args: List[str],
    logger: logging.Logger
) -> bool:
    """执行JCVI命令|Execute JCVI command

    Args:
        jcvi_python: JCVI Python路径|JCVI Python path
        module: JCVI模块名|JCVI module name
        action: JCVI动作名|JCVI action name
        args: 命令参数|Command arguments
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Whether successful
    """
    command = [jcvi_python, "-m", f"jcvi.{module}", action] + args

    # 设置环境变量，确保JCVI能找到LAST等依赖工具|Set env to ensure JCVI finds LAST etc.
    env = os.environ.copy()
    jcvi_bin = os.path.dirname(os.path.dirname(jcvi_python))  # .../envs/jcvi_v.1.5.7
    jcvi_bin = os.path.join(jcvi_bin, "bin")

    # 将JCVI的bin目录添加到PATH最前面|Add JCVI bin to the front of PATH
    if "PATH" in env:
        env["PATH"] = f"{jcvi_bin}:{env['PATH']}"
    else:
        env["PATH"] = jcvi_bin

    logger.debug(f"设置PATH|Set PATH: {env['PATH']}")

    return run_command(command, logger, env=env)


def check_file_exists(file_path: str, logger: logging.Logger) -> bool:
    """检查文件是否存在|Check if file exists

    Args:
        file_path: 文件路径|File path
        logger: 日志器|Logger

    Returns:
        bool: 是否存在|Whether exists
    """
    from pathlib import Path

    path = Path(file_path)
    if not path.exists():
        logger.warning(f"文件不存在|File not found: {file_path}")
        return False

    return True


def format_number(num: int) -> str:
    """格式化数字|Format number

    大数字使用M(百万)单位显示|Large numbers use M (million) unit

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.1f}K"
    return str(num)


def parse_gff_gene_ids(gff_file: str, logger: logging.Logger) -> List[str]:
    """从GFF文件中提取基因ID|Extract gene IDs from GFF file

    Args:
        gff_file: GFF文件路径|GFF file path
        logger: 日志器|Logger

    Returns:
        List[str]: 基因ID列表|List of gene IDs
    """
    gene_ids = []

    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                # 提取ID属性|Extract ID attribute
                attributes = parts[8]
                if 'ID=' in attributes:
                    gene_id = attributes.split('ID=')[1].split(';')[0]
                    gene_ids.append(gene_id)

        logger.info(f"从GFF文件中提取|Extracted from GFF: {format_number(len(gene_ids))} 个基因|genes")

    except Exception as e:
        logger.error(f"读取GFF文件失败|Failed to read GFF file: {str(e)}")

    return gene_ids


def find_gene_in_bed(gene_id: str, bed_file: str, logger: logging.Logger) -> Optional[tuple]:
    """在BED文件中查找基因|Find gene in BED file

    Args:
        gene_id: 基因ID|Gene ID
        bed_file: BED文件路径|BED file path
        logger: 日志器|Logger

    Returns:
        Optional[tuple]: (chr, start, end, strand) 或 None
    """
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 4 and parts[3] == gene_id:
                    chr_name = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand = parts[5] if len(parts) > 5 else '+'
                    return (chr_name, start, end, strand)

    except Exception as e:
        logger.error(f"读取BED文件失败|Failed to read BED file: {str(e)}")

    return None
