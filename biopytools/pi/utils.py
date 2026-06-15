"""
Pi计算工具函数模块|Pi Calculation Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class PiRow:
    """Pi计算结果行|Pi calculation result row"""
    population: str
    chromosome: str
    window_start: Optional[int]  # None for genome-wide mode
    window_end: Optional[int]    # None for genome-wide mode
    pi_value: float
    n_sites: int = 0
    source: str = ""  # "vcftools" or "pixy"


class PiLogger:
    """Pi计算日志管理器|Pi Calculation Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "pi_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / '99_logs' / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 清空旧日志|Clear old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler - INFO级别|Stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上|Stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler, stderr_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, return environment name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 如果传入的是完整路径，直接从路径提取
    # Method 2: If full path is passed, extract directly from path
    if os.path.exists(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件
    Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用，添加--no-capture-output避免内存问题
        # Use conda run with --no-capture-output to avoid memory issues
        cmd_name = Path(command).name
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', cmd_name] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "", timeout: int = 3600) -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
            timeout: 超时时间(秒)|Timeout in seconds
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 记录完整命令（INFO级别）|Log complete command (INFO level)
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                timeout=timeout,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:1000]}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command execution timeout: {description}")
            return False


def parse_population_file(pop_file: str) -> Dict[str, List[str]]:
    """
    解析群体文件，返回群体到样本列表的映射
    Parse population file, return mapping of population to sample list

    Args:
        pop_file: 群体文件路径|Population file path

    Returns:
        {群体名: [样本列表]}|{population_name: [sample_list]}
    """
    pop_to_samples: Dict[str, List[str]] = {}

    with open(pop_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # 自动检测分隔符|Auto-detect separator
            if '\t' in line:
                parts = line.split('\t')
            elif ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()

            if len(parts) >= 2:
                sample_id = parts[0]
                population = parts[1]
                if population not in pop_to_samples:
                    pop_to_samples[population] = []
                pop_to_samples[population].append(sample_id)

    return pop_to_samples


def get_software_version(tool_path: str, logger) -> str:
    """
    自动检测软件版本|Auto-detect software version

    Args:
        tool_path: 工具路径|Tool path
        logger: 日志器|Logger

    Returns:
        版本字符串|Version string
    """
    try:
        cmd = build_conda_command(tool_path, ['--version'])
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            return version
        else:
            return "unknown"
    except Exception as e:
        logger.warning(f"版本检测失败|Version detection failed: {e}")
        return "unknown"


def parse_fai_file(fai_file: str) -> Dict[str, int]:
    """
    解析samtools faidx生成的.fai文件，获取各染色体长度
    Parse .fai file generated by samtools faidx, get chromosome lengths

    Args:
        fai_file: .fai文件路径|.fai file path

    Returns:
        {染色体名: 碱基长度}|{chromosome_name: base_length}
    """
    chrom_lengths: Dict[str, int] = {}

    with open(fai_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])

    return chrom_lengths


def format_number(num: float) -> str:
    """格式化数字|Format number"""
    return f"{num:.6f}"


def is_step_completed(output_file: str) -> bool:
    """检查步骤是否已完成（通过输出文件存在性判断）|Check if step is done"""
    return Path(output_file).exists()


def ensure_tabix_index(vcf_file: str, vcftools_path: str, logger) -> bool:
    """
    检查VCF文件的tabix索引(.tbi)是否存在，不存在则自动创建
    Check if VCF tabix index (.tbi) exists, create if missing

    Args:
        vcf_file: VCF文件路径|VCF file path
        vcftools_path: vcftools路径（用于定位conda环境中的tabix）|vcftools path (to locate tabix in conda env)
        logger: 日志器|Logger

    Returns:
        True if index exists or was created successfully, False otherwise
    """
    tbi_file = vcf_file + '.tbi'

    if os.path.exists(tbi_file):
        logger.info(f"VCF索引已存在|VCF index already exists: {tbi_file}")
        return True

    logger.info(f"VCF索引不存在，正在创建|VCF index not found, creating: {tbi_file}")

    # 使用vcftools所在conda环境中的tabix（同一环境通常都有htslib/tabix）
    # Use tabix from the same conda env as vcftools
    tabix_path = vcftools_path.replace('vcftools', 'tabix')
    if not os.path.exists(tabix_path):
        # 回退到直接用vcftools路径所在目录找tabix
        # Fallback: look for tabix in the same directory
        tabix_path = 'tabix'

    cmd = build_conda_command(tabix_path, ['-p', 'vcf', vcf_file])
    try:
        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            timeout=7200  # 大文件建索引可能需要较长时间
        )

        if os.path.exists(tbi_file):
            logger.info(f"VCF索引创建成功|VCF index created successfully: {tbi_file}")
            return True
        else:
            logger.error(f"tabix命令执行成功但索引文件未生成|tabix command succeeded but index file not created")
            return False

    except subprocess.CalledProcessError as e:
        logger.error(f"创建VCF索引失败|Failed to create VCF index")
        logger.error(f"错误代码|Error code: {e.returncode}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr[:1000]}")
        return False
    except subprocess.TimeoutExpired:
        logger.error(f"创建VCF索引超时|VCF index creation timeout")
        return False
