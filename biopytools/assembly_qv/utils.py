"""
Merqury QV计算工具函数模块|Merqury QV Calculation Utility Functions Module
"""

import logging
import sys
import subprocess
import os
import re
import shutil
from pathlib import Path
from typing import List, Tuple, Optional


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


class MerquryQVLogger:
    """Merqury QV计算日志管理器|Merqury QV Calculation Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 配置日志格式|Configure logging format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 设置日志级别|Set log level
        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        # 配置logging|Configure logging
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


def find_r1_r2_pairs(fastq_files: List[str]) -> Tuple[List[Tuple[str, str]], List[str]]:
    """
    自动识别并配对R1和R2文件|Automatically identify and pair R1 and R2 files

    Args:
        fastq_files: FASTQ文件列表|List of FASTQ files

    Returns:
        tuple: (配对列表, 单端文件列表)|(pairs_list, single_end_files_list)
            - pairs: [(R1_file, R2_file), ...]
            - single_end: [file1, file2, ...]
    """
    # 分离单端和双端数据|Separate single-end and paired-end data
    r1_files = {}
    r2_files = {}
    single_end_files = []

    for file_path in fastq_files:
        filename = Path(file_path).name

        # 检测R1/R2标记|Detect R1/R2 markers
        if re.search(r'[_\.]R1[_\.]', filename) or re.search(r'[_\.]read1[_\.]', filename):
            # 提取样本名|Extract sample name
            sample = re.sub(r'[_\.](R1|read1)[_.].*$', '', filename)
            r1_files[sample] = file_path
        elif re.search(r'[_\.]R2[_\.]', filename) or re.search(r'[_\.]read2[_\.]', filename):
            sample = re.sub(r'[_\.](R2|read2)[_.].*$', '', filename)
            r2_files[sample] = file_path
        else:
            # 单端数据或无法识别|Single-end or unrecognized
            single_end_files.append(file_path)

    # 配对R1和R2|Pair R1 and R2
    pairs = []
    matched_samples = set()

    # 匹配双端数据|Match paired-end data
    for sample in r1_files.keys():
        if sample in r2_files:
            pairs.append((r1_files[sample], r2_files[sample]))
            matched_samples.add(sample)

    # 未匹配的R1或R2作为单端处理|Unmatched R1 or R2 as single-end
    for sample in r1_files.keys():
        if sample not in matched_samples:
            single_end_files.append(r1_files[sample])

    for sample in r2_files.keys():
        if sample not in matched_samples:
            single_end_files.append(r2_files[sample])

    return pairs, single_end_files


def detect_data_type(fastq_files: List[str]) -> str:
    """
    自动检测数据类型|Automatically detect data type

    Args:
        fastq_files: FASTQ文件列表|List of FASTQ files

    Returns:
        str: 数据类型|Data type (illumina, hifi, unknown)
    """
    if not fastq_files:
        return "unknown"

    # 检查文件名特征|Check filename patterns
    has_r1_r2 = any(re.search(r'[_\.]R[12][_\.]', f) for f in fastq_files)
    has_hifi = any(re.search(r'hifi|ccs|subreads', f, re.IGNORECASE) for f in fastq_files)

    if has_hifi:
        return "hifi"
    elif has_r1_r2:
        return "illumina"
    else:
        # 默认认为是Illumina|Default to illumina
        return "illumina"


def run_command(command: List[str], conda_env: str, logger, description: str = "") -> Tuple[bool, str]:
    """
    在conda环境中运行命令|Run command in conda environment

    Args:
        command: 命令列表|Command list
        conda_env: Conda环境路径|Conda environment path (保留用于兼容性|Kept for compatibility)
        logger: 日志器|Logger
        description: 命令描述|Command description

    Returns:
        tuple: (成功状态, 输出)|(success, output)
    """
    if description:
        logger.info(f"执行|Executing: {description}")

    # 设置环境变量|Set environment variables
    env = os.environ.copy()
    # 保留MERQURY环境变量设置|Keep MERQURY environment variable setting
    env['MERQURY'] = conda_env.replace("/bin/", "/share/merqury/")

    # 构建完整命令|Build full command
    full_command = " ".join(command)

    # 自动包装conda环境的命令|Automatically wrap conda environment commands
    # 提取命令名称|Extract command name
    first_part = command[0] if command else ""
    if " " in first_part:
        # 如果第一个参数包含空格（如"meryl count"），取第一个词|If first part contains space (like "meryl count"), take first word
        cmd_name = first_part.split()[0]
    else:
        # 否则取整个第一个参数|Otherwise take entire first argument
        cmd_name = os.path.basename(first_part) if first_part else ""

    # 使用conda wrapper包装命令|Use conda wrapper to wrap command
    if cmd_name:
        wrapped_command = build_conda_command_string(cmd_name, full_command[len(first_part):])
    else:
        wrapped_command = full_command

    logger.debug(f"命令|Command: {wrapped_command}")

    try:
        result = subprocess.run(
            wrapped_command,
            shell=True,
            capture_output=True,
            text=True,
            env=env,
            check=True
        )

        logger.debug(f"命令执行成功|Command executed successfully: {description}")
        return True, result.stdout

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {description}")
        logger.error(f"错误代码|Error code: {e.returncode}")
        logger.error(f"错误信息|Error message: {e.stderr}")
        return False, e.stderr


def parse_qv_output(qv_file: str) -> dict:
    """
    解析Merqury QV输出文件|Parse Merqury QV output file

    Merqury .qv文件格式|Merqury .qv file format:
    列1|Col1: 序列名称|Sequence name
    列2|Col2: N50或contig数量|N50 or contig count
    列3|Col3: 总碱基数量|Total base count
    列4|Col4: QV值|QV value
    列5|Col5: 错误率|Error rate

    Args:
        qv_file: QV文件路径|QV file path

    Returns:
        dict: QV统计信息|QV statistics
    """
    stats = {}

    if not os.path.exists(qv_file):
        return stats

    with open(qv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # 解析tab分隔的数据|Parse tab-separated data
            parts = line.split('\t')
            if len(parts) >= 5:
                stats['sequence_name'] = parts[0]
                stats['n50_or_contigs'] = parts[1]
                stats['total_bases'] = parts[2]
                stats['qv_value'] = parts[3]
                stats['error_rate'] = parts[4]

                # 添加原始行用于参考|Add raw line for reference
                stats['raw_line'] = line
                break  # 只取第一行数据|Only take first line

    return stats if stats else {'raw_output': f"无法解析QV文件|Failed to parse QV file: {qv_file}"}
