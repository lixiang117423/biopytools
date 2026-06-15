"""
序列提取工具函数模块|Sequence Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
import re
import os
from pathlib import Path
from typing import List, Tuple, Optional

class ExtractorLogger:
    """序列提取日志管理器|Sequence Extraction Logger Manager"""
    
    def __init__(self, output_file: str, log_name: str = "seq_extraction.log", verbose: bool = True):
        self.output_dir = Path(output_file).parent
        self.log_file = self.output_dir / log_name
        self.verbose = verbose
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 日志格式|Log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 创建logger|Create logger
        self.logger = logging.getLogger(f"seq_extractor_{id(self)}")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

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

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir.resolve() if working_dir else Path.cwd()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 |Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")
        
        self.logger.info(f"命令|Command: {cmd}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功 |Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 |Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    策略|Strategy:
    1. 首先尝试从which命令路径检测（优先级高）|First try detecting from which command path (high priority)
    2. 如果未找到，搜索所有conda环境（兜底方案）|If not found, search all conda environments (fallback)

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'samtools' or '/path/to/samtools')

    Returns:
        conda环境名称或None|conda environment name or None (e.g., 'GATK_v.4.6.2.0' or None)
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        # 例如: /miniforge3/envs/GATK_v.4.6.2.0/bin/samtools
        # e.g.: /miniforge3/envs/GATK_v.4.6.2.0/bin/samtools
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda|CONDA_EXE is usually /path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
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
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表 (适用于subprocess.run(shell=False))|Complete command list (for subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('samtools', ['view', 'file.bam'])
        ['conda', 'run', '-n', 'env_name', '--no-capture-output', 'samtools', 'view', 'file.bam']

    注意|Note:
        返回的列表应配合 subprocess.run(shell=False) 使用|The returned list must be used with subprocess.run(shell=False)

    ⚠️ 重要|IMPORTANT:
        必须使用--no-capture-output避免conda缓冲输出导致内存问题|Must use --no-capture-output to avoid conda buffering output causing memory issues
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用|Use conda run to invoke
        # 如果command是命令名，conda run会自动找到环境中的版本|If command is a name, conda run will auto-find the version in the environment
        # 添加--no-capture-output避免缓冲输出导致内存问题|Add --no-capture-output to avoid buffering output causing memory issues
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd


def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = []

    # DNA序列需要samtools|DNA sequences need samtools
    if config.sequence_type == "dna":
        dependencies.append((config.samtools_path, "samtools"))

    missing_deps = []

    for cmd, name in dependencies:
        try:
            # 使用conda wrapper构建命令|Use conda wrapper to build command
            wrapped_cmd = build_conda_command(cmd, ['--version'])
            result = subprocess.run(wrapped_cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f" {name} 可用|available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)

    if missing_deps:
        error_msg = f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True

def parse_regions_file(regions_file: str) -> List[Tuple[str, int, int, str]]:
    """解析区域文件 |Parse regions file"""
    regions = []
    
    with open(regions_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                parts = line.split()  # 尝试空格分隔|Try space separation
            
            if len(parts) < 3:
                raise ValueError(f"区域文件第{line_num}行格式错误，需要至少3列|"
                               f"Invalid format in regions file line {line_num}, need at least 3 columns: {line}")
            
            try:
                chrom = parts[0]
                start = int(parts[1])  # 1-based
                end = int(parts[2])    # 1-based
                
                # 检查第四列是否有链信息|Check if 4th column has strand info
                strand = '+'  # 默认正链|Default positive strand
                if len(parts) >= 4:
                    strand_col = parts[3].strip()
                    if strand_col in ['+', '-']:
                        strand = strand_col
                
                if start <= 0 or end <= 0:
                    raise ValueError(f"坐标必须为正整数|Coordinates must be positive integers")
                
                if start > end:
                    raise ValueError(f"起始位置不能大于终止位置|Start position cannot be greater than end position")
                
                regions.append((chrom, start, end, strand))
                
            except ValueError as e:
                raise ValueError(f"区域文件第{line_num}行坐标解析错误|"
                               f"Coordinate parsing error in regions file line {line_num}: {e}")
    
    return regions

def format_region_name(chrom: str, start: int, end: int, original_header: str = None) -> str:
    """格式化区域名称 |Format region name"""
    region_info = f"{chrom}:{start}-{end}"
    
    if original_header:
        # 移除原始header中的'>'符号|Remove '>' from original header
        clean_header = original_header.lstrip('>')
        
        # 检查是否已经包含区域信息，避免重复|Check if region info already exists to avoid duplication
        if region_info in clean_header:
            # 如果已经包含区域信息，直接返回|If region info already exists, return as is
            return f">{clean_header}"
        else:
            # 如果原始header是染色体名，则用区域信息替换|If original header is chromosome name, replace with region info
            if clean_header == chrom:
                return f">{region_info}"
            else:
                return f">{clean_header}_{region_info}"
    else:
        return f">region_{region_info}"
