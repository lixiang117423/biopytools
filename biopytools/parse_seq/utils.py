"""
序列提取工具函数模块|Sequence Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Tuple

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command

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
