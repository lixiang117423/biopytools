"""
SRA转换工具函数模块|SRA Conversion Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path


def format_number(num: int) -> str:
    """格式化数字|Format number

    Args:
        num: 要格式化的数字|Number to format

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.1f}K"
    return str(num)

class ConvertLogger:
    """转换日志管理器|Conversion Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "sra_conversion.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")
        
        self.logger.info(f"命令|Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功|Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")
    
    # 首先检查parallel-fastq-dump
    parallel_available = False
    fastq_dump_available = False
    
    try:
        result = subprocess.run(['parallel-fastq-dump', '-V'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info(f"parallel-fastq-dump 可用 (推荐)|available (recommended)")
            logger.info(f"   版本|Version: {result.stdout.strip()}")
            parallel_available = True
            config.use_parallel = True
            config.tool_path = 'parallel-fastq-dump'
    except (subprocess.TimeoutExpired, FileNotFoundError):
        logger.warning(f"parallel-fastq-dump 未找到|not found")
    
    # 检查fastq-dump作为备选
    try:
        result = subprocess.run(['fastq-dump', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info(f"fastq-dump 可用 (备选)|available (fallback)")
            fastq_dump_available = True
            if not parallel_available:
                config.use_parallel = False
                config.tool_path = 'fastq-dump'
                logger.warning(f"将使用fastq-dump (速度较慢)|Will use fastq-dump (slower)")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        logger.warning(f"fastq-dump 未找到|not found")
    
    if not parallel_available and not fastq_dump_available:
        error_msg = "缺少必需软件: 需要 parallel-fastq-dump 或 fastq-dump|Missing required software: need parallel-fastq-dump or fastq-dump"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True
