"""
K-mer计数工具函数模块|K-mer Count Utility Functions Module
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Optional


class KmerCountLogger:
    """K-mer计数日志管理器|K-mer Count Logger Manager"""

    def __init__(self, output_dir: Path, log_file: str = "kmer_count.log", verbose: bool = False):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            log_file: 日志文件名|Log file name
            verbose: 详细模式|Verbose mode
        """
        self.output_dir = output_dir
        self.log_file = output_dir / log_file
        self.verbose = verbose
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = logging.DEBUG if self.verbose else logging.INFO

        # 清除现有处理器|Clear existing handlers
        logger = logging.getLogger('kmer_count')
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        logger.setLevel(level)

        # 创建格式化器|Create formatter
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 控制台处理器|Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)

        # 文件处理器|File handler
        file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir=None):
        self.logger = logger
        self.working_dir = working_dir or os.getcwd()
    
    def run(self, cmd, description: str = "", check: bool = True) -> bool:
        """执行命令|Execute command

        注意|Note:
            cmd 可为列表或字符串; 字符串按空白拆分后以 shell=False 执行
            (等价于 shell 解析,但不支持管道/通配符/变量,更安全)
            cmd can be a list or string; a string is split on whitespace and
            executed with shell=False (safer; no pipes/globs/vars)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 统一为列表,使用shell=False更安全|Normalize to list, use shell=False (safer)
        cmd_list = cmd.split() if isinstance(cmd, str) else cmd
        self.logger.info(f"命令|Command: {' '.join(cmd_list)}")

        try:
            result = subprocess.run(
                cmd_list,
                capture_output=True,
                text=True,
                check=check,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功|Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"标准错误|Stderr: {result.stderr}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            if check:
                raise
            return False


def setup_logger(output_dir: Path, verbose: bool = False) -> logging.Logger:
    """设置日志记录器|Setup logger"""
    logger = logging.getLogger('kmer_count')
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # 清除现有处理器|Clear existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 控制台处理器|Console handler
    # console_handler = logging.StreamHandler()
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # 文件处理器(§12: 日志集中到 99_logs/)|File handler (§12: logs centralized to 99_logs/)
    log_dir = output_dir / '99_logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / 'kmer_count.log'
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    
    # 格式设置|Format setting
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    return logger


def check_dependencies(jellyfish_path: str, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")
    
    try:
        result = subprocess.run([jellyfish_path, "--version"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            logger.info(f"Jellyfish 可用|available: {version}")
            return True
        else:
            logger.error(f"Jellyfish 不可用|not available")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"Jellyfish 不可用|not available: {e}")
        return False

class FileIntegrityError(Exception):
    """文件完整性错误|File Integrity Error"""
    pass