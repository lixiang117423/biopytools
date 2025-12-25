"""
CNCB FTP工具函数模块 | CNCB FTP Utility Functions Module
"""

import time
import logging
from pathlib import Path
from typing import Optional

class Logger:
    """日志管理器 | Logger Manager"""
    
    def __init__(self, name: str = "cncb_ftp_finder"):
        self.logger = logging.getLogger(name)
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '[%(asctime)s] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)
    
    def info(self, message: str):
        """信息日志 | Info log"""
        self.logger.info(f"ℹ️  {message}")
    
    def warning(self, message: str):
        """警告日志 | Warning log"""
        self.logger.warning(f"⚠️  {message}")
    
    def error(self, message: str):
        """错误日志 | Error log"""
        self.logger.error(f"❌ {message}")
    
    def success(self, message: str):
        """成功日志 | Success log"""
        self.logger.info(f"✅ {message}")
    
    def debug(self, message: str):
        """调试日志 | Debug log"""
        self.logger.debug(f"🔍 {message}")

def log_message(message: str, logger: Optional[Logger] = None):
    """兼容旧版本的日志函数 | Compatible legacy log function"""
    if logger:
        logger.info(message)
    else:
        timestamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
        print(f"{timestamp} {message}")

class PathValidator:
    """路径验证器 | Path Validator"""
    
    @staticmethod
    def validate_input_file(file_path: str) -> bool:
        """验证输入文件 | Validate input file"""
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"❌ 输入文件不存在 | Input file not found: {file_path}")
        if not path.is_file():
            raise ValueError(f"❌ 路径不是文件 | Path is not a file: {file_path}")
        return True
    
    @staticmethod
    def ensure_output_dir(file_path: str) -> Path:
        """确保输出目录存在 | Ensure output directory exists"""
        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path
