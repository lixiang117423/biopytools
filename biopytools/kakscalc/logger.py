"""
📊 Ka/Ks Calculator日志管理模块
功能: 提供emoji增强的日志记录功能 | Emoji-enhanced logging functionality
"""

import logging
import sys
from typing import Optional
from datetime import datetime

class Logger:
    """📊 支持emoji的日志工具类 | Emoji-enabled logging utility"""
    
    def __init__(self, log_file: Optional[str] = None, verbose: bool = False, name: str = 'KaKsAnalysis'):
        """
        🏗️ 初始化日志器 | Initialize logger
        
        Args:
            log_file: 日志文件路径 | Log file path
            verbose: 详细模式 | Verbose mode
            name: 日志器名称 | Logger name
        """
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG if verbose else logging.INFO)
        
        # 🗑️ 清除现有处理器 | Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
        
        # 🖥️ 控制台处理器 | Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%H:%M:%S'
        )
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
        
        # 📄 文件处理器 | File handler
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_formatter = logging.Formatter(
                '%(asctime)s | %(levelname)s | %(funcName)s:%(lineno)d | %(message)s'
            )
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)
    
    def info(self, message: str, emoji: str = "ℹ️"):
        """📢 信息日志 | Info logging"""
        self.logger.info(f"{emoji} {message}")
    
    def success(self, message: str):
        """✅ 成功日志 | Success logging"""
        self.logger.info(f"✅ {message}")
    
    def warning(self, message: str):
        """⚠️ 警告日志 | Warning logging"""
        self.logger.warning(f"⚠️ {message}")
    
    def error(self, message: str):
        """❌ 错误日志 | Error logging"""
        self.logger.error(f"❌ {message}")
    
    def debug(self, message: str):
        """🔍 调试日志 | Debug logging"""
        self.logger.debug(f"🔍 {message}")
    
    def progress(self, message: str, current: int, total: int):
        """📊 进度日志 | Progress logging"""
        percentage = (current / total) * 100 if total > 0 else 0
        self.logger.info(f"📊 {message} [{current}/{total}] ({percentage:.1f}%)")
    
    def separator(self, title: str = ""):
        """🔹 分隔符日志 | Separator logging"""
        if title:
            self.logger.info(f"🔹 {title} " + "="*50)
        else:
            self.logger.info("="*60)
