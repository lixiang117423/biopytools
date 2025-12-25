"""
📝 K-mer分析日志管理模块 | K-mer Analysis Logger Management Module
"""

import logging
import sys
from pathlib import Path

class KmerLogger:
    """K-mer分析日志管理器 | K-mer Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "metagraph_kmer.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 创建格式化器 | Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # 文件处理器 | File handler
        file_handler = logging.FileHandler(self.log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # 控制台处理器 | Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        
        # 配置根日志器 | Configure root logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger
