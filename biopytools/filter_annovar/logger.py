"""
日志管理模块 | Logger Management Module
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

class FilterLogger:
    """变异提取日志管理器 | Variant Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "filter_annovar.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger
