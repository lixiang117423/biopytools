"""
PacBio SV专用日志模块 | PacBio SV Dedicated Logger Module
"""

import logging
import sys
from pathlib import Path


class PacBioSVLogger:
    """PacBio SV专用日志器 | PacBio SV Dedicated Logger"""
    
    def __init__(self, output_dir: Path, log_name: str = "pacbio_sv_analysis.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / log_name
        self.logger = None
        self._setup_logger()
    
    def _setup_logger(self):
        """设置日志器 | Setup logger"""
        # 确保输出目录存在
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 如果日志文件存在则删除
        if self.log_file.exists():
            self.log_file.unlink()
        
        # 创建logger
        self.logger = logging.getLogger('pacbio_sv')
        self.logger.setLevel(logging.DEBUG)
        
        # 清除现有的handlers
        if self.logger.handlers:
            self.logger.handlers.clear()
        
        # 创建formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # 文件handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # 控制台handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        
        # 添加handlers
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def get_logger(self):
        """获取logger实例 | Get logger instance"""
        return self.logger
