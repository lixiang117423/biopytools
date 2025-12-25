"""
基因组组装日志管理模块 📝 | Genome Assembly Logger Module
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

class AssemblyLogger:
    """组装日志管理器 | Assembly Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "assembly.log"):
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
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger
    
    def log_command(self, cmd: str, description: str = ""):
        """记录命令 | Log command"""
        if description:
            self.logger.info(f"🔧 {description}")
        self.logger.info(f"💻 执行命令 | Executing command: {cmd}")
    
    def log_step(self, step_name: str, details: str = ""):
        """记录步骤 | Log step"""
        self.logger.info(f"🚀 {step_name}")
        if details:
            self.logger.info(f"   {details}")
    
    def log_error(self, error_msg: str):
        """记录错误 | Log error"""
        self.logger.error(f"❌ {error_msg}")
    
    def log_success(self, success_msg: str):
        """记录成功 | Log success"""
        self.logger.info(f"✅ {success_msg}")
    
    def log_warning(self, warning_msg: str):
        """记录警告 | Log warning"""
        self.logger.warning(f"⚠️ {warning_msg}")
    
    def log_info(self, info_msg: str):
        """记录信息 | Log info"""
        self.logger.info(f"ℹ️ {info_msg}")
