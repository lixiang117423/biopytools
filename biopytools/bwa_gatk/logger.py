"""
📝 日志管理模块 | Logger Management Module
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

class PipelineLogger:
    """流程日志管理器 | Pipeline Logger Manager"""
    
    def __init__(self, logs_dir: Path, verbose: bool = False):
        self.logs_dir = logs_dir
        self.verbose = verbose
        self.log_file = logs_dir / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 确定日志级别 | Determine log level
        log_level = logging.DEBUG if self.verbose else logging.INFO
        
        # 配置日志格式 | Configure log format
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        
        # 清空现有handlers | Clear existing handlers
        logging.root.handlers = []
        
        # 配置logging | Configure logging
        logging.basicConfig(
            level=log_level,
            format=log_format,
            datefmt=date_format,
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        
        # 记录日志配置 | Log configuration
        self.logger.info("=" * 80)
        self.logger.info("🧬 BWA-GATK变异检测流程启动 | BWA-GATK Variant Calling Pipeline Started")
        self.logger.info("=" * 80)
        self.logger.info(f"📝 日志文件 | Log file: {self.log_file}")
        self.logger.info(f"📊 详细模式 | Verbose mode: {'开启' if self.verbose else '关闭'}")
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger
