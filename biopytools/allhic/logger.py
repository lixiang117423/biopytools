"""
ALLHiC流水线日志管理模块 📝 | ALLHiC Pipeline Logger Module
"""

import logging
import os
import sys
import time
from pathlib import Path
from datetime import datetime

class PipelineLogger:
    """流水线日志管理器 | Pipeline Logger Manager"""
    
    def __init__(self, work_dir: str, log_name: str = "allhic_pipeline"):
        self.work_dir = Path(work_dir)
        self.log_file = self.work_dir / "logs" / f"{log_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        self.step_start_time = 0
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        
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
    
    def log(self, msg: str):
        """记录消息 | Log message"""
        self.logger.info(msg)

    def log_info(self, msg: str):
        """记录信息 | Log info"""
        self.logger.info(msg)
    
    def log_section(self, title: str):
        """记录章节标题 | Log section title"""
        separator = "━" * 70
        self.logger.info(f"\n{separator}\n{title}\n{separator}")
        self.step_start_time = time.time()
    
    def log_step_time(self):
        """记录步骤耗时 | Log step time"""
        if self.step_start_time > 0:
            elapsed = int(time.time() - self.step_start_time)
            formatted_time = self.format_time(elapsed)
            self.logger.info(f"⏱️  步骤完成 | Step completed in {formatted_time}")
    
    def format_time(self, seconds: int) -> str:
        """格式化时间 | Format time"""
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)

        if hours > 0:
            return f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            return f"{minutes}m {seconds}s"
        else:
            return f"{seconds}s"

    def format_file_size(self, file_path: str) -> str:
        """格式化文件大小 | Format file size"""
        try:
            size = os.path.getsize(file_path)
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size < 1024.0:
                    return f"{size:.1f}{unit}"
                size /= 1024.0
            return f"{size:.1f}TB"
        except OSError:
            return "未知大小"
    
    def log_command(self, cmd: str, description: str = ""):
        """记录命令 | Log command"""
        if description:
            self.logger.info(f"🔧 {description}")
        self.logger.info(f"💻 执行命令 | Executing command: {cmd}")
    
    def log_error(self, error_msg: str):
        """记录错误 | Log error"""
        self.logger.error(f"❌ {error_msg}")
    
    def log_success(self, success_msg: str):
        """记录成功 | Log success"""
        self.logger.info(f"✅ {success_msg}")
    
    def log_warning(self, warning_msg: str):
        """记录警告 | Log warning"""
        self.logger.warning(f"⚠️ {warning_msg}")
    
    def run_command(self, cmd: str, description: str = "", work_dir: str = None) -> bool:
        """执行命令并记录 | Execute command and log"""
        self.log_command(cmd, description)
        
        try:
            import subprocess
            result = subprocess.run(
                cmd, 
                shell=True, 
                cwd=work_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                self.logger.debug(f"输出 | Output: {result.stdout}")
            
            self.log_success("命令执行成功 | Command executed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            self.log_error(f"命令执行失败 | Command execution failed: {e}")
            if e.stderr:
                self.log_error(f"错误信息 | Error message: {e.stderr}")
            return False
