# """
# RNA-seq分析工具函数模块 | RNA-seq Analysis Utility Functions Module
# """

# import os
# import logging
# import subprocess
# import sys
# from pathlib import Path

# class RNASeqLogger:
#     """RNA-seq分析日志管理器 | RNA-seq Analysis Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "rnaseq_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 | Setup logging"""
#         if self.log_file.exists():
#             self.log_file.unlink()
        
#         logging.basicConfig(
#             level=logging.INFO,
#             format='%(asctime)s - %(levelname)s - %(message)s',
#             handlers=[
#                 logging.FileHandler(self.log_file),
#                 logging.StreamHandler(sys.stdout)
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def get_logger(self):
#         """获取日志器 | Get logger"""
#         return self.logger

# class CommandRunner:
#     """命令执行器 | Command Runner"""
    
#     def __init__(self, logger):
#         self.logger = logger
    
#     def run(self, cmd: str, description: str = "") -> bool:
#         """执行命令 | Execute command"""
#         if description:
#             self.logger.info(f"运行 | Running: {description}")
        
#         self.logger.info(f"命令 | Command: {cmd}")
        
#         try:
#             result = subprocess.run(
#                 cmd, 
#                 shell=True, 
#                 check=True, 
#                 capture_output=True, 
#                 text=True
#             )
            
#             self.logger.info(f"✓ {description} 完成 | completed")
#             return True
            
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"✗ {description} 失败 | failed")
#             self.logger.error(f"错误信息 | Error message: {e.stderr}")
#             sys.exit(1)

# class FileValidator:
#     """文件验证器 | File Validator"""
    
#     def __init__(self, logger):
#         self.logger = logger
    
#     def check_file_exists(self, file_path: str, description: str = "") -> bool:
#         """检查文件是否存在 | Check if file exists"""
#         if os.path.exists(file_path):
#             if description:
#                 self.logger.info(f"✓ {description}已存在，跳过 | already exists, skipping: {file_path}")
#             return True
#         return False

"""
RNA-seq分析工具函数模块 | RNA-seq Analysis Utility Functions Module
"""

import os
import logging
import subprocess
import sys
from pathlib import Path

class RNASeqLogger:
    """RNA-seq分析日志管理器 | RNA-seq Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "rnaseq_analysis.log"):
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

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 运行 | Running: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                check=True, 
                capture_output=True, 
                text=True
            )
            
            self.logger.info(f"✅ {description} 完成 | completed")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ {description} 失败 | failed")
            self.logger.error(f"💥 错误信息 | Error message: {e.stderr}")
            sys.exit(1)

class FileValidator:
    """文件验证器 | File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在 | Check if file exists"""
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"⏭️  {description}已存在，跳过 | already exists, skipping: {file_path}")
            return True
        return False