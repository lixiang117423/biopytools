# """
# FASTP质控工具函数模块 | FASTP Quality Control Utility Functions Module
# """

# import logging
# import subprocess
# import sys
# from pathlib import Path

# class FastpLogger:
#     """FASTP质控日志管理器 | FASTP Quality Control Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "fastp_processing.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 | Setup logging"""
#         # 确保输出目录存在 | Ensure output directory exists
#         self.output_dir.mkdir(parents=True, exist_ok=True)
        
#         # 如果日志文件存在，删除旧的 | Remove old log file if exists
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
    
#     def run(self, cmd: list, description: str = "") -> bool:
#         """执行命令 | Execute command"""
#         if description:
#             self.logger.info(f"执行步骤 | Executing step: {description}")
        
#         # 将命令列表转换为字符串用于日志 | Convert command list to string for logging
#         cmd_str = " ".join(str(arg) for arg in cmd)
#         self.logger.info(f"命令 | Command: {cmd_str}")
        
#         try:
#             result = subprocess.run(
#                 cmd, 
#                 check=True, 
#                 capture_output=True, 
#                 text=True
#             )
            
#             self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
#             # 如果有输出，记录到debug级别 | Log output to debug level if exists
#             if result.stdout.strip():
#                 self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
#             return True
            
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"命令执行失败 | Command execution failed: {description}")
#             self.logger.error(f"错误代码 | Error code: {e.returncode}")
#             self.logger.error(f"错误信息 | Error message: {e.stderr}")
#             return False
    
#     def check_executable(self, executable_path: str) -> bool:
#         """检查可执行文件是否可用 | Check if executable is available"""
#         try:
#             result = subprocess.run(
#                 [executable_path, "--version"], 
#                 capture_output=True, 
#                 text=True, 
#                 timeout=10
#             )
#             return result.returncode == 0
#         except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
#             return False

"""
🔧 FASTP质控工具函数模块 | FASTP Quality Control Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class FastpLogger:
    """📝 FASTP质控日志管理器 | FASTP Quality Control Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "fastp_processing.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """⚙️ 设置日志 | Setup logging"""
        # 📂 确保输出目录存在 | Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 🗑️ 如果日志文件存在，删除旧的 | Remove old log file if exists
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
        """📝 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """💻 命令执行器 | Command Runner"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def run(self, cmd: list, description: str = "") -> bool:
        """⚡ 执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        # 📝 将命令列表转换为字符串用于日志 | Convert command list to string for logging
        cmd_str = " ".join(str(arg) for arg in cmd)
        self.logger.info(f"💻 命令 | Command: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            # 📊 如果有输出，记录到debug级别 | Log output to debug level if exists
            if result.stdout.strip():
                self.logger.debug(f"📊 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            return False
    
    def check_executable(self, executable_path: str) -> bool:
        """🔍 检查可执行文件是否可用 | Check if executable is available"""
        try:
            result = subprocess.run(
                [executable_path, "--version"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            return False