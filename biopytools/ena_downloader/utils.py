"""
ENA下载工具辅助函数模块 | ENA Downloader Utilities Module
"""

import os
import logging
import sys
from pathlib import Path
from typing import Dict, Any

class DownloadLogger:
    """下载日志管理器 | Download Logger Manager"""
    
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.log_file = output_path / "ena_download.log"
        self._setup_logging()
    
    def _setup_logging(self):
        """设置日志配置 | Setup logging configuration"""
        # 创建logger | Create logger
        self.logger = logging.getLogger('ena_downloader')
        self.logger.setLevel(logging.DEBUG)
        
        # 避免重复添加handler | Avoid adding duplicate handlers
        if self.logger.handlers:
            return
        
        # 创建formatter | Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # 文件handler | File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # 控制台handler | Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        
        # 添加handlers | Add handlers
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def get_logger(self):
        """获取logger实例 | Get logger instance"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, output_path: Path):
        self.logger = logger
        self.output_path = output_path
    
    def run_command(self, cmd: str, shell: bool = True) -> bool:
        """执行系统命令 | Execute system command"""
        import subprocess
        
        self.logger.info(f"执行命令 | Executing command: {cmd}")
        
        try:
            # 改变到输出目录 | Change to output directory
            original_cwd = os.getcwd()
            os.chdir(self.output_path)
            
            result = subprocess.run(
                cmd,
                shell=shell,
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时 | 1 hour timeout
            )
            
            # 恢复原目录 | Restore original directory
            os.chdir(original_cwd)
            
            if result.returncode == 0:
                self.logger.info("命令执行成功 | Command executed successfully")
                if result.stdout:
                    self.logger.debug(f"标准输出 | STDOUT: {result.stdout}")
                return True
            else:
                self.logger.error(f"命令执行失败 | Command failed with return code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误输出 | STDERR: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            os.chdir(original_cwd)
            self.logger.error("命令执行超时 | Command execution timed out")
            return False
        except Exception as e:
            os.chdir(original_cwd)
            self.logger.error(f"命令执行异常 | Command execution error: {str(e)}")
            return False

def check_dependencies(config, logger) -> bool:
    """检查依赖软件 | Check dependencies"""
    dependencies = []
    
    # 检查wget (对于FTP下载) | Check wget (for FTP downloads)
    if config.protocol == "ftp":
        dependencies.append("wget")
    
    # 检查ascp (对于Aspera下载) | Check ascp (for Aspera downloads)
    if config.protocol == "aspera":
        dependencies.append("ascp")
    
    missing_deps = []
    for dep in dependencies:
        if not check_command_exists(dep):
            missing_deps.append(dep)
    
    if missing_deps:
        logger.error(f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}")
        return False
    
    logger.info("所有依赖软件检查通过 | All dependencies check passed")
    return True

def check_command_exists(command: str) -> bool:
    """检查命令是否存在 | Check if command exists"""
    import shutil
    return shutil.which(command) is not None
