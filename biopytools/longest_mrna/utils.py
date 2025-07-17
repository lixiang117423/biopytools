"""
最长转录本提取工具函数模块 | Longest mRNA Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
import tempfile
import os
from pathlib import Path

class LongestMRNALogger:
    """最长转录本提取日志管理器 | Longest mRNA Extraction Logger Manager"""
    
    def __init__(self, log_name: str = "longest_mrna_extraction.log"):
        self.log_name = log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
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
    
    def run(self, cmd: str, description: str = "", check: bool = True) -> subprocess.CompletedProcess:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行 | Executing: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                encoding='utf-8',
                check=check
            )
            
            if result.returncode == 0:
                self.logger.info(f"✓ 命令执行成功 | Command executed successfully: {description}")
                if result.stdout.strip():
                    self.logger.debug(f"标准输出 | Stdout: {result.stdout.strip()}")
            else:
                self.logger.error(f"✗ 命令执行失败 | Command execution failed: {description}")
                self.logger.error(f"返回码 | Return code: {result.returncode}")
                if result.stderr.strip():
                    self.logger.error(f"错误信息 | Error: {result.stderr.strip()}")
                if check:
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            
            return result
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            raise

class TempFileManager:
    """临时文件管理器 | Temporary File Manager"""
    
    def __init__(self, logger):
        self.logger = logger
        self.temp_files = []
    
    def create_temp_file(self, mode='w+', delete=False, suffix='', encoding='utf-8'):
        """创建临时文件 | Create temporary file"""
        temp_file = tempfile.NamedTemporaryFile(
            mode=mode, 
            delete=delete, 
            suffix=suffix, 
            encoding=encoding
        )
        self.temp_files.append(temp_file.name)
        return temp_file
    
    def cleanup(self):
        """清理临时文件 | Cleanup temporary files"""
        for temp_file in self.temp_files:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    self.logger.debug(f"清理临时文件 | Cleaned up temp file: {temp_file}")
            except Exception as e:
                self.logger.warning(f"清理临时文件失败 | Failed to cleanup temp file {temp_file}: {e}")
