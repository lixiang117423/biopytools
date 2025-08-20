"""
🔧 工具函数 | Utility Functions
"""

import os
import sys
import logging
import subprocess
from pathlib import Path

class CommandRunner:
    """🚀 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir=None):
        self.logger = logger
        self.working_dir = working_dir or os.getcwd()
    
    def run(self, cmd: str, description: str = "", check: bool = True) -> bool:
        """🚀 执行命令 | Execute command"""
        if description:
            self.logger.info(f"▶️ 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=check,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"⚠️ 标准错误 | Stderr: {result.stderr}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"⚠️ 错误信息 | Error message: {e.stderr}")
            if check:
                raise
            return False


def setup_logger(output_dir: Path, verbose: bool = False) -> logging.Logger:
    """📝 设置日志记录器 | Setup logger"""
    logger = logging.getLogger('kmer_count')
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # 清除现有处理器 | Clear existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 控制台处理器 | Console handler
    # console_handler = logging.StreamHandler()
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # 文件处理器 | File handler
    log_file = output_dir / 'kmer_count.log'
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    
    # 格式设置 | Format setting
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    return logger


def check_dependencies(jellyfish_path: str, logger) -> bool:
    """🔍 检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    try:
        result = subprocess.run([jellyfish_path, "--version"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            logger.info(f"✅ Jellyfish 可用 | available: {version}")
            return True
        else:
            logger.error(f"❌ Jellyfish 不可用 | not available")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"❌ Jellyfish 不可用 | not available: {e}")
        return False

class FileIntegrityError(Exception):
    """💥 文件完整性错误"""
    pass