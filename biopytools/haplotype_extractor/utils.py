"""
VCF单体型提取工具函数模块 | VCF Haplotype Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
import tempfile
import os
from pathlib import Path

class HaplotypeLogger:
    """单体型分析日志管理器 | Haplotype Analysis Logger Manager"""
    
    def __init__(self, output_file: str, log_name: str = "haplotype_extraction.log"):
        self.output_dir = Path(output_file).parent
        self.log_file = self.output_dir / log_name
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
    
    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir.resolve() if working_dir else Path.cwd()
    
    def run(self, cmd: list, description: str = "") -> tuple:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.debug(f"命令 | Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.debug(f"命令执行成功 | Command executed successfully: {description}")
            return True, result.stdout, result.stderr
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False, e.stdout, e.stderr
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到 | Command not found: {cmd[0]}")
            return False, "", str(e)

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    try:
        result = subprocess.run(
            [config.bcftools_path, "--version"], 
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            logger.info(f"✓ bcftools 可用 | bcftools available")
            return True
        else:
            raise RuntimeError("bcftools不可用 | bcftools not available")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        error_msg = "缺少依赖软件: bcftools | Missing dependency: bcftools"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

class TempFileManager:
    """临时文件管理器 | Temporary File Manager"""
    
    def __init__(self, logger):
        self.logger = logger
        self.temp_dir = tempfile.mkdtemp()
        self.temp_files = []
    
    def create_temp_file(self, suffix: str = "") -> str:
        """创建临时文件 | Create temporary file"""
        temp_file = os.path.join(self.temp_dir, f"temp_{len(self.temp_files)}{suffix}")
        self.temp_files.append(temp_file)
        return temp_file
    
    def cleanup(self):
        """清理临时文件 | Cleanup temporary files"""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                self.logger.debug(f"删除临时文件 | Deleted temporary file: {temp_file}")
        
        if os.path.exists(self.temp_dir):
            os.rmdir(self.temp_dir)
            self.logger.debug(f"删除临时目录 | Deleted temporary directory: {self.temp_dir}")
