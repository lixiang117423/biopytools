"""
INDEL PAV分析工具函数模块 | INDEL PAV Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import gzip
from pathlib import Path
from typing import Union

class PAVLogger:
    """PAV分析日志管理器 | PAV Analysis Logger Manager"""
    
    def __init__(self, output_path: Path, log_name: str = "indel_pav_analysis.log"):
        self.output_path = output_path
        self.log_file = output_path / log_name
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

class FileHandler:
    """文件处理器 | File Handler"""
    
    @staticmethod
    def is_compressed(file_path: str) -> bool:
        """检查文件是否压缩 | Check if file is compressed"""
        return file_path.endswith('.gz')
    
    @staticmethod
    def open_file(file_path: str, mode: str = 'r'):
        """智能打开文件 | Smart file opening"""
        if FileHandler.is_compressed(file_path):
            if 'b' not in mode:
                mode = mode.replace('r', 'rt').replace('w', 'wt')
            return gzip.open(file_path, mode)
        else:
            return open(file_path, mode)
    
    @staticmethod
    def count_samples(vcf_file: str) -> int:
        """统计VCF文件中的样本数 | Count samples in VCF file"""
        with FileHandler.open_file(vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # VCF格式：前9列是固定的，从第10列开始是样本
                    return len(line.strip().split('\t')) - 9
        return 0

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"📝 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.bcftools_path, "BCFtools")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✅ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True
