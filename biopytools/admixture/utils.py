"""
🛠️ ADMIXTURE分析工具函数模块 | ADMIXURE Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path
from typing import Optional, List

class AdmixtureLogger:
    """📝 ADMIXTURE分析日志管理器 | ADMIXTURE Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "admixture_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """⚙️ 设置日志 | Setup logging"""
        # 🧹 清理旧的处理器，避免日志重复打印
        logger = logging.getLogger()
        if logger.hasHandlers():
            logger.handlers.clear()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, mode='w'), # ✍️ 使用 'w' 模式确保每次运行都是新日志
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """👉 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """⚡️ 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> str:
        """🚀 执行命令 | Execute command"""
        if description:
            self.logger.info(f"▶️ 开始 | Starting: {description}")
        
        cleaned_cmd = " ".join(cmd.strip().split())
        self.logger.info(f"💻 执行命令 | Executing command: {cleaned_cmd}")
        
        try:
            result = subprocess.run(
                cleaned_cmd, 
                shell=True, 
                check=True, 
                capture_output=True, 
                text=True, 
                cwd=self.working_dir
            )
            
            if result.stdout.strip():
                self.logger.info(f"📋 命令输出 | Command output:\n{result.stdout.strip()}")
            
            if result.stderr.strip():
                self.logger.warning(f"⚠️ 命令警告 | Command warnings:\n{result.stderr.strip()}")
            
            self.logger.info(f"✅ 完成 | Completed: {description}")
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"💥 命令执行失败 | Command execution failed: {cleaned_cmd}\n"
                f"   - 🔢 返回码 | Return code: {e.returncode}\n"
                f"   - 📋 标准输出 | Stdout: {e.stdout.strip()}\n"
                f"   - 🚨 标准错误 | Stderr: {e.stderr.strip()}"
            )
            raise

class SoftwareChecker:
    """🔎 软件环境检查器 | Software Environment Checker"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def check_software(self, software: str) -> bool:
        """🧐 检查软件是否安装 | Check if software is installed"""
        if shutil.which(software):
            self.logger.info(f"✅ {software} 已安装 | {software} is installed")
            return True
        else:
            self.logger.warning(f"❌ {software} 未安装 | {software} is not installed")
            return False
    
    def check_dependencies(self) -> bool:
        """🛠️ 检查所有依赖软件 | Check all dependencies"""
        required_software = ["plink", "admixture", "bcftools"]
        optional_software = ["Rscript"]
        
        self.logger.info("=== 🛠️  软件环境检查 | Software Environment Check ===")
        
        all_required = True
        for software in required_software:
            if not self.check_software(software):
                all_required = False
        
        for software in optional_software:
            self.check_software(software)
        
        if not all_required:
            self.logger.error("🛑 缺少必需软件，请安装后再运行 | Missing required software, please install and try again")
            return False
        
        return True

def check_file_format(file_path: str, expected_format: str) -> bool:
    """🧐📄 检查文件格式 | Check file format"""
    if expected_format == "vcf":
        return file_path.endswith(('.vcf', '.vcf.gz'))
    elif expected_format == "bed":
        return file_path.endswith('.bed')
    return False