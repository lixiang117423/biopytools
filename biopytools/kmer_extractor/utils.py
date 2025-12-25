"""
🛠️ K-mer提取工具函数模块 | K-mer Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class KmerLogger:
    """📊 K-mer提取日志管理器 | K-mer Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "kmer_extraction.log"):
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
    """⚡ 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    """🔍 检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = []
    
    # 添加工具依赖
    if hasattr(config, 'file_type'):
        if config.file_type == 'fasta':
            dependencies.append((config.unikmer_path, "Unikmer"))
        elif config.file_type == 'fastq':
            dependencies.append((config.jellyfish_path, "Jellyfish"))
    else:
        # 如果文件类型未确定，检查两个工具
        dependencies.extend([
            (config.unikmer_path, "Unikmer"),
            (config.jellyfish_path, "Jellyfish")
        ])
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            if name == "Unikmer":
                result = subprocess.run([cmd, "version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "Jellyfish":
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
