"""
🛠️ parabricks WGS分析工具函数模块 | parabricks WGS Analysis Utility Functions Module 🛠️
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

class parabricksLogger:
    """📝 parabricks分析日志管理器 | parabricks Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "parabricks_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 ✍️ | Setup logging"""
        if self.log_file.exists():
            # 备份现有日志文件 💾 | Backup existing log file
            backup_name = f"parabricks_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
            self.log_file.rename(self.output_dir / backup_name)
        
        logging.basicConfig(
            level=logging.INFO,
            format='[%(asctime)s] %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 📢 | Get logger"""
        return self.logger

class CommandRunner:
    """🚀 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 ▶️ | Execute command"""
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
                self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
            return False

class FileProcessor:
    """📁 文件处理器 | File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def find_fastq_files(self):
        """查找FASTQ文件对 🔍 | Find FASTQ file pairs"""
        self.logger.info("🔍 搜索输入文件 | Searching input files")
        
        input_path = Path(self.config.input_dir)
        r1_files = list(input_path.glob(self.config.read1_pattern))
        
        if not r1_files:
            raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到 {self.config.read1_pattern} 文件")
        
        # 按文件名排序 🔠 | Sort by filename
        r1_files.sort()
        
        # 验证R2文件存在 ✅ | Validate R2 files exist
        file_pairs = []
        for r1_file in r1_files:
            # 提取样品名 🏷️ | Extract sample name
            sample_name = r1_file.name.replace("_1.clean.fq.gz", "")
            
            # 构建R2文件路径 🏗️ | Build R2 file path
            r2_file = input_path / f"{sample_name}_2.clean.fq.gz"
            
            if not r2_file.exists():
                self.logger.warning(f"⚠️ 找不到对应的R2文件 | Cannot find corresponding R2 file: {r2_file}")
                continue
            
            file_pairs.append((sample_name, str(r1_file), str(r2_file)))
        
        total_samples = len(file_pairs)
        self.logger.info(f"✅ 找到 {total_samples} 个样品需要处理 | Found {total_samples} samples to process")
        
        return file_pairs
    
    def check_output_exists(self, sample_name: str) -> bool:
        """检查输出文件是否已存在 🧐 | Check if output files already exist"""
        vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
        bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        return vcf_file.exists() and bam_file.exists()
    
    def get_file_size(self, file_path: str) -> str:
        """获取文件大小 📏 | Get file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size_bytes < 1024.0:
                    return f"{size_bytes:.1f} {unit}"
                size_bytes /= 1024.0
            return f"{size_bytes:.1f} TB"
        except:
            return "Unknown"

def check_dependencies(config, logger):
    """检查依赖软件 🧩 | Check dependencies"""
    logger.info("🧩 检查依赖软件 | Checking dependencies")
    
    # 检查parabricks程序 💻 | Check parabricks program
    if not os.path.exists(config.parabricks_path):
        error_msg = f"❌ parabricks程序不存在 | parabricks program does not exist: {config.parabricks_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    # 检查parabricks程序是否可执行 🏃 | Check if parabricks program is executable
    if not os.access(config.parabricks_path, os.X_OK):
        error_msg = f"❌ parabricks程序不可执行 | parabricks program is not executable: {config.parabricks_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    logger.info("✅ ✓ parabricks程序检查通过 | parabricks program check passed")
    return True
