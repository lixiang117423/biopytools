"""
Bismark流程工具函数模块 | Bismark Pipeline Utility Functions Module
"""
# (此文件无改动 | No changes in this file)
import logging
import subprocess
import sys
import shutil
from pathlib import Path

class BismarkLogger:
    def __init__(self, output_dir: Path):
        self.log_dir = output_dir / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = self.log_dir / f"bismark_pipeline_{timestamp}.log"
        self.logger = logging.getLogger('bismark_pipeline')
        self.logger.setLevel(logging.INFO)
        if not self.logger.handlers:
            formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setFormatter(formatter)
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
            self.logger.addHandler(console_handler)
    
    def get_logger(self):
        return self.logger

class CommandRunner:
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> bool:
        if description: self.logger.info(f"🔧 执行步骤 | Executing step: {description}")
        self.logger.info(f"  命令 | Command: {cmd}")
        self.logger.info(f"  工作目录 | Working directory: {self.working_dir}")
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True, cwd=self.working_dir)
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            if result.stdout: self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"  错误代码 | Error code: {e.returncode}")
            self.logger.error(f"  错误信息 | Error message: {e.stderr}")
            if e.stdout: self.logger.error(f"  标准输出 | Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    logger.info("🔍 检查依赖软件 | Checking dependencies...")
    required = [
        (config.bismark_path, "bismark"),
        (config.bismark_genome_preparation_path, "bismark_genome_preparation"),
        (config.bowtie2_path, "bowtie2"),
        (config.bismark_methylation_extractor_path, "bismark_methylation_extractor"),
    ]
    missing = [name for cmd, name in required if not shutil.which(cmd)]
    if missing:
        raise RuntimeError(f"缺少必需软件 | Missing required software: {', '.join(missing)}")
    logger.info("✅ 依赖软件检查完成 | Dependency check complete.")
    return True
