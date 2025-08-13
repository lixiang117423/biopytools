"""
BAM统计分析工具函数模块 🔧 | BAM Statistics Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path
from typing import List, Dict, Any

class BAMStatsLogger:
    """BAM统计分析日志管理器 📝 | BAM Statistics Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "bam_stats.log"):
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
    """命令执行器 ⚡ | Command Runner"""
    
    def __init__(self, logger, working_dir: Path, threads: int = 88):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.threads = threads
    
    def run(self, cmd: str, description: str = "", use_threads: bool = True) -> tuple[bool, str]:
        """执行命令 ⚡ | Execute command"""
        if description:
            self.logger.info(f"🔄 执行步骤 | Executing step: {description}")
        
        # 为支持多线程的工具添加线程参数 | Add thread parameters for multi-threading tools
        if use_threads and self.threads > 1:
            if 'samtools' in cmd and '--threads' not in cmd and '-@' not in cmd:
                # 为samtools命令添加线程参数 | Add thread parameter for samtools commands
                if 'samtools view' in cmd:
                    cmd = cmd.replace('samtools view', f'samtools view -@ {self.threads}')
                elif 'samtools sort' in cmd:
                    cmd = cmd.replace('samtools sort', f'samtools sort -@ {self.threads}')
                elif 'samtools index' in cmd:
                    cmd = cmd.replace('samtools index', f'samtools index -@ {self.threads}')
                elif 'samtools depth' in cmd:
                    cmd = cmd.replace('samtools depth', f'samtools depth -@ {self.threads}')
                elif 'samtools flagstat' in cmd:
                    cmd = cmd.replace('samtools flagstat', f'samtools flagstat -@ {self.threads}')
        
        self.logger.debug(f"🖥️ 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.debug(f"✅ 命令执行成功 | Command executed successfully: {description}")
            return True, result.stdout
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            return False, e.stderr

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.samtools_path, "samtools"),
        (config.bedtools_path, "bedtools")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True

def get_sample_name(bam_file: str) -> str:
    """从BAM文件路径提取样品名称 🧪 | Extract sample name from BAM file path"""
    return Path(bam_file).stem
