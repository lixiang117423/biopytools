"""
MSA工具函数模块 | MSA Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class MSALogger:
    """MSA日志管理器 | MSA Logger Manager"""
    
    def __init__(self, output_prefix: str, log_name: str = "msa_analysis.log"):
        self.log_file = Path(f"{output_prefix}.log")
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
    
    def __init__(self, logger):
        self.logger = logger
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🔄 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    # 根据选择的方法检查对应工具 | Check tool based on selected method
    tool_map = {
        'mafft': config.mafft_path,
        'clustalo': config.clustalo_path,
        'muscle': config.muscle_path,
        't_coffee': config.tcoffee_path
    }
    
    tool_cmd = tool_map.get(config.method)
    
    try:
        # 不同工具的版本命令不同 | Different tools have different version commands
        if config.method == 'mafft':
            result = subprocess.run([tool_cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
        elif config.method == 'clustalo':
            result = subprocess.run([tool_cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
        elif config.method == 'muscle':
            result = subprocess.run([tool_cmd, "-version"], 
                                  capture_output=True, text=True, timeout=10)
        elif config.method == 't_coffee':
            result = subprocess.run([tool_cmd, "-version"], 
                                  capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0 or config.method in ['mafft', 'clustalo']:
            logger.info(f"✅ {config.method.upper()} 可用 | available")
            return True
        else:
            raise FileNotFoundError
            
    except (subprocess.TimeoutExpired, FileNotFoundError):
        error_msg = f"❌ 缺少依赖软件 | Missing dependency: {config.method.upper()}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

def count_sequences(fasta_file: str) -> int:
    """统计FASTA文件中的序列数 | Count sequences in FASTA file"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count
