"""
K-mer分析工具函数模块 | K-mer Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List

class KmerLogger:
    """K-mer分析日志管理器 | K-mer Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "kmer_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
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
        self.working_dir = working_dir
    
    def run(self, cmd: List[str], description: str = "", shell: bool = False) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        if isinstance(cmd, list):
            cmd_str = ' '.join(cmd)
        else:
            cmd_str = cmd
        
        self.logger.info(f"命令 | Command: {cmd_str}")
        
        try:
            if shell:
                result = subprocess.run(
                    cmd_str,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                    cwd=self.working_dir
                )
            else:
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    cwd=self.working_dir
                )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常 | Command execution exception: {e}")
            return False
    
    def run_with_realtime_output(self, cmd: List[str], description: str = "") -> bool:
        """执行命令并实时输出 | Execute command with real-time output"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        self.logger.info(f"命令 | Command: {cmd_str}")
        
        # 调试信息：显示工作目录和文件存在性
        if self.working_dir:
            self.logger.info(f"工作目录 | Working directory: {self.working_dir}")
        
        # 检查关键文件是否存在
        for item in cmd:
            if isinstance(item, str) and (item.endswith('.fof') or 'samples.fof' in item):
                fof_path = Path(item)
                if fof_path.exists():
                    self.logger.info(f"✓ FOF文件存在 | FOF file exists: {fof_path.resolve()}")
                else:
                    self.logger.error(f"✗ FOF文件不存在 | FOF file not found: {fof_path.resolve()}")
                    # 尝试在工作目录中查找
                    if self.working_dir:
                        alt_path = Path(self.working_dir) / fof_path.name
                        if alt_path.exists():
                            self.logger.info(f"✓ 在工作目录中找到FOF文件 | Found FOF file in working directory: {alt_path}")
                        else:
                            self.logger.error(f"✗ 工作目录中也没有FOF文件 | FOF file not in working directory either: {alt_path}")
        
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1,
                cwd=self.working_dir
            )
            
            for line in process.stdout:
                self.logger.info(f"输出 | Output: {line.strip()}")
            
            process.wait()
            
            if process.returncode != 0:
                self.logger.error(f"命令执行失败 | Command execution failed: {description}")
                return False
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            return True
            
        except Exception as e:
            self.logger.error(f"命令执行异常 | Command execution exception: {e}")
            return False

class FileValidator:
    """文件验证器 | File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def check_file_exists(self, file_path: Path, description: str = "") -> bool:
        """检查文件是否存在 | Check if file exists"""
        if file_path.exists():
            if description:
                self.logger.info(f"✓ {description}已存在，跳过 | already exists, skipping: {file_path}")
            return True
        return False
    
    def check_tool_available(self, tool_name: str) -> bool:
        """检查工具是否可用 | Check if tool is available"""
        try:
            result = subprocess.run([tool_name, '--help'], capture_output=True, text=True, timeout=10)
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.logger.error(f"{tool_name}未安装或不在PATH中 | {tool_name} not installed or not in PATH")
            return False

class SequenceUtils:
    """序列处理工具 | Sequence Processing Utilities"""
    
    def __init__(self):
        self.complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    
    def get_reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列 | Get reverse complement sequence"""
        sequence = sequence.upper()
        reverse_sequence = sequence[::-1]
        complement = ''.join(self.complement_dict.get(base, base) for base in reverse_sequence)
        return complement
    
    def get_canonical_kmer(self, kmer: str) -> str:
        """获取canonical k-mer | Get canonical k-mer"""
        rc_kmer = self.get_reverse_complement(kmer)
        return min(kmer, rc_kmer)

def check_dependencies():
    """检查依赖 | Check dependencies"""
    required_packages = [
        ('pyfastx', 'pip install pyfastx'),
        ('rocksdb', 'pip install python-rocksdb'),
    ]
    
    optional_packages = [
        ('sklearn', 'pip install scikit-learn'),
        ('matplotlib', 'pip install matplotlib'),
        ('seaborn', 'pip install seaborn'),
        ('pandas', 'pip install pandas'),
        ('numpy', 'pip install numpy'),
    ]
    
    missing_required = []
    missing_optional = []
    
    for package, install_cmd in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_required.append((package, install_cmd))
    
    for package, install_cmd in optional_packages:
        try:
            __import__(package)
        except ImportError:
            missing_optional.append((package, install_cmd))
    
    if missing_required:
        print("错误 | Error: 缺少必需的依赖包 | Missing required packages:")
        for package, install_cmd in missing_required:
            print(f"  - {package}: {install_cmd}")
        return False
    
    if missing_optional:
        print("警告 | Warning: 缺少可选的依赖包 | Missing optional packages:")
        for package, install_cmd in missing_optional:
            print(f"  - {package}: {install_cmd}")
    
    return True
