# ===== FILE: orthofinder_pangenome/utils.py =====
"""
OrthoFinder泛基因组分析工具函数模块 | OrthoFinder Pangenome Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Dict, Tuple

class PangenomeLogger:
    """泛基因组分析日志管理器 | Pangenome Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "pangenome_analysis.log"):
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
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        self.logger.info(f"工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
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
            self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.orthofinder_path, "OrthoFinder")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            # 检查文件是否存在且可执行
            if not os.path.isfile(cmd):
                logger.error(f"文件不存在 | File does not exist: {cmd}")
                missing_deps.append(name)
                continue
                
            if not os.access(cmd, os.X_OK):
                logger.error(f"文件不可执行 | File not executable: {cmd}")
                missing_deps.append(name)
                continue
            
            # 使用绝对路径和shell=True执行
            result = subprocess.run(f'"{cmd}" -h', shell=True,
                                  capture_output=True, text=True, timeout=10)
            
            # 检查输出中是否包含OrthoFinder版本信息
            if "OrthoFinder version" in result.stdout or "OrthoFinder version" in result.stderr:
                logger.info(f"{name} 可用 | available")
            else:
                logger.error(f"程序响应异常 | Program response abnormal: {cmd}")
                logger.error(f"标准输出: {result.stdout}")
                logger.error(f"标准错误: {result.stderr}")
                missing_deps.append(name)
                
        except subprocess.TimeoutExpired:
            logger.error(f"程序响应超时 | Program response timeout: {cmd}")
            missing_deps.append(name)
        except Exception as e:
            logger.error(f"检查程序时出错 | Error checking program {cmd}: {e}")
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True

def count_fasta_files(input_dir: str) -> Tuple[int, List[str]]:
    """统计FASTA文件数量和列表 | Count FASTA files and get file list"""
    input_path = Path(input_dir)
    
    # 支持的FASTA文件扩展名 | Supported FASTA file extensions
    fasta_extensions = ['*.fa', '*.faa', '*.fas', '*.fasta', '*.pep', '*.protein']
    
    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(input_path.glob(ext))
    
    # 获取文件名（不含扩展名）作为样本名 | Get filenames without extensions as sample names
    sample_names = [f.stem for f in fasta_files]
    
    return len(fasta_files), sample_names

def validate_fasta_files(input_dir: str, logger) -> bool:
    """验证FASTA文件格式 | Validate FASTA file format"""
    logger.info("验证输入文件格式 | Validating input file format")
    
    count, sample_names = count_fasta_files(input_dir)
    
    if count == 0:
        logger.error("未找到FASTA文件 | No FASTA files found")
        return False
    
    logger.info(f"找到 {count} 个FASTA文件 | Found {count} FASTA files")
    logger.info(f"样本名称 | Sample names: {', '.join(sample_names[:5])}{'...' if count > 5 else ''}")
    
    # 检查文件是否为空 | Check if files are empty
    input_path = Path(input_dir)
    empty_files = []
    
    for sample_name in sample_names:
        # 查找对应的文件 | Find corresponding file
        possible_files = list(input_path.glob(f"{sample_name}.*"))
        if possible_files:
            file_path = possible_files[0]
            if file_path.stat().st_size == 0:
                empty_files.append(sample_name)
    
    if empty_files:
        logger.warning(f"发现空文件 | Found empty files: {', '.join(empty_files)}")
    
    return True

# ===== END FILE =====