"""
基因组组装工具函数模块 | Genome Assembly Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class AssemblyLogger:
    """基因组组装日志管理器 | Genome Assembly Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "assembly_analysis.log"):
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
    
    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
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
                cwd=self.working_dir,
                timeout=timeout
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时 | Command execution timeout: {description}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    # 基本依赖 | Basic dependencies
    basic_deps = [
        (config.verkko_path, "Verkko"),
        (config.hifiasm_path, "hifiasm"),
        (config.minimap2_path, "minimap2"),
        (config.mashmap_path, "mashmap")
    ]
    
    # 可选依赖 | Optional dependencies
    optional_deps = []
    if config.run_contamination_screen:
        optional_deps.append((config.fcs_path, "NCBI FCS"))
    if config.run_flagger:
        optional_deps.append((config.flagger_path, "Flagger"))
    if config.run_merqury:
        optional_deps.append((config.merqury_path, "Merqury"))
    if config.run_inspector:
        optional_deps.append((config.inspector_path, "Inspector"))
    if config.run_deepvariant:
        optional_deps.append((config.deepvariant_path, "DeepVariant"))
    if config.run_compleasm:
        optional_deps.append((config.compleasm_path, "compleasm"))
    
    missing_deps = []
    optional_missing = []
    
    # 检查基本依赖 | Check basic dependencies
    for cmd, name in basic_deps:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    # 检查可选依赖 | Check optional dependencies
    for cmd, name in optional_deps:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                optional_missing.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            optional_missing.append(name)
    
    if missing_deps:
        error_msg = f"缺少必需依赖软件 | Missing required dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    if optional_missing:
        logger.warning(f"缺少可选依赖软件 | Missing optional dependencies: {', '.join(optional_missing)}")
        logger.warning("相关功能将被跳过 | Related functions will be skipped")
    
    return True
