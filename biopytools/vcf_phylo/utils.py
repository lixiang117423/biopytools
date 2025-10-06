"""
VCF系统发育分析工具函数模块 | VCF Phylogenetic Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class PhyloLogger:
    """系统发育分析日志管理器 | Phylogenetic Analysis Logger Manager"""
    
    def __init__(self, output_prefix: str, log_name: str = "phylo_analysis.log"):
        self.output_prefix = output_prefix
        self.log_file = f"{output_prefix}.log"
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 如果日志文件存在则删除 | Remove log file if exists
        if Path(self.log_file).exists():
            Path(self.log_file).unlink()
        
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
        (config.vcf2dis_path, "VCF2Dis")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 or "usage" in result.stdout.lower() or "usage" in result.stderr.lower():
                logger.info(f"✓ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    # 检查Python依赖 | Check Python dependencies
    python_deps = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("scikit-bio", "skbio")  # 注意：包名是scikit-bio，但导入名是skbio
    ]
    
    python_missing = []
    for dep_name, import_name in python_deps:
        try:
            __import__(import_name)
            logger.info(f"✓ Python包 {dep_name} 可用 | Python package {dep_name} available")
        except ImportError:
            python_missing.append(dep_name)
    
    if missing_deps:
        error_msg = f"缺少外部依赖软件 | Missing external dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    if python_missing:
        error_msg = f"缺少Python依赖包 | Missing Python dependencies: {', '.join(python_missing)}"
        logger.error(error_msg)
        
        # 提供安装建议 | Provide installation suggestions
        logger.error("安装建议 | Installation suggestions:")
        for dep in python_missing:
            if dep == "scikit-bio":
                logger.error(f"  pip install {dep}")
                logger.error(f"  或者 | Or: conda install -c conda-forge {dep}")
                logger.error(f"  或者 | Or: conda install -c bioconda {dep}")
            else:
                logger.error(f"  pip install {dep}")
        
        raise RuntimeError(error_msg)
    
    return True
