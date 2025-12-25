"""
🔧 系统发育树构建工具函数模块 | Phylogenetic Tree Builder Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

class PhyloLogger:
    """系统发育树分析日志管理器 | Phylogenetic Tree Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "mafft_fasttree.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.DEBUG,  # 使用DEBUG级别输出详细日志
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
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
    
    def run(self, cmd: str, description: str = "", allow_fail: bool = False) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"📝 命令 | Command: {cmd}")
        self.logger.debug(f"📂 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=not allow_fail,
                cwd=self.working_dir
            )
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout:\n{result.stdout}")
            
            if result.stderr:
                # 有些程序会把正常信息输出到stderr
                if result.returncode == 0:
                    self.logger.debug(f"📋 标准错误 | Stderr:\n{result.stderr}")
                else:
                    self.logger.error(f"❌ 错误信息 | Error message:\n{result.stderr}")
            
            if result.returncode == 0:
                self.logger.info(f"✅ 命令执行成功 | Command executed successfully")
                return True
            elif allow_fail:
                self.logger.warning(f"⚠️ 命令执行失败但已忽略 | Command failed but ignored")
                return False
            else:
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"❌ 错误代码 | Error code: {e.returncode}")
            if not allow_fail:
                raise
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.mafft_path, "MAFFT"),
        (config.fasttree_path, "FastTree")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            # MAFFT使用--version
            if 'mafft' in cmd.lower():
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            # FastTree使用-expert然后立即退出
            else:
                result = subprocess.run([cmd], 
                                      capture_output=True, text=True, timeout=10, input="\n")
            
            logger.info(f"✅ {name} 可用 | available")
            
        except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
            logger.error(f"❌ {name} 不可用 | not available: {e}")
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True
