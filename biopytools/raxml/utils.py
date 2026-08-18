"""
RAxML分析工具函数模块|RAxML Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path
from typing import Optional

from ..common.conda_runner import build_conda_command, check_tools

class RAxMLLogger:
    """RAxML分析日志管理器|RAxML Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "raxml_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd, description: str = "", timeout: Optional[int] = None) -> bool:
        """执行命令|Execute command

        cmd 可为列表(shell=False, 由 build_conda_command 构建)或字符串(shell=True)
        |cmd may be a list (shell=False, built by build_conda_command) or a string (shell=True)
        """
        if description:
            self.logger.info(f" 执行步骤|Executing step: {description}")

        if isinstance(cmd, (list, tuple)):
            use_shell = False
            display = ' '.join(str(c) for c in cmd)
        else:
            use_shell = True
            display = str(cmd)
        self.logger.info(f" 命令|Command: {display}")
        self.logger.info(f" 工作目录|Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=use_shell,
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                timeout=timeout
            )
            
            self.logger.info(f" 命令执行成功|Command executed successfully: {description}")
            
            if result.stdout:
                # 只记录关键输出，避免日志过长|Only log key output to avoid overly long logs
                stdout_lines = result.stdout.strip().split('\n')
                if len(stdout_lines) <= 10:
                    self.logger.debug(f" 标准输出|Stdout: {result.stdout}")
                else:
                    self.logger.debug(f" 标准输出|Stdout (前5行|first 5 lines): {chr(10).join(stdout_lines[:5])}...")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f" 命令执行失败|Command execution failed: {description}")
            self.logger.error(f" 错误代码|Error code: {e.returncode}")
            self.logger.error(f" 错误信息|Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f" 标准输出|Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f" 命令执行超时|Command execution timeout: {description}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    missing = check_tools([(config.raxml_path, "RAxML", ["-v"])], logger)
    if missing:
        logger.error(" 请确保RAxML已安装并在PATH中|Please ensure RAxML is installed and in PATH")
        return False
    return True

def get_raxml_version(raxml_path: str) -> str:
    """获取RAxML版本信息|Get RAxML version information"""
    try:
        cmd = build_conda_command(raxml_path, ["-v"])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return result.stdout.strip()
        return "Unknown version"
    except Exception:
        return "Version detection failed"
