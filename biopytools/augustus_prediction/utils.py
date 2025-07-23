"""
Augustus流水线工具函数模块 | Augustus Pipeline Utilities Module
"""

import os
import logging
import subprocess
import sys
from pathlib import Path

class PipelineLogger:
    """流水线日志管理器 | Pipeline Logger Manager"""
    
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.log_file = output_path / "augustus_pipeline.log"
        self._setup_logging()
    
    def _setup_logging(self):
        """设置日志配置 | Setup logging configuration"""
        # 创建logger | Create logger
        self.logger = logging.getLogger('augustus_pipeline')
        self.logger.setLevel(logging.DEBUG)
        
        # 避免重复添加handler | Avoid adding duplicate handlers
        if self.logger.handlers:
            return
        
        # 创建formatter | Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # 文件handler | File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # 控制台handler | Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        
        # 添加handlers | Add handlers
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def get_logger(self):
        """获取logger实例 | Get logger instance"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, output_path: Path):
        self.logger = logger
        self.output_path = output_path
    
    def run_command(self, command: str, description: str = "", timeout: int = 3600) -> subprocess.CompletedProcess:
        """执行系统命令 | Execute system command"""
        self.logger.info(f"执行: {description} | Executing: {description}")
        self.logger.debug(f"命令: {command} | Command: {command}")
        
        try:
            # 改变到输出目录 | Change to output directory
            original_cwd = os.getcwd()
            os.chdir(self.output_path)
            
            result = subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8',
                timeout=timeout
            )
            
            # 恢复原目录 | Restore original directory
            os.chdir(original_cwd)
            
            if result.stdout:
                self.logger.debug(f"标准输出: {result.stdout[:500]} | STDOUT: {result.stdout[:500]}")
            
            self.logger.info(f"命令执行成功 | Command executed successfully")
            return result
            
        except subprocess.CalledProcessError as e:
            os.chdir(original_cwd)
            self.logger.error(f"命令执行失败: {e} | Command execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr} | STDERR: {e.stderr}")
            raise
        except subprocess.TimeoutExpired:
            os.chdir(original_cwd)
            self.logger.error(f"命令执行超时 | Command execution timed out")
            raise
        except Exception as e:
            os.chdir(original_cwd)
            self.logger.error(f"命令执行异常: {str(e)} | Command execution error: {str(e)}")
            raise

def check_dependencies(config, logger) -> bool:
    """检查依赖软件 | Check dependencies"""
    dependencies = [
        ('perl', 'Perl interpreter'),
        ('augustus', 'Augustus gene prediction tool')
    ]
    
    missing_deps = []
    
    # 检查基本依赖 | Check basic dependencies
    for dep, desc in dependencies:
        if not check_command_exists(dep):
            missing_deps.append(f"{dep} ({desc})")
    
    # 检查Augustus特定工具 | Check Augustus specific tools
    augustus_tools = [
        'etraining',
        'new_species.pl',
        'gff2gbSmallDNA.pl',
        'randomSplit.pl',
        'optimize_augustus.pl'
    ]
    
    for tool in augustus_tools:
        tool_path = os.path.join(config.augustus_path, tool)
        if not os.path.exists(tool_path):
            missing_deps.append(f"{tool} (Augustus tool at {tool_path})")
    
    if missing_deps:
        logger.error(f"缺少依赖软件: {', '.join(missing_deps)} | Missing dependencies: {', '.join(missing_deps)}")
        return False
    
    logger.info("所有依赖软件检查通过 | All dependencies check passed")
    return True

def check_command_exists(command: str) -> bool:
    """检查命令是否存在 | Check if command exists"""
    import shutil
    return shutil.which(command) is not None
