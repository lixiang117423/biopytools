"""
HiFiasm工具模块 | HiFiasm Utilities Module
"""

import os
import sys
import logging
import subprocess
import shutil
import time
import signal
from pathlib import Path
from typing import Optional, List, Dict, Any, Union
from datetime import datetime

class HifiasmLogger:
    """HiFiasm日志管理器 | HiFiasm Logger Manager"""
    
    def __init__(self, output_dir: str = '.', log_level: str = 'INFO', verbose: int = 0):
        self.output_dir = Path(output_dir)
        self.log_dir = self.output_dir / 'logs'
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置日志级别
        if verbose >= 3:
            log_level = 'DEBUG'
        elif verbose >= 2:
            log_level = 'INFO'
        elif verbose >= 1:
            log_level = 'WARNING'
        
        self.log_level = getattr(logging, log_level.upper())
        
        # 创建主日志器
        self.base_logger = self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:
        """设置日志器 | Setup logger"""
        logger = logging.getLogger('hifiasm')
        logger.setLevel(self.log_level)
        
        # 清除已有处理器
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
        
        # 创建格式化器
        formatter = logging.Formatter(
            fmt='%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # 文件处理器
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = self.log_dir / f'hifiasm_{timestamp}.log'
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # 控制台处理器
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.log_level)
        
        # 彩色格式化器（如果支持）
        if self._supports_color():
            console_formatter = ColoredFormatter()
        else:
            console_formatter = formatter
        console_handler.setFormatter(console_formatter)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        
        return logger
    
    def _supports_color(self) -> bool:
        """检查是否支持彩色输出 | Check if color output is supported"""
        return (
            hasattr(sys.stdout, 'isatty') and sys.stdout.isatty() and
            os.environ.get('TERM') != 'dumb'
        )
    
    def get_logger(self) -> 'EnhancedLogger':
        """获取增强日志器 | Get enhanced logger"""
        return EnhancedLogger(self.base_logger)

class EnhancedLogger:
    """增强的日志器类 | Enhanced Logger Class"""
    
    def __init__(self, base_logger: logging.Logger):
        self.base_logger = base_logger
    
    def start(self, message: str):
        """开始步骤日志 | Start step log"""
        self.base_logger.info(f"🚀 开始 | Starting: {message}")
    
    def success(self, message: str):
        """成功日志 | Success log"""
        self.base_logger.info(f"✅ 成功 | Success: {message}")
    
    def error(self, message: str):
        """错误日志 | Error log"""
        self.base_logger.error(f"❌ 错误 | Error: {message}")
    
    def warning(self, message: str):
        """警告日志 | Warning log"""
        self.base_logger.warning(f"⚠️  警告 | Warning: {message}")
    
    def info(self, message: str):
        """信息日志 | Info log"""
        self.base_logger.info(f"ℹ️  信息 | Info: {message}")
    
    def debug(self, message: str):
        """调试日志 | Debug log"""
        self.base_logger.debug(f"🐛 调试 | Debug: {message}")
    
    def progress(self, message: str):
        """进度日志 | Progress log"""
        self.base_logger.info(f"📊 进度 | Progress: {message}")
    
    # 代理标准Logger方法
    def __getattr__(self, name):
        """代理其他方法到基础logger | Proxy other methods to base logger"""
        return getattr(self.base_logger, name)

class ColoredFormatter(logging.Formatter):
    """彩色日志格式化器 | Colored Log Formatter"""
    
    COLORS = {
        'DEBUG': '\033[36m',     # 青色
        'INFO': '\033[32m',      # 绿色
        'WARNING': '\033[33m',   # 黄色
        'ERROR': '\033[31m',     # 红色
        'CRITICAL': '\033[35m',  # 紫色
    }
    RESET = '\033[0m'
    
    def format(self, record):
        log_color = self.COLORS.get(record.levelname, self.RESET)
        record.levelname = f"{log_color}{record.levelname}{self.RESET}"
        
        formatter = logging.Formatter(
            fmt='%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        return formatter.format(record)

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger: logging.Logger, working_dir: str, timeout: Optional[int] = None):
        self.logger = logger
        self.working_dir = Path(working_dir).resolve()
        self.timeout = timeout
        self.processes = []  # 跟踪运行的进程
    
    def run(self, cmd: Union[str, List[str]], description: str = "", 
            check: bool = True, capture_output: bool = True, 
            timeout: Optional[int] = None) -> subprocess.CompletedProcess:
        """执行命令 | Execute command"""
        
        if isinstance(cmd, list):
            cmd_str = ' '.join(str(c) for c in cmd)
            cmd_list = cmd
        else:
            cmd_str = cmd
            cmd_list = cmd.split()
        
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd_str}")
        self.logger.info(f"工作目录 | Working directory: {self.working_dir}")
        
        effective_timeout = timeout or self.timeout
        
        try:
            # 创建进程
            process = subprocess.Popen(
                cmd_list,
                cwd=str(self.working_dir),
                stdout=subprocess.PIPE if capture_output else None,
                stderr=subprocess.PIPE if capture_output else None,
                text=True,
                encoding='utf-8'
            )
            
            self.processes.append(process)
            
            try:
                stdout, stderr = process.communicate(timeout=effective_timeout)
                return_code = process.returncode
            except subprocess.TimeoutExpired:
                self.logger.error(f"命令超时 | Command timeout: {description}")
                process.kill()
                stdout, stderr = process.communicate()
                return_code = -1
            finally:
                if process in self.processes:
                    self.processes.remove(process)
            
            # 创建返回对象
            result = subprocess.CompletedProcess(
                args=cmd_list,
                returncode=return_code,
                stdout=stdout,
                stderr=stderr
            )
            
            if return_code == 0:
                self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
                if stdout and self.logger.level <= logging.DEBUG:
                    self.logger.debug(f"标准输出 | Stdout:\n{stdout}")
            else:
                error_msg = f"命令执行失败 | Command execution failed: {description}"
                self.logger.error(error_msg)
                self.logger.error(f"返回码 | Return code: {return_code}")
                if stderr:
                    self.logger.error(f"错误信息 | Error message:\n{stderr}")
                if stdout:
                    self.logger.error(f"标准输出 | Stdout:\n{stdout}")
                
                if check:
                    raise subprocess.CalledProcessError(return_code, cmd_list, stdout, stderr)
            
            return result
            
        except FileNotFoundError as e:
            error_msg = f"命令未找到 | Command not found: {cmd_list[0]}"
            self.logger.error(error_msg)
            if check:
                raise RuntimeError(error_msg) from e
            return subprocess.CompletedProcess(cmd_list, 127, "", str(e))
        
        except Exception as e:
            error_msg = f"命令执行异常 | Command execution exception: {e}"
            self.logger.error(error_msg)
            if check:
                raise
            return subprocess.CompletedProcess(cmd_list, 1, "", str(e))
    
    def run_with_progress(self, cmd: Union[str, List[str]], description: str = "",
                         progress_pattern: Optional[str] = None) -> subprocess.CompletedProcess:
        """执行带进度监控的命令 | Execute command with progress monitoring"""
        
        if isinstance(cmd, list):
            cmd_str = ' '.join(str(c) for c in cmd)
            cmd_list = cmd
        else:
            cmd_str = cmd
            cmd_list = cmd.split()
        
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd_str}")
        
        try:
            process = subprocess.Popen(
                cmd_list,
                cwd=str(self.working_dir),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding='utf-8',
                bufsize=1,
                universal_newlines=True
            )
            
            self.processes.append(process)
            
            output_lines = []
            
            try:
                for line in iter(process.stdout.readline, ''):
                    line = line.rstrip()
                    if line:
                        output_lines.append(line)
                        self.logger.debug(f"输出 | Output: {line}")
                        
                        # 检查进度模式
                        if progress_pattern and progress_pattern in line:
                            self.logger.info(f"进度 | Progress: {line}")
                
                process.stdout.close()
                return_code = process.wait()
                
            finally:
                if process in self.processes:
                    self.processes.remove(process)
            
            output = '\n'.join(output_lines)
            
            result = subprocess.CompletedProcess(
                args=cmd_list,
                returncode=return_code,
                stdout=output,
                stderr=""
            )
            
            if return_code == 0:
                self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            else:
                self.logger.error(f"命令执行失败 | Command execution failed: {description}")
                self.logger.error(f"返回码 | Return code: {return_code}")
                raise subprocess.CalledProcessError(return_code, cmd_list, output, "")
            
            return result
            
        except Exception as e:
            self.logger.error(f"命令执行异常 | Command execution exception: {e}")
            raise
    
    def cleanup(self):
        """清理进程 | Cleanup processes"""
        for process in self.processes[:]:
            if process.poll() is None:
                self.logger.warning(f"终止进程 | Terminating process: {process.pid}")
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.logger.warning(f"强制杀死进程 | Force killing process: {process.pid}")
                    process.kill()

def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.hifiasm_path, "HiFiasm", ["--version"]),
        (config.python_path, "Python3", ["--version"]),
    ]
    
    if not config.skip_busco:
        dependencies.append((config.busco_path, "BUSCO", ["--version"]))
    
    if not config.skip_quast:
        dependencies.append((config.quast_path, "QUAST", ["--version"]))
    
    missing_deps = []
    
    for cmd, name, version_args in dependencies:
        try:
            result = subprocess.run(
                [cmd] + version_args,
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                version_info = result.stdout.strip() or result.stderr.strip()
                logger.info(f"✓ {name} 可用 | available: {version_info.split()[0] if version_info else 'unknown version'}")
            else:
                missing_deps.append(name)
                logger.error(f"✗ {name} 版本检查失败 | version check failed")
                
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            missing_deps.append(name)
            logger.error(f"✗ {name} 未找到或不可执行 | not found or not executable: {cmd}")
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        return False
    
    logger.info("✅ 所有依赖软件检查通过 | All dependencies check passed")
    return True

def estimate_resources(input_file: str, logger: logging.Logger) -> Dict[str, Any]:
    """估计资源需求 | Estimate resource requirements"""
    logger.info("估计资源需求 | Estimating resource requirements")
    
    try:
        file_size = os.path.getsize(input_file)
        file_size_gb = file_size / (1024**3)
        
        # 基于文件大小估计资源需求
        estimated_genome_size_gb = file_size_gb / 50  # 假设50x覆盖度
        
        # 估计内存需求（通常需要基因组大小的10-20倍）
        estimated_memory_gb = max(64, int(estimated_genome_size_gb * 15))
        
        # 估计运行时间（小时）
        estimated_runtime_hours = max(8, int(file_size_gb / 10))
        
        resources = {
            'input_file_size_gb': round(file_size_gb, 2),
            'estimated_genome_size_gb': round(estimated_genome_size_gb, 2),
            'recommended_memory_gb': estimated_memory_gb,
            'estimated_runtime_hours': estimated_runtime_hours,
            'recommended_threads': min(64, max(16, int(file_size_gb / 5)))
        }
        
        logger.info(f"输入文件大小 | Input file size: {resources['input_file_size_gb']} GB")
        logger.info(f"估计基因组大小 | Estimated genome size: {resources['estimated_genome_size_gb']} GB")
        logger.info(f"推荐内存 | Recommended memory: {resources['recommended_memory_gb']} GB")
        logger.info(f"估计运行时间 | Estimated runtime: {resources['estimated_runtime_hours']} hours")
        logger.info(f"推荐线程数 | Recommended threads: {resources['recommended_threads']}")
        
        return resources
        
    except Exception as e:
        logger.warning(f"资源估计失败 | Resource estimation failed: {e}")
        return {
            'input_file_size_gb': 0,
            'estimated_genome_size_gb': 1.4,
            'recommended_memory_gb': 64,
            'estimated_runtime_hours': 12,
            'recommended_threads': 32
        }

def setup_signal_handlers(cmd_runner: CommandRunner, logger: logging.Logger):
    """设置信号处理器 | Setup signal handlers"""
    def signal_handler(signum, frame):
        logger.warning(f"接收到信号 | Received signal: {signum}")
        logger.info("正在清理进程 | Cleaning up processes...")
        cmd_runner.cleanup()
        logger.info("进程清理完成 | Process cleanup completed")
        sys.exit(1)
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

def format_duration(seconds: float) -> str:
    """格式化持续时间 | Format duration"""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = int(seconds % 60)
    
    if hours > 0:
        return f"{hours}小时{minutes}分钟{seconds}秒 | {hours}h{minutes}m{seconds}s"
    elif minutes > 0:
        return f"{minutes}分钟{seconds}秒 | {minutes}m{seconds}s"
    else:
        return f"{seconds}秒 | {seconds}s"

def format_file_size(size_bytes: int) -> str:
    """格式化文件大小 | Format file size"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"

def create_directory_structure(output_dir: str) -> Dict[str, str]:
    """创建目录结构 | Create directory structure"""
    base_dir = Path(output_dir)
    
    directories = {
        'base': str(base_dir),
        'logs': str(base_dir / 'logs'),
        'tmp': str(base_dir / 'tmp'),
        'assembly': str(base_dir / 'assembly'),
        'quality_assessment': str(base_dir / 'quality_assessment'),
        'statistics': str(base_dir / 'statistics'),
        'plots': str(base_dir / 'plots'),
        'final_results': str(base_dir / 'final_results')
    }
    
    for dir_path in directories.values():
        os.makedirs(dir_path, exist_ok=True)
    
    return directories