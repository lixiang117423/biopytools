"""
Swave工具函数模块|Swave Utility Functions Module
"""

import logging
import sys
import subprocess
import os
from typing import Optional, List, Tuple

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class SwaveLogger:
    """Swave日志管理器|Swave Logger Manager"""

    def __init__(self, output_dir: str, log_file: str = 'swave.log', log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            log_file: 日志文件名|Log file name
            log_level: 日志级别|Log level
        """
        self.log_file = log_file
        self.log_level = log_level
        self.output_dir = output_dir

        # 创建日志目录|Create log directory
        log_dir = os.path.join(output_dir, '99_logs')
        os.makedirs(log_dir, exist_ok=True)
        self.log_path = os.path.join(log_dir, log_file)

        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 清除现有的handlers|Clear existing handlers
        root_logger = logging.getLogger()
        root_logger.handlers.clear()

        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 配置root logger|Configure root logger
        root_logger.setLevel(level)
        root_logger.addHandler(stdout_handler)
        root_logger.addHandler(stderr_handler)
        root_logger.addHandler(file_handler)

        # 禁止传播|Disable propagation
        root_logger.propagate = False

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger
            output_dir: 输出目录|Output directory
        """
        self.logger = logger
        self.output_dir = output_dir

    def run(self, cmd: List[str], description: str = "",
            cwd: Optional[str] = None, env: Optional[dict] = None) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description
            cwd: 工作目录|Working directory
            env: 环境变量|Environment variables

        Returns:
            (成功状态, 标准输出, 标准错误)|(Success status, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        try:
            # 完整命令记录到INFO（规范2.2.1：调试/可复现必需，不能用DEBUG）|Full command at INFO (spec 2.2.1)
            self.logger.info(f"命令|Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=cwd,
                env=env
            )

            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {description}")
                self.logger.error(f"错误输出|Error output: {result.stderr}")
                return False, result.stdout, result.stderr

            if description:
                self.logger.info(f"完成|Completed: {description}")

            return True, result.stdout, result.stderr

        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {description}")
            self.logger.error(f"异常信息|Exception: {str(e)}")
            return False, "", str(e)


def check_dependencies(config, logger: logging.Logger) -> bool:
    """
    检查依赖软件|Check dependencies

    Args:
        config: Swave配置对象|Swave config object
        logger: 日志器|Logger

    Returns:
        是否所有依赖都可用|Whether all dependencies are available
    """
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = []

    # 检查Python环境|Check Python environment
    swave_py = os.path.join(config.swave_path, 'Swave.py')
    if os.path.exists(swave_py):
        logger.info(f"Swave主程序|Swave main script: {swave_py}")
        dependencies.append(True)
    else:
        logger.error(f"Swave主程序不存在|Swave main script not found: {swave_py}")
        dependencies.append(False)

    # 检查minigraph（如果需要）|Check minigraph (if needed)
    if config.gfa_source == 'minigraph':
        minigraph_cmd = build_conda_command(config.minigraph_path, ['--version'])
        logger.info(f"命令|Command: {' '.join(minigraph_cmd)}")
        try:
            result = subprocess.run(minigraph_cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"minigraph可用|minigraph available: {config.minigraph_path}")
                dependencies.append(True)
            else:
                logger.error(f"minigraph不可用|minigraph not available: {config.minigraph_path}")
                dependencies.append(False)
        except Exception as e:
            logger.error(f"minigraph检测失败|minigraph check failed: {e}")
            dependencies.append(False)

    # 检查gfatools（可选）|Check gfatools (optional)
    try:
        gfatools_cmd = build_conda_command(config.gfatools_path, [])
        logger.info(f"命令|Command: {' '.join(gfatools_cmd)}")
        result = subprocess.run(gfatools_cmd, capture_output=True, text=True, timeout=10)
        # gfatools不支持--version，无参数时输出usage（含"gfatools"）并返回非0
        # gfatools doesn't support --version, outputs usage containing "gfatools" without args
        if result.returncode != 127 and 'gfatools' in result.stderr:
            logger.info(f"gfatools可用|gfatools available: {config.gfatools_path}")
        else:
            logger.warning(f"gfatools不可用（可选）|gfatools not available (optional)")
    except Exception:
        logger.warning(f"gfatools检测失败（可选）|gfatools check failed (optional)")

    all_ok = all(dependencies)
    if all_ok:
        logger.info("依赖检查完成|Dependency check completed")
    else:
        logger.error("部分依赖缺失|Some dependencies missing")

    return all_ok
