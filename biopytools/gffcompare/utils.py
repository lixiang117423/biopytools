"""
GFFcompare两两比较工具函数模块|GFFcompare Pairwise Comparison Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Tuple, Optional

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class GffCompareLogger:
    """GFFcompare分析日志管理器|GFFcompare Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "gffcompare.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / "99_logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logger = logging.getLogger("gffcompare")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()

    def run(self, cmd: list, description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False, e.stdout, e.stderr


def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        cmd = build_conda_command(config.gffcompare_path, ['--version'])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)

        if result.returncode == 0:
            version_info = result.stdout.strip()
            logger.info(f"gffcompare 可用|gffcompare available: {version_info}")
        else:
            raise RuntimeError(f"gffcompare版本检查失败|gffcompare version check failed: {result.stderr}")

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"gffcompare不可用|gffcompare not available: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True
