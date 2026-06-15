"""
覆盖度过滤工具函数模块|Coverage Filter Utility Functions Module
"""

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path


class CoverageFilterLogger:
    """覆盖度过滤日志管理器|Coverage Filter Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "coverage_filter.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除旧日志文件|Delete old log file
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies(logger):
    """检查依赖软件|Check dependencies"""
    dependencies = ['samtools', 'seqtk', 'seqkit']
    missing = []

    for cmd in dependencies:
        # 使用 shutil.which 查找命令|Use shutil.which to find command
        cmd_path = shutil.which(cmd)
        if cmd_path:
            logger.debug(f"依赖检查通过|Dependency check passed: {cmd} ({cmd_path})")
        else:
            logger.error(f"依赖软件未找到|Dependency not found: {cmd}")
            missing.append(cmd)

    if missing:
        logger.error(f"缺少依赖|Missing dependencies: {', '.join(missing)}")
        logger.error("请安装缺少的软件或确保它们在PATH中|Please install missing software or ensure they are in PATH")
        return False

    logger.info("所有依赖检查通过|All dependencies check passed")
    return True


def run_command(cmd: str, logger, description: str = "") -> bool:
    """执行命令|Execute command"""
    if description:
        logger.info(f"执行步骤|Executing step: {description}")

    logger.debug(f"命令|Command: {cmd}")

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )

        logger.debug(f"命令执行成功|Command executed successfully")

        if result.stdout:
            logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")  # 只记录前500字符

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {description}")
        logger.error(f"错误代码|Error code: {e.returncode}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False
