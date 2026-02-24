"""
BAM统计分析工具函数模块|BAM Statistics Analysis Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import shutil
from pathlib import Path
from typing import List, Dict, Any

class BAMStatsLogger:
    """BAM统计分析日志管理器|BAM Statistics Analysis Logger Manager"""

    def __init__(self, output_dir, log_name: str = "bam_stats.log"):
        # 确保output_dir是Path对象|Ensure output_dir is Path object
        self.output_dir = Path(output_dir).expanduser() if not isinstance(output_dir, Path) else output_dir.expanduser()
        self.log_file = self.output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 确保输出目录存在|Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if self.log_file.exists():
            self.log_file.unlink()

        # 配置日志格式|Configure logging format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            datefmt=date_format,
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path, threads: int = 88):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.threads = threads

    def run(self, cmd: str, description: str = "", use_threads: bool = True) -> tuple[bool, str]:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 为支持多线程的工具添加线程参数|Add thread parameters for multi-threading tools
        if use_threads and self.threads > 1:
            if 'samtools' in cmd and '--threads' not in cmd and '-@' not in cmd:
                # 为samtools命令添加线程参数|Add thread parameter for samtools commands
                if 'samtools view' in cmd:
                    cmd = cmd.replace('samtools view', f'samtools view -@ {self.threads}')
                elif 'samtools sort' in cmd:
                    cmd = cmd.replace('samtools sort', f'samtools sort -@ {self.threads}')
                elif 'samtools index' in cmd:
                    cmd = cmd.replace('samtools index', f'samtools index -@ {self.threads}')
                elif 'samtools depth' in cmd:
                    cmd = cmd.replace('samtools depth', f'samtools depth -@ {self.threads}')
                elif 'samtools flagstat' in cmd:
                    cmd = cmd.replace('samtools flagstat', f'samtools flagstat -@ {self.threads}')

        self.logger.debug(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.debug(f"命令执行成功|Command executed successfully: {description}")
            return True, result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False, e.stderr

def check_dependencies(logger=None):
    """检查依赖软件是否已安装|Check if required software is installed"""

    if logger is None:
        # 创建临时logger|Create temporary logger
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

    # 检查Python库|Check Python libraries
    missing_libs = []
    try:
        import pandas
    except ImportError:
        missing_libs.append('pandas')
    try:
        import openpyxl
    except ImportError:
        missing_libs.append('openpyxl')
    try:
        import tqdm
    except ImportError:
        missing_libs.append('tqdm')

    if missing_libs:
        logger.error(f"错误:缺少必要的Python库|Error:Missing required Python libraries. 请安装|Please install: {' '.join(missing_libs)}")
        logger.info(f"   运行|Run: pip install {' '.join(missing_libs)}")
        return False

    # 检查samtools|Check samtools
    try:
        subprocess.run(['samtools', '--version'], check=True, capture_output=True)
        logger.info("samtools已找到，准备开始分析|samtools found, ready to start analysis...")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("错误:未找到'samtools'命令|Error:'samtools' command not found")
        logger.info("   请确保samtools已安装并且在您的系统PATH环境变量中|Please ensure samtools is installed and in your PATH")
        logger.info("   您可以通过'conda install -c bioconda samtools'进行安装|You can install with: conda install -c bioconda samtools")
        return False

def get_sample_name(bam_file: str) -> str:
    """从BAM文件路径提取样品名称|Extract sample name from BAM file path"""
    return Path(bam_file).stem
