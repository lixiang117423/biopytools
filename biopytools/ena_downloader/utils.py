"""
ENA下载工具辅助函数模块|ENA Downloader Utilities Module
"""

import os
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Any, List

#  ENA编号格式: 大写字母前缀+纯数字, 覆盖 PRJNA/SRR/ERR/DRR/SAMN 等|ENA accession format: uppercase prefix + digits
ACCESSION_PATTERN = re.compile(r'^[A-Z]{2,6}\d+$')

class DownloadLogger:
    """下载日志管理器|Download Logger Manager"""

    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.log_file = output_path / "ena_download.log"
        self._setup_logging()

    def _setup_logging(self):
        """设置日志配置|Setup logging configuration"""
        # 创建logger|Create logger
        self.logger = logging.getLogger('ena_downloader')
        self.logger.setLevel(logging.DEBUG)

        # 避免重复添加handler|Avoid adding duplicate handlers
        if self.logger.handlers:
            return

        # 创建formatter|Create formatter
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # 控制台handler|Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)

        # 添加handlers|Add handlers
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, output_path: Path):
        self.logger = logger
        self.output_path = output_path

    def run_command(self, cmd: str, shell: bool = True) -> bool:
        """执行系统命令|Execute system command"""
        import subprocess

        self.logger.info(f"执行命令|Executing command: {cmd}")

        try:
            # 改变到输出目录|Change to output directory
            original_cwd = os.getcwd()
            os.chdir(self.output_path)

            result = subprocess.run(
                cmd,
                shell=shell,
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时|1 hour timeout
            )

            # 恢复原目录|Restore original directory
            os.chdir(original_cwd)

            if result.returncode == 0:
                self.logger.info("命令执行成功|Command executed successfully")
                if result.stdout:
                    self.logger.debug(f"标准输出|STDOUT: {result.stdout}")
                return True
            else:
                self.logger.error(f"命令执行失败|Command failed with return code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误输出|STDERR: {result.stderr}")
                return False

        except subprocess.TimeoutExpired:
            os.chdir(original_cwd)
            self.logger.error("命令执行超时|Command execution timed out")
            return False
        except Exception as e:
            os.chdir(original_cwd)
            self.logger.error(f"命令执行异常|Command execution error: {str(e)}")
            return False

def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    dependencies = []

    # 检查wget (对于FTP下载)|Check wget (for FTP downloads)
    if config.protocol == "ftp":
        dependencies.append("wget")

    # 检查ascp (对于Aspera下载)|Check ascp (for Aspera downloads)
    if config.protocol == "aspera":
        dependencies.append("ascp")

    missing_deps = []
    for dep in dependencies:
        if not check_command_exists(dep):
            missing_deps.append(dep)

    if missing_deps:
        logger.error(f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}")
        return False

    logger.info("所有依赖软件检查通过|All dependencies check passed")
    return True

def check_command_exists(command: str) -> bool:
    """检查命令是否存在|Check if command exists"""
    import shutil
    return shutil.which(command) is not None

def is_accession(value: str) -> bool:
    """判断字符串是否为合法ENA编号格式|Check whether a string looks like a valid ENA accession"""
    return bool(ACCESSION_PATTERN.match(value))

def classify_accession(accession: str) -> str:
    """按前缀分类编号层级, 用于日志提示|Classify accession level by prefix for log hints"""
    prefix = accession[:3]
    if prefix in ("PRJ", "SRP", "ERP", "DRP"):
        return "项目级|study"
    if prefix in ("SRR", "ERR", "DRR"):
        return "运行级|run"
    if prefix in ("SRS", "ERS", "DRS") or accession.startswith(("SAMN", "SAME", "SAMD")):
        return "样本级|sample"
    return "未知类型|unknown"

def read_id_file(file_path: Path) -> List[str]:
    """读取ID文件, 每行一个编号, 支持空行和#注释行|Read ID file with one accession per line; blank and # lines skipped"""
    accessions = []
    seen = set()
    with open(file_path, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            entry = line.strip()
            if not entry or entry.startswith('#'):
                continue
            # 同文件内去重, 避免重复下载|Deduplicate within one file to avoid repeated downloads
            if entry in seen:
                continue
            seen.add(entry)
            accessions.append(entry)
    return accessions

def expand_input(raw_value: str) -> List[str]:
    """展开输入为编号列表: 已存在的文件路径按行读取, 否则视为单个编号|Expand input into an accession list: an existing file path is read line by line, otherwise treated as a single accession"""
    candidate = Path(os.path.expanduser(raw_value))
    if candidate.is_file():
        return read_id_file(candidate)
    return [raw_value]
