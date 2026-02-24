"""
BWA比对工具函数模块|BWA Alignment Utility Functions Module
"""

import logging
import subprocess
import sys
import os
import glob
from pathlib import Path
from typing import List, Tuple

class AlignLogger:
    """比对分析日志管理器|Alignment Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "bwament.log"):
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
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = [
        (config.bwa_path, "BWA"),
        (config.samtools_path, "SAMtools")
    ]

    missing_deps = []

    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--version"],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 or name == "BWA":  # BWA returns non-zero for --version
                version_info = result.stdout.strip() or result.stderr.strip()
                logger.info(f"{name} 可用|{name} available: {version_info.split()[0] if version_info else 'installed'}")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)

    if missing_deps:
        error_msg = f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True

def find_fastq_pairs(input_dir: str, pattern: str, logger) -> List[Tuple[str, str, str]]:
    """
    查找FASTQ配对文件|Find FASTQ paired files

    Returns:
        List of tuples: (sample_name, read1_path, read2_path)
    """
    logger.info(f"查找FASTQ文件|Finding FASTQ files")
    logger.info(f"输入目录|Input directory: {input_dir}")
    logger.info(f"匹配模式|Pattern: {pattern}")

    # 构建read1搜索模式|Build read1 search pattern
    search_pattern = os.path.join(input_dir, f"*{pattern}")
    read1_files = glob.glob(search_pattern)

    if not read1_files:
        raise ValueError(f"未找到匹配的FASTQ文件|No FASTQ files found matching: {search_pattern}")

    pairs = []
    for read1 in read1_files:
        # 提取样品名称|Extract sample name
        basename = os.path.basename(read1)
        sample_name = basename.replace(pattern, "")

        # 构建read2路径|Build read2 path
        read2_pattern = pattern.replace("_1.", "_2.")
        read2 = read1.replace(pattern, read2_pattern)

        if not os.path.exists(read2):
            logger.warning(f"未找到配对文件|Paired file not found: {read2}")
            continue

        pairs.append((sample_name, read1, read2))
        logger.info(f"找到样品|Found sample: {sample_name}")

    logger.info(f"共找到 {len(pairs)} 个样品对|Found {len(pairs)} sample pairs")

    return pairs
