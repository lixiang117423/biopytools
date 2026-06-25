"""插入检测工具函数|Insert detection utility functions"""

import logging
import sys
import subprocess
import os
import re
from pathlib import Path
from typing import List, Tuple, Dict
import pandas as pd


class InsertDetectionLogger:
    """插入检测日志管理器|Insert detection logger manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers,
            force=True
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def identify_samples(fastq_dir: str, read1_suffix: str, read2_suffix: str, logger) -> Dict[str, Dict]:
    """识别样品信息和配对关系|Identify sample information and pairing

    Args:
        fastq_dir: FASTQ文件目录|FASTQ file directory
        read1_suffix: R1后缀（包括扩展名，如 _1.fq.gz 或 _1.clean.fq.gz）|R1 suffix with extension
        read2_suffix: R2后缀（包括扩展名，如 _2.fq.gz 或 _2.clean.fq.gz）|R2 suffix with extension
        logger: 日志器|Logger

    Returns:
        dict: 样品信息字典|Sample information dictionary
        {
            'OV8-276': {'r1': 'path/to/OV8-276_1.clean.fq.gz', 'r2': 'path/to/OV8-276_2.clean.fq.gz'},
            'OV8-277': {'r1': 'path/to/OV8-277_1.clean.fq.gz', 'r2': 'path/to/OV8-277_2.clean.fq.gz'}
        }
    """
    fastq_path = Path(fastq_dir)
    # 使用通配符获取所有文件，然后根据后缀过滤|Use wildcard to get all files, then filter by suffix
    fastq_files = list(fastq_path.glob("*"))

    logger.info(f"找到|Found {len(fastq_files)} 个文件|files in directory")
    logger.debug(f"R1后缀|R1 suffix: {read1_suffix}")
    logger.debug(f"R2后缀|R2 suffix: {read2_suffix}")

    samples = {}

    for fq_file in fastq_files:
        filename = fq_file.name

        # 判断是R1还是R2并提取样品名|Determine if R1/R2 and extract sample name
        sample_name = None
        read_type = None

        if filename.endswith(read1_suffix):
            sample_name = filename[:-len(read1_suffix)]
            read_type = 'r1'
        elif filename.endswith(read2_suffix):
            sample_name = filename[:-len(read2_suffix)]
            read_type = 'r2'
        else:
            # 跳过不匹配的文件|Skip non-matching files
            continue

        if not sample_name:
            logger.warning(f"样品名为空，跳过|Empty sample name, skipping: {filename}")
            continue

        file_str = str(fq_file)

        # 初始化或更新样品信息|Initialize or update sample info
        if sample_name not in samples:
            samples[sample_name] = {}

        samples[sample_name][read_type] = file_str

    # 验证配对完整性|Validate pairing completeness
    valid_samples = {}
    for sample_name, files in samples.items():
        if 'r1' in files and 'r2' in files:
            valid_samples[sample_name] = files
            logger.info(f"识别到样品|Identified sample: {sample_name}")
        else:
            logger.warning(f"样品配对不完整，跳过|Sample pairing incomplete, skipping: {sample_name}")
            logger.warning(f"  R1: {files.get('r1', 'Missing')}")
            logger.warning(f"  R2: {files.get('r2', 'Missing')}")

    logger.info(f"有效样品数|Valid samples: {len(valid_samples)}")

    return valid_samples


def check_dependencies(config, logger):
    """检查依赖工具|Check dependency tools

    Args:
        config: InsertDetectionConfig配置对象|InsertDetectionConfig object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有工具都可用|Whether all tools are available
    """
    logger.info("检查依赖工具|Checking dependency tools...")

    tools = {
        'bowtie2': config.bowtie2_path,
        'samtools': config.samtools_path,
    }

    missing = []

    for tool_name, tool_path in tools.items():
        try:
            result = subprocess.run(
                [tool_path, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                logger.debug(f" {tool_name}: 可用|available")
            else:
                missing.append(tool_name)
                logger.warning(f" {tool_name}: 未找到|not found")

        except FileNotFoundError:
            missing.append(tool_name)
            logger.warning(f" {tool_name}: 未找到|not found at {tool_path}")
        except Exception as e:
            logger.warning(f" {tool_name}: 检查失败|check failed: {e}")
            missing.append(tool_name)

    if missing:
        logger.error(f"缺少依赖工具|Missing required tools: {', '.join(missing)}")
        return False

    logger.info("所有依赖工具可用|All dependency tools available")
    return True


def run_command(cmd, logger, description=""):
    """运行命令并记录日志|Run command and log

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        description: 命令描述|Command description

    Returns:
        bool: 是否成功|Whether successful
    """
    try:
        if description:
            logger.info(f"{description}...|{description}...")

        logger.debug(f"运行|Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )

        if description:
            logger.info(f"{description}完成|{description} completed")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"{description}失败|{description} failed")
        if e.stderr:
            logger.error(f"错误信息|Error: {e.stderr}")
        return False

    except Exception as e:
        logger.error(f"命令执行异常|Command execution error: {e}")
        return False
