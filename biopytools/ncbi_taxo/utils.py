"""
NCBI分类学注释工具函数模块|NCBI Taxonomy Annotation Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import shutil
from pathlib import Path
from typing import List, Tuple


class NCBITaxoLogger:
    """NCBI分类学注释日志管理器|NCBI Taxonomy Annotation Logger Manager"""

    def __init__(self, output_prefix: str, log_name: str = "ncbi_taxo.log"):
        self.log_file = f"{output_prefix}.{log_name}"
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除旧日志文件|Delete old log file
        if os.path.exists(self.log_file):
            os.remove(self.log_file)

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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = "."):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> Tuple[bool, str]:
        """执行命令|Execute command

        Returns:
            Tuple[bool, str]: (成功状态, 输出信息)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

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

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

            return True, result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:500]}")
            return False, e.stderr if e.stderr else str(e)


def check_dependencies():
    """检查必需的外部工具|Check required external tools"""
    required_tools = {
        'taxonkit': 'taxonkit',
        'zgrep': 'zgrep',
        'zcat': 'zcat'
    }

    missing_tools = []
    for tool_name, cmd in required_tools.items():
        if not shutil.which(cmd):
            missing_tools.append(tool_name)

    if missing_tools:
        raise RuntimeError(
            f"缺少必需的工具，请先安装|Missing required tools, please install first: {', '.join(missing_tools)}\n"
            f"安装方法|Installation: conda install -c bioconda taxonkit taxonkit"
        )

    return True


def detect_input_type(file_path: str) -> str:
    """自动检测输入文件类型|Auto-detect input file type

    Returns:
        str: 'blast' 或 'accession'
    """
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()

            # 检查是否是制表符分隔的多列文件（BLAST格式）
            fields = first_line.split('\t')
            if len(fields) >= 2:
                # 检查第2列是否像accession（包含字母数字和点）
                if len(fields) > 1:
                    second_col = fields[1].strip()
                    # 简单的accession格式检测
                    if any(c.isalpha() for c in second_col) and any(c.isdigit() for c in second_col):
                        return 'blast'

            # 如果是单列或符合accession格式
            if len(fields) == 1 or (len(fields) >= 1 and any(c.isalpha() for c in fields[0]) and any(c.isdigit() for c in fields[0])):
                return 'accession'

            # 默认返回blast
            return 'blast'

    except Exception as e:
        raise ValueError(f"无法检测输入文件类型|Cannot detect input file type: {e}")


def format_number(num: int) -> str:
    """格式化数字，大数字使用M单位|Format number, use M for large numbers"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def parse_lineage(lineage_str: str, level: str) -> str:
    """从lineage字符串中提取指定层级的名称|Extract name of specified level from lineage string

    Args:
        lineage_str: rank:name格式的lineage字符串，如 "superkingdom:Eukaryota;phylum:Stramenopiles;..."
        level: 目标层级，如 'kingdom', 'phylum', 'genus', 'species' 等

    Returns:
        str: 该层级的名称，如果未找到则返回 'unclassified'
    """
    if not lineage_str or lineage_str.strip() == "":
        return "unclassified"

    parts = lineage_str.split(';')

    # rank别名映射|Rank alias mapping
    # NCBI可能使用不同的rank名称，需要映射到标准名称
    rank_aliases = {
        'kingdom': ['kingdom', 'superkingdom'],
        'phylum': ['phylum'],
        'class': ['class'],
        'order': ['order'],
        'family': ['family'],
        'genus': ['genus'],
        'species': ['species']
    }

    # 获取该层级的所有可能rank名称
    possible_ranks = rank_aliases.get(level, [level])

    # 遍历lineage，查找匹配的rank
    for part in parts:
        part = part.strip()
        if ':' in part:
            rank, name = part.split(':', 1)
            rank = rank.strip().lower()
            name = name.strip()

            # 检查是否匹配目标rank
            for possible_rank in possible_ranks:
                if rank == possible_rank.lower():
                    return name

    # 如果没找到，返回unclassified
    return "unclassified"
