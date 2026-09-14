"""a-liner pipeline 工具模块|aliner pipeline utility module"""

import os
import re
import sys
import logging
import subprocess
from typing import List, Optional, Tuple, Dict

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import conda_env_run_prefix  # §13 同源conda绝对路径+run -p


def parse_seq_spec(spec: str) -> Tuple[str, Optional[int], Optional[int]]:
    """
    解析序列规格|Parse sequence spec

    格式|format:
        'chrZ'            -> ('chrZ', None, None)        整条|full length
        'chrZ:1-30000000' -> ('chrZ', 1, 30000000)       1-based 区段|1-based region
    """
    spec = spec.strip()
    if not spec:
        raise ValueError(f"非法序列规格|Invalid sequence spec: '{spec}'。"
                         f"格式|format: 'chrZ' 或|or 'chrZ:start-end'")
    m = re.match(r'^([^:]+)(?::(\d+)-(\d+))?$', spec)
    if not m:
        raise ValueError(f"非法序列规格|Invalid sequence spec: '{spec}'。"
                         f"格式|format: 'chrZ' 或|or 'chrZ:start-end'")
    seq_id = m.group(1)
    if m.group(2) is not None:
        start, end = int(m.group(2)), int(m.group(3))
        if start < 1 or end < start:
            raise ValueError(f"非法区段|Invalid region: '{spec}'（start>=1 且|and start<=end）")
        return (seq_id, start, end)
    return (seq_id, None, None)


class AlinerLogger:
    """a-liner pipeline 日志管理器|aliner pipeline logger manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        level = getattr(logging, log_level.upper(), logging.INFO)
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))
        logging.basicConfig(level=level, format=log_format, datefmt=date_format, handlers=handlers)
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def extract_fasta_lengths(fasta_path: str, samtools_path: str, logger=None) -> Dict[str, int]:
    """
    用samtools faidx获取FASTA各序列长度|Get seq lengths via samtools faidx

    若无.fai则先建索引（samtools faidx）|build index if .fai absent
    """
    fai_path = f"{fasta_path}.fai"
    if not os.path.exists(fai_path):
        cmd = build_conda_command(samtools_path, ['faidx', fasta_path])
        if logger:
            logger.info(f"执行|Executing: samtools faidx 索引|samtools faidx index")
            logger.info(f"命令|Command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    lengths = {}
    with open(fai_path, 'r') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def check_dependencies(config, logger) -> bool:
    """检查依赖工具存在|Check dependency tools exist"""
    logger.info("检查依赖软件|Checking dependencies")
    for name, path in (('minimap2', config.minimap2_path), ('samtools', config.samtools_path)):
        if not os.path.exists(path):
            raise RuntimeError(f"{name} 不存在|{name} not found: {path}")
    logger.info("依赖检查通过|Dependencies OK")
    return True


def get_tool_version(tool_path: str, args: Optional[List[str]] = None) -> str:
    """获取工具版本|Get tool version (via build_conda_command)"""
    args = args or ['--version']
    try:
        cmd = build_conda_command(tool_path, args)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return (result.stdout.strip() or result.stderr.strip()) or 'unknown'
    except Exception:
        return 'unknown'


def get_aliner_version(aliner_env: str) -> str:
    """获取a-liner版本（固定环境）|Get a-liner version (fixed env)"""
    try:
        cmd = conda_env_run_prefix(aliner_env).split() + ['a-liner', '--version']
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return (result.stdout.strip() or result.stderr.strip()) or 'unknown'
    except Exception:
        return 'unknown'


def format_number(num: int) -> str:
    """格式化数字（M/K单位）|Format number (M/K units)"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    if num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
