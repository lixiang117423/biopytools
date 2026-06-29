"""a-liner pipeline 工具模块|aliner pipeline utility module"""

import os
import re
import sys
import shutil
import logging
import subprocess
from typing import List, Optional, Tuple, Dict


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令所在conda环境|Detect conda env of command

    策略|strategy:
        1. 从命令完整路径的 /envs/<name>/ 检测|from /envs/<name>/ in full path
        2. 从 CONDA_EXE 搜索所有环境兜底|fallback: search all envs via CONDA_EXE
    """
    # 方法1: 完整路径直接正则|method 1: regex on full path
    match = re.search(r'/envs/([^/]+)/bin/', command)
    if match:
        return match.group(1)
    # 方法1b: which 解析命令名|resolve via which
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)
    # 方法2: 搜索所有conda环境|method 2: search all envs
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, 'bin', command)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令|Build conda run command

    注意：command 必须是完整路径（含 /envs/<env>/bin/），禁止用 os.path.basename 提取命令名（§13.6）
    Note: command must be a full path; never reduce via os.path.basename
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


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
        cmd = ['conda', 'run', '-n', aliner_env, '--no-capture-output', 'a-liner', '--version']
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
