"""gene_table 工具函数|gene_table utilities (logging, conda wrap, FASTA io, seq ops)"""

import gzip
import logging
import os
import re
import shutil
import sys
from typing import Dict, List, Optional, Tuple

_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


class GeneTableLogger:
    """三分日志:stdout=INFO / stderr=WARNING+ / file=DEBUG+|3-way split logger"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = 'INFO',
                 verbose: bool = False):
        self.log_file = log_file
        self.log_level = logging.DEBUG if verbose else getattr(
            logging, str(log_level).upper(), logging.INFO)

    def get_logger(self) -> logging.Logger:
        """构建并返回配置好的 logger|Build and return a configured logger"""
        logger = logging.getLogger('gene_table')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False  # 避免重复输出|avoid duplicate output

        fmt = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')

        stdout_h = logging.StreamHandler(sys.stdout)
        stdout_h.setLevel(self.log_level)
        stdout_h.setFormatter(fmt)
        logger.addHandler(stdout_h)

        stderr_h = logging.StreamHandler(sys.stderr)
        stderr_h.setLevel(logging.WARNING)
        stderr_h.setFormatter(fmt)
        logger.addHandler(stderr_h)

        if self.log_file:
            file_h = logging.FileHandler(self.log_file)
            file_h.setLevel(logging.DEBUG)
            file_h.setFormatter(fmt)
            logger.addHandler(file_h)
        return logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所属 conda 环境(从完整路径 /envs/<name>/ 提取)|Detect conda env from full path

    策略:完整路径下只看 /envs/<name>/ 命中;只有裸命令名(无路径分隔符)才扫描所有环境兜底。
    传完整路径但不在 /envs/ 下 → 返回 None(直接调用),避免把 ~/.local/bin 的独立二进制
    误包进某个 conda 环境。|A full path not under /envs/ returns None (direct call), so a
    standalone ~/.local/bin binary is never mis-wrapped into some conda env.
    """
    cmd_path = shutil.which(command) or command
    match = re.search(r'/envs/([^/]+)/', cmd_path)
    if match:
        return match.group(1)
    # 仅对裸命令名(无路径分隔符)扫描所有 conda 环境兜底|Scan envs only for bare command names
    if os.sep not in command and '/' not in command:
        conda_exe = os.environ.get('CONDA_EXE')
        if conda_exe:
            envs_dir = os.path.join(os.path.dirname(os.path.dirname(conda_exe)), 'envs')
            if os.path.isdir(envs_dir):
                for env_name in os.listdir(envs_dir):
                    if os.path.exists(os.path.join(envs_dir, env_name, 'bin', command)):
                        return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """包装 conda 环境内命令;非 conda 则直接调用|Wrap conda-env command; direct call otherwise

    必须传完整路径,禁止 basename|Must pass the full path, never the basename.
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def reverse_complement(seq: str) -> str:
    """反向互补(ACGT/N,其余碱基→N)|Reverse complement"""
    return ''.join(_COMPLEMENT.get(b.upper(), 'N') for b in reversed(seq))


def extract_sequence_region(genome_seq: str, start: int, end: int, strand: str) -> str:
    """按 GFF 1-based 坐标切片,负链反向互补|Slice by GFF 1-based coords, revcomp on minus strand"""
    seq = genome_seq[start - 1:end]
    return reverse_complement(seq) if strand == '-' else seq


def format_number(num: int) -> str:
    """大数字(≥1M)用 M 单位|Format big numbers with M unit"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


def _open_text(path: str):
    """透明打开 .gz 文本|Transparently open .gz text"""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def read_fasta_to_dict(path: str) -> Dict[str, str]:
    """读 FASTA→{id: seq},id 取 header 首个空白 token|Read FASTA into dict, id = first header token"""
    seqs: Dict[str, str] = {}
    cur_id, chunks = None, []
    with _open_text(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(chunks)
                cur_id = line[1:].split()[0] if len(line) > 1 else ''
                chunks = []
            elif line:
                chunks.append(line.strip())
    if cur_id is not None:
        seqs[cur_id] = ''.join(chunks)
    return seqs


def write_fasta(pairs: List[Tuple[str, str]], path: str, line_width: int = 60) -> None:
    """写 FASTA (id, seq) 列表,按行宽折行|Write FASTA pairs, wrap at line_width"""
    with open(path, 'w') as f:
        for seq_id, seq in pairs:
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + '\n')
