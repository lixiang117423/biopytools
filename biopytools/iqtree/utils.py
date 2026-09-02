"""
IQ-TREE分析工具函数模块|IQ-TREE Analysis Utility Functions Module
"""

import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in a conda environment

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    import re
    import os

    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建conda run命令来运行conda环境中的软件|Build conda run command for conda env software

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


# FASTA常见后缀|Common FASTA suffixes
_FASTA_SUFFIXES = {'.fasta', '.fa', '.fna', '.faa', '.fsa'}
# PHYLIP常见后缀|Common PHYLIP suffixes
_PHYLIP_SUFFIXES = {'.phy', '.phylip', '.ph'}
# 蛋白特异字符(不在DNA IUPAC字母表中)|Protein-only chars absent from DNA IUPAC alphabet
_PROTEIN_ONLY_CHARS = set('EFILPQZJ*')
# 不参与判定的字符(缺口/两可字符/数字)|Chars that don't vote (gaps/ambiguous/digits)
_NEUTRAL_CHARS = set('-?.~ 0123456789')


def sniff_sequence_type(alignment_file: str) -> Optional[str]:
    """嗅探比对文件的序列类型|Sniff sequence type of an alignment file

    Why: IQ-TREE 3.x自动检测对简并码/N富集(非ACGT字符全比对>=10%)的DNA比对
    误报 "Unknown sequence type" 后退出(vcf2tree模块2026-08-18在3.1.3实测,
    2026-09-02 psoja_365再次复现), 故DNA比对需显式-st。按字母表嗅探: 含任一
    蛋白特异字符即AA, 否则视为DNA(蛋白检测本身不受此bug影响, -st AA同样安全)。
    仅嗅探FASTA/PHYLIP; NEXUS等复杂格式返回None, 交回IQ-TREE自动检测(现状行为)。
    |Why: IQ-TREE 3.x auto-detection aborts with "Unknown sequence type" on
    ambiguity-rich DNA alignments (>=10% non-ACGT alignment-wide; reproduced on
    3.1.3), so DNA needs explicit -st. Sniff by alphabet: any protein-only char
    means AA, otherwise DNA. Only FASTA/PHYLIP are sniffed; NEXUS etc. return
    None and fall back to IQ-TREE auto-detection (previous behavior).

    Args:
        alignment_file: 比对文件路径|Alignment file path

    Returns:
        'DNA'或'AA'; 无法判断时None|'DNA' or 'AA'; None if undeterminable
    """
    path = Path(alignment_file)
    suffix = path.suffix.lower()
    is_fasta = suffix in _FASTA_SUFFIXES
    is_phylip = suffix in _PHYLIP_SUFFIXES
    if not (is_fasta or is_phylip):
        return None

    chars = set()
    first_line = True
    try:
        with open(path, 'r', errors='ignore') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                if is_fasta:
                    if line.startswith('>'):
                        continue
                    seq = line
                else:
                    # PHYLIP头行形如"366 147677"|PHYLIP header like "366 147677"
                    tokens = line.split()
                    if first_line and len(tokens) == 2 and all(t.isdigit() for t in tokens):
                        first_line = False
                        continue
                    first_line = False
                    # 名字与序列之间必有空白(顺序格式首块); 交错格式后续块无名字
                    # (无空白, 切不开)→整行即序列, 不会误剥序列前10列
                    # |Name and sequence are separated by whitespace (sequential
                    # first block); interleaved later blocks have no name and no
                    # whitespace, so the whole line is sequence
                    parts = line.split(None, 1)
                    seq = parts[1] if len(parts) == 2 else parts[0]

                chars.update(seq.upper())
    except OSError:
        return None

    chars -= _NEUTRAL_CHARS
    if not chars:
        return None
    if chars & _PROTEIN_ONLY_CHARS:
        return 'AA'
    return 'DNA'


class TreeLogger:
    """系统发育树分析日志管理器|Phylogenetic Tree Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "iqtree_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

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

    def run(self, cmd: list, description: str = "") -> bool:
        """执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

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

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False


def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        cmd = build_conda_command(config.iqtree_path, ['--version'])
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            logger.info(f"IQ-TREE 可用|IQ-TREE available: {version_line}")
            return True
        else:
            raise RuntimeError("IQ-TREE 未正确安装|IQ-TREE not properly installed")

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"缺少依赖软件|Missing dependency: IQ-TREE"
        logger.error(error_msg)
        logger.error(f"请确保IQ-TREE已安装并在PATH中|Please ensure IQ-TREE is installed and in PATH")
        raise RuntimeError(error_msg)
