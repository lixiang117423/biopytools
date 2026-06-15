"""泛基因组Block构建工具函数模块|Pan-Genome Block Construction Utility Functions Module"""

import logging
import sys
import subprocess
import os
import shutil
import re
import tempfile
from typing import Optional, List, Tuple, Dict


class PanBlocksLogger:
    """Pan-Blocks日志管理器|Pan-Blocks Logger Manager"""

    def __init__(self, output_dir: str, log_file: str = 'pan_blocks.log', log_level: str = "INFO"):
        self.log_file = log_file
        self.log_level = log_level
        self.output_dir = output_dir

        log_dir = os.path.join(output_dir, '99_logs')
        os.makedirs(log_dir, exist_ok=True)
        self.log_path = os.path.join(log_dir, log_file)

        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        root_logger = logging.getLogger()
        root_logger.handlers.clear()

        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        file_handler = logging.FileHandler(self.log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        root_logger.setLevel(level)
        root_logger.addHandler(stdout_handler)
        root_logger.addHandler(stderr_handler)
        root_logger.addHandler(file_handler)
        root_logger.propagate = False

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name"""
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
    """构建conda run命令|Build conda run command"""
    conda_env = get_conda_env(command)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def run(self, cmd: List[str], description: str = "",
            cwd: Optional[str] = None) -> Tuple[bool, str, str]:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, shell=False, capture_output=True, text=True,
                check=False, cwd=cwd
            )
            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {description}")
                self.logger.error(f"错误输出|Error output: {result.stderr}")
                return False, result.stdout, result.stderr
            if description:
                self.logger.info(f"完成|Completed: {description}")
            return True, result.stdout, result.stderr
        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {description}: {e}")
            return False, "", str(e)

    def run_to_file(self, cmd: List[str], output_file: str, description: str = "",
                    cwd: Optional[str] = None) -> bool:
        """执行命令，输出到文件|Execute command with output to file"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)} > {output_file}")

        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd, stdout=f, stderr=subprocess.PIPE,
                    text=False, check=False, cwd=cwd
                )
            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Error output: {result.stderr.decode('utf-8', errors='ignore')}")
                return False
            if description:
                self.logger.info(f"完成|Completed: {description}")
            return True
        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {description}: {e}")
            return False

    def run_to_file_with_input(self, cmd: List[str], input_file: str, output_file: str,
                               description: str = "", cwd: Optional[str] = None) -> bool:
        """执行命令，从文件读取stdin，输出到文件|Execute command with file input and file output"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)} < {input_file} > {output_file}")

        try:
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                result = subprocess.run(
                    cmd, stdin=fin, stdout=fout, stderr=subprocess.PIPE,
                    text=False, check=False, cwd=cwd
                )
            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Error output: {result.stderr.decode('utf-8', errors='ignore')}")
                return False
            if description:
                self.logger.info(f"完成|Completed: {description}")
            return True
        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {description}: {e}")
            return False


def parse_fasta_lengths(fasta_path: str) -> Dict[str, int]:
    """解析FASTA文件获取序列长度|Parse FASTA file to get sequence lengths"""
    lengths = {}
    current_name = None
    current_length = 0

    if fasta_path.endswith('.gz'):
        import gzip
        opener = gzip.open
    else:
        opener = open

    with opener(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if current_name is not None:
                    lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line.strip())
        if current_name is not None:
            lengths[current_name] = current_length
    return lengths


def parse_chrlen_file(path: str) -> Dict[str, Dict[str, int]]:
    """解析chrlen文件|Parse chromosome length file

    格式|Format: genome_name<TAB>chr_name<TAB>chr_length (每行一个基因组-染色体对)
    """
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            genome, chr_name, chr_len = parts[0], parts[1], int(parts[2])
            if genome not in data:
                data[genome] = {}
            data[genome][chr_name] = chr_len
    return data


def parse_coords_file(path: str) -> List[dict]:
    """解析show-coords -TrHcl输出文件|Parse show-coords -TrHcl output file

    MUMmer 4.x show-coords -TrHcl 输出列 (tab分隔|tab-separated):
    [0]REF_START [1]REF_END [2]QRY_START [3]QRY_END [4]%IDY
    [5]REF_LEN [6]QRY_LEN [7]%COV_REF [8]%COV_QRY
    [9]TAGS [10]FRM [11]REF_CHR [12]QRY_CHR
    """
    alignments = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('[R'):
                continue
            parts = line.split('\t')
            if len(parts) < 13:
                continue
            try:
                aln = {
                    'ref_start': int(parts[0]),
                    'ref_end': int(parts[1]),
                    'qry_start': int(parts[2]),
                    'qry_end': int(parts[3]),
                    'identity': float(parts[4]),
                    'ref_len': int(parts[5]),
                    'qry_len': int(parts[6]),
                    'ref_chr': parts[11],
                    'qry_chr': parts[12],
                }
                alignments.append(aln)
            except (ValueError, IndexError):
                continue
    return alignments


def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    tools = [
        ('nucmer', config.nucmer_path),
        ('delta-filter', config.delta_filter_path),
        ('show-coords', config.show_coords_path),
        ('bedtools', config.bedtools_path),
    ]

    all_ok = True
    for name, path in tools:
        if os.path.isfile(path):
            cmd = build_conda_command(path, ['--version'] if name != 'bedtools' else ['--version'])
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    version = result.stdout.strip().split('\n')[0]
                    logger.info(f"{name} 可用|{name} available: {version}")
                else:
                    logger.warning(f"{name} 可执行但版本检查失败|{name} executable but version check failed")
            except Exception:
                logger.warning(f"{name} 可执行但版本检查异常|{name} executable but version check error")
        else:
            logger.error(f"{name} 不存在|{name} not found: {path}")
            all_ok = False

    if all_ok:
        logger.info("依赖检查完成|Dependency check completed")
    else:
        logger.error("部分依赖缺失|Some dependencies missing")

    return all_ok


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def write_bed_file(regions: List[Tuple[str, int, int]], output_path: str) -> None:
    """写入BED文件|Write BED file (0-based half-open)"""
    with open(output_path, 'w') as f:
        for chr_name, start, end in sorted(regions, key=lambda x: (x[0], x[1], x[2])):
            f.write(f"{chr_name}\t{start}\t{end}\n")


def read_bed_file(path: str) -> List[Tuple[str, int, int]]:
    """读取BED文件|Read BED file"""
    regions = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions
