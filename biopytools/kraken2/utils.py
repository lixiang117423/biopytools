"""
Kraken2工具函数模块|Kraken2 Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class Kraken2Logger:
    """Kraken2日志管理器|Kraken2 Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        self.logger = logging.getLogger('kraken2')
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler|File handler
        if self.log_file:
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect conda environment name from command path

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
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
    """
    构建conda run命令|Build conda run command

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        full_cmd = [
            'conda', 'run', '-n', conda_env,
            '--no-capture-output', command
        ] + args
    else:
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Optional[str] = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: List[str], description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description

        Returns:
            (success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        cmd_str = ' '.join(cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

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
            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False, e.stdout, e.stderr

        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, '', str(e)


class PairFinder:
    """配对FASTQ文件检测器|Paired FASTQ File Finder"""

    def __init__(self, input_dir: str, r1_suffix: str, r2_suffix: str):
        self.input_dir = input_dir
        self.r1_suffix = r1_suffix
        self.r2_suffix = r2_suffix

    def find_pairs(self) -> Dict[str, Tuple[str, str]]:
        """
        检测配对的FASTQ文件|Detect paired FASTQ files

        Returns:
            {sample_name: (r1_path, r2_path)}
        """
        r1_files = {}
        r2_files = {}

        for fname in sorted(os.listdir(self.input_dir)):
            if fname.endswith(self.r1_suffix):
                sample_name = fname[:-len(self.r1_suffix)]
                r1_files[sample_name] = os.path.join(self.input_dir, fname)
            elif fname.endswith(self.r2_suffix):
                sample_name = fname[:-len(self.r2_suffix)]
                r2_files[sample_name] = os.path.join(self.input_dir, fname)

        pairs = {}
        r1_only = set(r1_files.keys()) - set(r2_files.keys())
        r2_only = set(r2_files.keys()) - set(r1_files.keys())

        for sample_name in sorted(r1_files.keys()):
            if sample_name in r2_files:
                pairs[sample_name] = (r1_files[sample_name], r2_files[sample_name])

        return pairs


def format_number(num: int) -> str:
    """格式化大数字|Format large numbers"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def parse_kraken2_report(report_file: str) -> dict:
    """
    解析Kraken2报告文件|Parse Kraken2 report file

    Kraken2 report格式|Kraken2 report format:
    percentage  clade_reads  taxid  rank_code  name

    Args:
        report_file: Kraken2报告文件路径|Kraken2 report file path

    Returns:
        解析结果字典|Parsed result dictionary
    """
    result = {
        'total_reads': 0,
        'classified_reads': 0,
        'unclassified_reads': 0,
        'classified_pct': 0.0,
        'unclassified_pct': 0.0,
        'species': []
    }

    if not os.path.exists(report_file):
        return result

    with open(report_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            pct = float(fields[0])
            reads = int(fields[1])
            taxid = fields[2]
            rank = fields[3]
            name = fields[4].strip()

            if taxid == '0':
                result['unclassified_reads'] = reads
                result['unclassified_pct'] = pct
            elif rank == 'U':
                result['classified_reads'] = reads
                result['classified_pct'] = pct
                result['total_reads'] = int(reads / pct * 100) if pct > 0 else reads
            elif rank == 'S':
                result['species'].append({
                    'name': name,
                    'reads': reads,
                    'percentage': pct,
                    'taxid': taxid
                })

    result['species'].sort(key=lambda x: x['reads'], reverse=True)

    return result
