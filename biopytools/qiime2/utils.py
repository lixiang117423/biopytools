"""
QIIME2工具函数模块|QIIME2 Utility Functions Module

包含日志管理、conda命令包装、命令执行器、双端FASTQ配对检测、引物反向互补等
|Includes logging, conda command wrapping, command runner, paired FASTQ detection, primer reverse complement
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# IUPAC互补碱基映射(DNA: A↔T, C↔G)|IUPAC complementary base mapping (DNA: A↔T, C↔G)
IUPAC_COMPLEMENT = str.maketrans(
    'ACGTURYSWKMBDHVNacgturyswkmbdhvn',
    'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn'
)


class Qiime2Logger:
    """QIIME2日志管理器|QIIME2 Logger Manager

    三handler分离(§2.3): INFO→stdout, WARNING+→stderr, DEBUG+→文件
    |Three-handler separation: INFO→stdout, WARNING+→stderr, DEBUG+→file
    """

    def __init__(self, log_file: str):
        self.log_file = Path(log_file)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 创建日志目录|Create log directory
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 删除旧日志|Delete old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 标准日志格式|Standard log format (§2.1)
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 配置日志器|Configure logger
        logger = logging.getLogger('qiime2')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False  # 避免重复|Avoid duplicates

        logger.addHandler(file_handler)
        logger.addHandler(stdout_handler)
        logger.addHandler(stderr_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中,返回环境名称|Detect conda environment name from command path

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    # 方法1: 从命令路径检测|Method 1: detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 命令本身是绝对路径时也检测|Also detect when command is an absolute path
    if command.startswith('/'):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: search all conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            cmd_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', cmd_name)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来调用conda环境中的软件|Build conda run command

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表(配合subprocess.run(shell=False))|Complete command list (for shell=False)

    Note:
        必须传递完整路径(§13.6),否则无法提取环境名
        |Must pass full path, otherwise env name cannot be extracted
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # --no-capture-output避免conda缓冲输出导致内存问题(§13.2.0)
        # |--no-capture-output avoids conda buffering output causing memory issues
        full_cmd = [
            'conda', 'run', '-n', conda_env,
            '--no-capture-output', command
        ] + args
    else:
        # 非conda环境,直接调用|Non-conda environment, direct call
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner

    始终使用shell=False + 列表(§13),执行前记录完整命令到INFO级别(§2.2.1)
    |Always shell=False + list, log full command at INFO before execution
    """

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

        # 记录完整命令到INFO级别(§2.2.1)|Log full command at INFO level
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                self.logger.error(f"错误代码|Error code: {result.returncode}")
                if result.stderr:
                    # 记录最后3000字符|Log last 3000 chars
                    stderr_tail = result.stderr[-3000:] if len(result.stderr) > 3000 else result.stderr
                    self.logger.error(f"错误信息|Error message: {stderr_tail}")

                    # 保存完整stderr到文件|Save full stderr to file
                    stderr_file = os.path.join(
                        os.path.dirname(self.working_dir) if self.working_dir else '.',
                        "99_logs",
                        f"stderr_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
                    )
                    try:
                        os.makedirs(os.path.dirname(stderr_file), exist_ok=True)
                        with open(stderr_file, 'w') as f:
                            f.write(result.stderr)
                        self.logger.error(f"完整错误日志已保存|Full error log saved: {stderr_file}")
                    except Exception:
                        pass
                return False, result.stdout, result.stderr

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[-1000:]}")

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True, result.stdout, result.stderr

        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, '', str(e)
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {description}")
            self.logger.error(f"异常信息|Exception: {str(e)}")
            return False, '', str(e)


class PairFinder:
    """双端FASTQ文件配对检测器|Paired FASTQ File Finder

    按R1/R2后缀配对样品|Pair samples by R1/R2 suffixes
    """

    def __init__(self, input_dir: str, r1_suffix: str, r2_suffix: str):
        self.input_dir = input_dir
        self.r1_suffix = r1_suffix
        self.r2_suffix = r2_suffix

    def find_pairs(self) -> Dict[str, Tuple[str, str]]:
        """
        检测配对的FASTQ文件|Detect paired FASTQ files

        Returns:
            {sample_name: (r1_path, r2_path)}, 按样品名排序|sorted by sample name
        """
        r1_files: Dict[str, str] = {}
        r2_files: Dict[str, str] = {}

        for fname in sorted(os.listdir(self.input_dir)):
            fpath = os.path.join(self.input_dir, fname)
            if not os.path.isfile(fpath):
                continue

            if fname.endswith(self.r1_suffix):
                sample_name = fname[:-len(self.r1_suffix)]
                r1_files[sample_name] = fpath
            elif fname.endswith(self.r2_suffix):
                sample_name = fname[:-len(self.r2_suffix)]
                r2_files[sample_name] = fpath

        pairs: Dict[str, Tuple[str, str]] = {}
        for sample_name in sorted(r1_files.keys()):
            if sample_name in r2_files:
                pairs[sample_name] = (r1_files[sample_name], r2_files[sample_name])

        return pairs


def generate_manifest(pairs: Dict[str, Tuple[str, str]], manifest_path: str) -> None:
    """
    生成QIIME2 V2 manifest文件|Generate QIIME2 V2 manifest file

    V2格式: 带表头TSV,列为sample-id/forward-absolute-filepath/reverse-absolute-filepath
    |V2 format: TSV with header, columns sample-id/forward/reverse-absolute-filepath
    文件路径必须为绝对路径|Filepaths must be absolute

    Args:
        pairs: {sample_name: (r1_path, r2_path)}
        manifest_path: manifest输出路径|manifest output path
    """
    with open(manifest_path, 'w') as f:
        f.write('sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')
        for sample_name in sorted(pairs.keys()):
            r1, r2 = pairs[sample_name]
            f.write(f"{sample_name}\t{os.path.abspath(r1)}\t{os.path.abspath(r2)}\n")


def reverse_complement(seq: str) -> str:
    """
    计算IUPAC序列的反向互补|Compute reverse complement of an IUPAC sequence

    用于cutadapt的--p-front-r参数(需传反向引物的反向互补序列)
    |Used for cutadapt --p-front-r (needs reverse complement of reverse primer)

    Args:
        seq: 核苷酸序列(可含IUPAC简并码)|nucleotide sequence (may include IUPAC degenerate codes)

    Returns:
        反向互补序列|reverse complement sequence

    Examples:
        >>> reverse_complement('GACTACHVGGGTATCTAATCC')
        'GGATTAGATACCCBDGTAGTC'
    """
    return seq.translate(IUPAC_COMPLEMENT)[::-1]


def format_number(num: int) -> str:
    """格式化大数字(§5.3)|Format large numbers with M/K units"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def generate_software_versions_yml(
    output_dir: str,
    config,
    classifier_source: str,
    sampling_depth: int,
    start_time: datetime
) -> None:
    """
    生成software_versions.yml文件|Generate software_versions.yml file

    Args:
        output_dir: 输出目录|Output directory
        config: Qiime2Config实例|Qiime2Config instance
        classifier_source: 分类器来源描述|Classifier source description
        sampling_depth: 使用的抽样深度|Sampling depth used
        start_time: 开始时间|Start time
    """
    try:
        import yaml
    except ImportError:
        logging.getLogger('qiime2').warning(
            "PyYAML未安装,跳过生成software_versions.yml|PyYAML not installed, skipping software_versions.yml"
        )
        return

    end_time = datetime.now()
    runtime_seconds = int((end_time - start_time).total_seconds())

    # 获取qiime版本|Get qiime version
    qiime_version = 'unknown'
    try:
        cmd = build_conda_command(config.qiime_bin, ['--version'])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            qiime_version = result.stdout.strip().splitlines()[0] if result.stdout.strip() else 'unknown'
    except Exception:
        pass

    info = {
        'pipeline': {
            'name': 'biopytools qiime2',
            'version': '1.0.0'
        },
        'tools': {
            'qiime2': {
                'version': qiime_version,
                'path': config.qiime_bin,
                'conda_env': get_conda_env(config.qiime_bin)
            }
        },
        'classifier': {
            'amplicon': config.amplicon,
            'source': classifier_source,
            'database_dir': config.database_dir
        },
        'parameters': {
            'input_dir': config.input_dir,
            'output_dir': config.output_dir,
            'amplicon': config.amplicon,
            'method': config.method,
            'fwd_primer': config.fwd_primer,
            'rev_primer': config.rev_primer,
            'trunc_len_f': config.trunc_len_f,
            'trunc_len_r': config.trunc_len_r,
            'trim_left_f': config.trim_left_f,
            'trim_left_r': config.trim_left_r,
            'sampling_depth': sampling_depth,
            'perc_identity': config.perc_identity,
            'confidence': config.confidence,
            'threads': config.threads,
            'skip_cutadapt': config.skip_cutadapt,
            'skip_phylogeny': config.skip_phylogeny,
            'r1_suffix': config.r1_suffix,
            'r2_suffix': config.r2_suffix
        },
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': runtime_seconds
        }
    }

    output_file = Path(output_dir) / '00_pipeline_info' / 'software_versions.yml'
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
