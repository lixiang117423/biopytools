"""
NLR-Annotator工具函数模块|NLR-Annotator Utility Functions Module
"""

import glob
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path
    """
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

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
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class NLRLogger:
    """NLR-Annotator日志管理器|NLR-Annotator Logger Manager"""

    def __init__(self, output_dir: str = "./output", log_name: str = "nlr_annotator.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / "99_logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logger = logging.getLogger("nlr_annotator")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        class _MaxLevelFilter(logging.Filter):
            def __init__(self, max_level):
                super().__init__()
                self.max_level = max_level

            def filter(self, record):
                return record.levelno <= self.max_level

        # stdout handler - 仅INFO|stdout handler - INFO only
        # INFO→stdout→.out, WARNING+→stderr→.err
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(_MaxLevelFilter(logging.INFO))
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


def clean_output(output_file: str, logger: logging.Logger):
    """
    清洗NLR-Annotator输出：加表头、去重并排序motif列表|Clean NLR-Annotator output: add header, deduplicate and sort motifs
    """
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        logger.warning(f"输出文件为空，跳过清洗|Output file is empty, skipping clean")
        return

    with open(output_file, 'r') as f:
        lines = f.read().strip().split('\n')

    cleaned = ['gene_id\tnlr_id\ttype\tstart\tend\tstrand\tmotifs']
    for line in lines:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 7:
            cleaned.append(line.strip())
            continue
        # motif列：去重、排序、去掉motif_前缀|Deduplicate, sort, remove motif_ prefix
        motifs = sorted(set(m.replace('motif_', '') for m in fields[6].split(',')))
        fields[6] = ','.join(motifs)
        cleaned.append('\t'.join(fields))

    with open(output_file, 'w') as f:
        f.write('\n'.join(cleaned) + '\n')

    logger.info(f"输出已清洗(含表头)|Output cleaned (with header): {output_file}")


def extract_sample_name(filename: str, sample_suffix: str) -> str:
    """
    从文件名提取样本名|Extract sample name from filename

    Args:
        filename: 文件名（含扩展名）|Filename with extension
        sample_suffix: 匹配后缀，如"*.cds.fa"|Match suffix, e.g. "*.cds.fa"
    """
    escaped = re.escape(sample_suffix).replace(r'\*', '(.*)')
    match = re.match(escaped, filename)
    if match and match.group(1):
        return match.group(1)
    return Path(filename).stem


def collect_input_files(input_path: str, sample_suffix: str, logger: logging.Logger) -> List[Tuple[str, str]]:
    """
    收集输入文件列表|Collect input file list

    支持单文件或目录批处理|Supports single file or directory batch processing

    Args:
        input_path: 输入文件或目录路径|Input file or directory path
        sample_suffix: 目录模式下文件匹配后缀|File match suffix in directory mode
        logger: 日志器|Logger

    Returns:
        [(file_path, sample_name), ...]
    """
    path = Path(input_path)

    if path.is_file():
        sample_name = path.stem
        logger.info(f"发现单个输入文件|Found single input file: {path.name}")
        return [(str(path), sample_name)]

    if path.is_dir():
        pattern = str(path / sample_suffix)
        files = sorted(glob.glob(pattern))
        if not files:
            raise ValueError(f"目录中未找到匹配'{sample_suffix}'的文件|No files matching '{sample_suffix}' found in directory")

        results = []
        for f in files:
            name = extract_sample_name(Path(f).name, sample_suffix)
            results.append((f, name))

        logger.info(f"发现批处理文件|Found batch files: {len(results)} 个文件|files")
        return results

    raise ValueError(f"输入路径无效|Invalid input path: {input_path}")


def generate_summary(sample_results: List[Tuple[str, str]], output_path: Path, logger: logging.Logger):
    """
    生成多样本汇总文件|Generate multi-sample summary file

    Args:
        sample_results: [(sample_name, tsv_path), ...]
        output_path: 输出目录|Output directory
        logger: 日志器|Logger
    """
    summary_file = output_path / "nlr_annotator_summary.tsv"

    all_rows = ['gene_id\tnlr_id\tsample\ttype\tstart\tend\tstrand\tmotifs']
    for sample_name, tsv_path in sample_results:
        if not os.path.exists(tsv_path):
            continue
        with open(tsv_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('gene_id'):
                    continue
                fields = line.split('\t')
                if len(fields) >= 6:
                    row = f"{fields[0]}\t{fields[1]}\t{sample_name}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}"
                    if len(fields) >= 7:
                        row += f"\t{fields[6]}"
                    all_rows.append(row)

    with open(summary_file, 'w') as f:
        f.write('\n'.join(all_rows) + '\n')

    logger.info(f"汇总文件已生成|Summary file generated: {summary_file}")
    logger.info(f"汇总NLR记录数|Total NLR records in summary: {len(all_rows) - 1}")
