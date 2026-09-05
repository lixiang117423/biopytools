"""fastANI工具函数模块|fastANI Utilities Module

日志管理器(三文件分离) + 基因组输入收集(目录/列表/单fasta)
|Logger manager (3-file separation) + genome input collection (dir/list/single fasta)
"""

import logging
import os
import sys
from typing import List

from biopytools.common.paths import expand_path

from .config import FASTA_GZ_SUFFIXES, FASTA_SUFFIXES


class FastaniLogger:
    """fastANI日志管理器|fastANI Logger Manager

    stdout(<=INFO) + stderr(>=WARNING) + 三个日志文件:
    fastani.log(全量) fastani_out.log(<=INFO) fastani_err.log(>=WARNING)
    |stdout (<=INFO) + stderr (>=WARNING) + three log files
    """

    def __init__(self, logs_dir: str, log_level: str = 'INFO'):
        self.log_file = os.path.join(logs_dir, 'fastani.log')
        self.out_log_file = os.path.join(logs_dir, 'fastani_out.log')
        self.err_log_file = os.path.join(logs_dir, 'fastani_err.log')
        os.makedirs(logs_dir, exist_ok=True)
        self.logger = self._setup_logging(log_level)

    def _setup_logging(self, log_level: str) -> logging.Logger:
        """设置日志|Setup logging(named logger,不污染root|no root pollution)"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger('biopytools.fastani')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        # stdout: <=INFO(受--log-level控制)|stdout: <=INFO (respects --log-level)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr: >=WARNING|stderr: >=WARNING
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 三文件|Three files
        specs = [(self.log_file, None),
                 (self.out_log_file, lambda r: r.levelno <= logging.INFO),
                 (self.err_log_file, lambda r: r.levelno >= logging.WARNING)]
        for path, level_filter in specs:
            handler = logging.FileHandler(path, encoding='utf-8')
            handler.setLevel(logging.DEBUG)
            if level_filter:
                handler.addFilter(level_filter)
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


def is_fasta_file(path: str) -> bool:
    """判断是否基因组FASTA(含gz)|Whether genome FASTA (incl. gz)"""
    name = os.path.basename(path).lower()
    return name.endswith(FASTA_SUFFIXES) or name.endswith(FASTA_GZ_SUFFIXES)


def genome_name(path: str) -> str:
    """基因组名=basename剥gz/fa/fna/fasta全部层|Genome name = basename, all suffix layers stripped"""
    name = os.path.basename(path)
    while True:
        root, ext = os.path.splitext(name)
        if ext.lower() in ('.gz',) + FASTA_SUFFIXES:
            name = root
        else:
            break
    return name


def collect_genome_files(path: str) -> List[str]:
    """收集基因组文件(三形态)|Collect genome files (dir/list/single fasta)

    目录→glob fasta;列表文件(非fasta后缀)→逐行读路径;fasta文件→自身
    |dir→glob fasta; list file (non-fasta suffix)→paths per line; fasta file→itself
    """
    errors = []
    genomes = []

    if os.path.isdir(path):
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            if os.path.isfile(full) and is_fasta_file(name):
                if os.path.getsize(full) == 0:
                    errors.append(f"基因组文件为空|Empty genome file: {full}")
                else:
                    genomes.append(os.path.abspath(full))
        if not genomes:
            errors.append(f"目录中未找到基因组文件|No genome files (.fa/.fna/.fasta[.gz]) "
                          f"in directory: {path}")
    elif os.path.isfile(path):
        if is_fasta_file(path):
            if os.path.getsize(path) == 0:
                errors.append(f"基因组文件为空|Empty genome file: {path}")
            else:
                genomes.append(os.path.abspath(path))
        else:
            # 列表文件|List file
            with open(path, encoding='utf-8') as fh:
                for line_no, raw in enumerate(fh, 1):
                    entry = raw.strip()
                    if not entry or entry.startswith('#'):
                        continue
                    entry = expand_path(entry)
                    if not os.path.isfile(entry):
                        errors.append(f"列表第{line_no}行文件不存在|List line {line_no} "
                                      f"not found: {entry}")
                    elif not is_fasta_file(entry):
                        errors.append(f"列表第{line_no}行非FASTA|List line {line_no} "
                                      f"not FASTA: {entry}")
                    elif os.path.getsize(entry) == 0:
                        errors.append(f"列表第{line_no}行: 基因组文件为空|Empty genome "
                                      f"file: {entry}")
                    else:
                        genomes.append(os.path.abspath(entry))
            if not genomes and not errors:
                errors.append(f"列表文件为空|Empty list file: {path}")
    else:
        errors.append(f"输入路径不存在|Input path not found: {path}")

    if errors:
        raise ValueError('\n'.join(errors))
    return genomes
