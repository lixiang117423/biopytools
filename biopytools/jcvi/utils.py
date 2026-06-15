"""
JCVI共享工具函数|JCVI Shared Utility Functions

从allelic_genes/utils.py迁移并泛化而来
Migrated and generalized from allelic_genes/utils.py
"""

import os
import subprocess
import logging
import shutil
import re
import sys
from pathlib import Path
from typing import Optional, List


class JcviLogger:
    """JCVI日志管理器|JCVI Logger Manager"""

    def __init__(self, log_file, name: str = "JCVI", log_level: str = "INFO"):
        self.log_file = Path(log_file)
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []
        self.logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        return self.logger


# 向后兼容别名|Backward compatibility alias
AllelicGenesLogger = JcviLogger


def _get_conda_base() -> Optional[str]:
    """获取conda base路径|Get conda base path

    按优先级尝试:
    1. CONDA_EXE 环境变量 (最可靠)
    2. shutil.which('conda') (次优)
    3. conda info --base 子进程调用 (兜底)
    """
    # 方法1: CONDA_EXE 环境变量
    conda_exe = os.environ.get('CONDA_EXE', '')
    if conda_exe:
        conda_base = os.path.dirname(os.path.dirname(conda_exe))
        if os.path.isdir(conda_base):
            return conda_base

    # 方法2: shutil.which
    conda_exe = shutil.which('conda')
    if conda_exe:
        conda_base = os.path.dirname(os.path.dirname(conda_exe))
        if os.path.isdir(conda_base):
            return conda_base

    # 方法3: 子进程调用
    try:
        result = subprocess.run(
            ['conda', 'info', '--base'],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            return result.stdout.strip().strip("'\"")
    except Exception:
        pass
    return None


_SENTINEL = object()
_conda_base_cache = _SENTINEL


def _get_conda_base_cached() -> Optional[str]:
    """获取conda base路径(带缓存)|Get conda base path (cached)"""
    global _conda_base_cache
    if _conda_base_cache is _SENTINEL:
        _conda_base_cache = _get_conda_base()
    return _conda_base_cache


def _get_conda_env_python(conda_env: str) -> Optional[str]:
    """获取conda环境内的python绝对路径|Get absolute python path within conda env"""
    conda_base = _get_conda_base_cached()
    if conda_base:
        python_path = os.path.join(conda_base, 'envs', conda_env, 'bin', 'python')
        if os.path.isfile(python_path):
            return python_path
    return None


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在conda环境|Detect conda environment for command"""
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)
    return None


def build_conda_command(command: str, args: List[str],
                        conda_env: Optional[str] = None) -> List[str]:
    """构建conda run命令|Build conda run command"""
    env_name = conda_env or get_conda_env(command)
    if env_name:
        return ['conda', 'run', '-n', env_name, '--no-capture-output', command] + args
    return [command] + args


def build_jcvi_command(module: str, args: List[str], conda_env: str) -> List[str]:
    """
    构建JCVI python -m调用命令|Build JCVI python -m invocation command

    优先使用conda环境内python的绝对路径, 避免conda run因PATH中存在
    同名python(如~/.local/bin/python)而调用错误解释器
    """
    env_python = _get_conda_env_python(conda_env)
    if env_python:
        return [env_python, '-m', module] + args
    conda_base = _get_conda_base_cached()
    if conda_base:
        env_python_fallback = os.path.join(conda_base, 'envs', conda_env, 'bin', 'python')
        if os.path.isfile(env_python_fallback):
            return [env_python_fallback, '-m', module] + args
    return ['conda', 'run', '-n', conda_env, '--no-capture-output',
            'python', '-m', module] + args


def discover_samples(input_dir: str, gff_ext: str, fa_ext: str,
                     logger: logging.Logger) -> List[str]:
    """
    在输入目录中发现样本前缀|Discover sample prefixes in input directory

    通过查找.gff文件确定样本,并检查对应.fa文件是否存在
    """
    gff_files = sorted(Path(input_dir).glob(f"*{gff_ext}"))
    if not gff_files:
        gff_files = sorted(Path(input_dir).rglob(f"*{gff_ext}"))

    samples = []
    for gff_file in gff_files:
        prefix = str(gff_file)
        if not prefix.endswith(gff_ext):
            continue
        prefix = prefix[:-len(gff_ext)]

        fa_file = prefix + fa_ext
        if Path(fa_file).exists():
            samples.append(prefix)
            logger.info(f"  发现样本|Found sample: {Path(prefix).name}")
        else:
            logger.warning(
                f"  样本 {Path(prefix).name} 缺少 {fa_ext} 文件, 跳过|"
                f"Sample {Path(prefix).name} missing {fa_ext}, skipped"
            )

    return samples


def get_sample_name(prefix_path: str) -> str:
    """从完整路径提取样本名(保留完整文件名)|Extract sample name from full path (keep full filename)"""
    return Path(prefix_path).name


def get_jcvi_stem(name: str) -> str:
    """
    获取JCVI内部使用的stem名称|Get stem name as used internally by JCVI tools

    JCVI工具(如bed uniq, align last等)用Path.stem处理文件名生成输出,
    会去掉最后一个点号后缀(如ZH13.v2→ZH13), 本函数复现此行为
    """
    return Path(name).stem
