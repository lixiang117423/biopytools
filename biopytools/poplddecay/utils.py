"""
PopLDdecay工具函数模块|PopLDdecay Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Tuple


class PopLDdecayLogger:
    """PopLDdecay日志管理器|PopLDdecay Logger Manager"""

    def __init__(self, log_file: Path):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
        """
        self.log_file = log_file

        # 创建logger|Create logger
        self.logger = logging.getLogger('PopLDdecay')
        self.logger.setLevel(logging.INFO)

        # 清除现有处理器|Clear existing handlers
        self.logger.handlers = []

        # 文件处理器|File handler
        fh = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        fh.setLevel(logging.INFO)

        # 控制台处理器|Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)

        # 格式化器（带毫秒）|Formatter with milliseconds
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

    def get_logger(self):
        """获取logger对象|Get logger object

        Returns:
            logger: Logger对象|Logger object
        """
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中|Detect if command is in a conda environment"""
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
    """构建命令|Build command

    绝对路径直接调用，命令名通过conda run包装|Absolute paths called directly, command names wrapped with conda run
    """
    if os.path.isabs(command):
        return [command] + args

    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def check_vcf_index(vcf_file: Path, logger: Optional[logging.Logger] = None) -> bool:
    """检查VCF文件是否有索引|Check if VCF file has index

    Args:
        vcf_file: VCF文件路径|VCF file path
        logger: 日志器|Logger

    Returns:
        bool: 是否有索引|Whether index exists
    """
    # 检查.tbi索引|Check .tbi index
    if vcf_file.with_suffix('.vcf.gz.tbi').exists():
        return True

    # 检查.csi索引|Check .csi index
    if vcf_file.with_suffix('.csi').exists():
        return True

    # 检查.vcf.gz + .tbi组合|Check .vcf.gz + .tbi combination
    if vcf_file.suffix == '.gz' and vcf_file.with_suffix('.tbi').exists():
        return True

    return False


def build_vcf_index(vcf_file: Path, logger: Optional[logging.Logger] = None) -> bool:
    """构建VCF索引|Build VCF index

    Args:
        vcf_file: VCF文件路径|VCF file path
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Success
    """
    import subprocess

    try:
        if logger:
            logger.info(f"构建VCF索引|Building VCF index: {vcf_file}")

        # 使用bcftools构建索引|Use bcftools to build index
        cmd = ['bcftools', 'index', str(vcf_file)]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            if logger:
                logger.info(f"VCF索引构建成功|VCF index built successfully")
            return True
        else:
            if logger:
                logger.error(f"VCF索引构建失败|Failed to build VCF index: {result.stderr}")
            return False

    except Exception as e:
        if logger:
            logger.error(f"构建VCF索引时出错|Error building VCF index: {e}")
        return False


def validate_subpop_file(subpop_file: Path, logger: Optional[logging.Logger] = None) -> bool:
    """验证子群体文件|Validate subpopulation file

    Args:
        subpop_file: 子群体文件路径|Subpopulation file path
        logger: 日志器|Logger

    Returns:
        bool: 是否有效|Valid
    """
    if not subpop_file.exists():
        if logger:
            logger.error(f"子群体文件不存在|Subpopulation file does not exist: {subpop_file}")
        return False

    try:
        with open(subpop_file, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]

        if len(samples) == 0:
            if logger:
                logger.error(f"子群体文件为空|Subpopulation file is empty: {subpop_file}")
            return False

        if logger:
            logger.info(f"子群体文件包含{len(samples)}个样本|Subpopulation file contains {len(samples)} samples")

        return True

    except Exception as e:
        if logger:
            logger.error(f"读取子群体文件时出错|Error reading subpopulation file: {e}")
        return False


def parse_subpop_file(subpop_file: Path, logger: Optional[logging.Logger] = None) -> dict:
    """解析子群体文件，提取群体信息|Parse subpopulation file and extract population information

    Args:
        subpop_file: 子群体文件路径|Subpopulation file path
        logger: 日志器|Logger

    Returns:
        dict: 群体字典 {群体名: [样本列表]}|Population dictionary {population_name: [sample_list]}
    """
    from collections import defaultdict

    populations = defaultdict(list)

    try:
        with open(subpop_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # 分割行（支持tab或空格分隔）|Split line (tab or space separated)
                parts = line.split('\t') if '\t' in line else line.split()

                if len(parts) >= 2:
                    sample_name = parts[0]
                    population_name = parts[1]
                    populations[population_name].append(sample_name)

        if logger:
            logger.info(f"识别到{len(populations)}个子群体|Found {len(populations)} subpopulations:")
            for pop_name, samples in sorted(populations.items()):
                logger.info(f"  {pop_name}: {len(samples)}个样本|{len(samples)} samples")

        return dict(populations)

    except Exception as e:
        if logger:
            logger.error(f"解析子群体文件时出错|Error parsing subpopulation file: {e}")
        return {}


def check_perl_modules(logger: Optional[logging.Logger] = None) -> bool:
    """检查Perl模块是否可用|Check if Perl modules are available

    Args:
        logger: 日志器|Logger

    Returns:
        bool: 是否可用|Available
    """
    import subprocess

    required_modules = ['Data::Dumper', 'Getopt::Long']

    for module in required_modules:
        try:
            cmd = ['perl', '-M' + module, '-e', '1']
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                if logger:
                    logger.warning(f"Perl模块{module}不可用|Perl module {module} is not available")
                return False

        except Exception as e:
            if logger:
                logger.warning(f"检查Perl模块{module}时出错|Error checking Perl module {module}: {e}")
            return False

    return True
