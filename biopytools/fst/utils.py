"""
Fst计算工具函数模块|Fst Calculation Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, List


class FstLogger:
    """Fst计算日志管理器|Fst Calculation Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "fst_calculation.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / '99_logs' / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 清空旧日志|Clear old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler - INFO级别|Stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, return environment name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs
        # Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件
    Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    # 提取命令的基本名称（去掉路径）|Extract command basename (remove path)
    cmd_name = Path(command).name

    if conda_env:
        # 使用conda run调用，添加--no-capture-output避免内存问题
        # Use conda run with --no-capture-output to avoid memory issues
        # 直接使用命令名，conda会自动处理PATH|Use command name directly, conda handles PATH
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', cmd_name] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 记录完整命令（INFO级别）|Log complete command (INFO level)
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=None  # 不设置工作目录，使用当前目录
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")  # 只显示前500字符

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:500]}")
            return False


def parse_population_file(pop_file: str) -> dict:
    """
    解析群体文件，自动检测分隔符
    Parse population file, auto-detect separator

    Args:
        pop_file: 群体文件路径|Population file path

    Returns:
        样本到群体的映射字典|Sample to population mapping dictionary
    """
    sample_to_pop = {}

    with open(pop_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # 自动检测分隔符|Auto-detect separator
            if '\t' in line:
                parts = line.split('\t')
            elif ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()

            if len(parts) >= 2:
                sample_id = parts[0]
                population = parts[1]
                sample_to_pop[sample_id] = population

    return sample_to_pop


def create_plink_population_file(sample_to_pop: dict, output_file: str):
    """
    创建PLINK格式的群体文件
    Create PLINK format population file

    Args:
        sample_to_pop: 样本到群体的映射|Sample to population mapping
        output_file: 输出文件路径|Output file path
    """
    with open(output_file, 'w') as f:
        for sample_id, population in sample_to_pop.items():
            # PLINK --within格式: FID IID POP
            # PLINK --within format: FID IID POP
            f.write(f"{sample_id}\t{sample_id}\t{population}\n")


def analyze_populations(sample_to_pop: dict, min_samples: int,
                        exclude_pops: Optional[str] = None) -> dict:
    """
    分析群体信息，过滤小样本群体
    Analyze population information, filter small populations

    Args:
        sample_to_pop: 样本到群体的映射|Sample to population mapping
        min_samples: 最小样本数阈值|Minimum sample count threshold
        exclude_pops: 手动指定要排除的群体（逗号分隔）|Manually specify populations to exclude

    Returns:
        群体统计信息字典|Population statistics dictionary
    """
    from collections import Counter

    # 统计各群体样本数|Count samples per population
    pop_counts = Counter(sample_to_pop.values())

    # 解析手动排除的群体|Parse manually excluded populations
    exclude_list = []
    if exclude_pops:
        exclude_list = [p.strip() for p in exclude_pops.split(',')]

    # 过滤群体|Filter populations
    filtered_pops = {}
    for pop, count in pop_counts.items():
        if pop in exclude_list:
            continue
        if count < min_samples:
            continue
        filtered_pops[pop] = count

    # 计算统计信息|Calculate statistics
    counts = list(filtered_pops.values())
    min_n = min(counts) if counts else 0
    max_n = max(counts) if counts else 0
    ratio = max_n / min_n if min_n > 0 else 0

    return {
        'pop_counts': filtered_pops,
        'min_n': min_n,
        'max_n': max_n,
        'ratio': ratio,
        'n_pops': len(filtered_pops)
    }


def decide_bootstrap_strategy(ratio: float, enable_bootstrap: bool = False) -> tuple:
    """
    根据样本量比值决定bootstrap策略
    Decide bootstrap strategy based on sample size ratio

    Args:
        ratio: 最大/最小样本量比值|Max/min sample size ratio
        enable_bootstrap: 是否强制启用bootstrap|Whether to force enable bootstrap

    Returns:
        (是否bootstrap, 迭代次数)|(should_bootstrap, n_iterations)
    """
    if enable_bootstrap:
        return True, 100

    if ratio < 3:
        return False, 1
    elif ratio < 5:
        return True, 50
    else:
        return True, 100


def decide_thinning_strategy(snp_count: int, thin_threshold: Optional[float] = None) -> Optional[float]:
    """
    根据SNP数量决定抽稀策略
    Decide thinning strategy based on SNP count

    Args:
        snp_count: SNP数量|SNP count
        thin_threshold: 手动指定的抽稀阈值|Manually specified thinning threshold

    Returns:
        抽稀比例（None表示不抽稀）|Thinning ratio (None means no thinning)
    """
    if thin_threshold is not None:
        return thin_threshold

    if snp_count < 500000:
        return None
    elif snp_count < 2000000:
        return None  # 建议先做LD pruning
    elif snp_count < 5000000:
        return 0.5
    else:
        return 0.2


def filter_population_file(pop_file: str, output_file: str,
                           keep_pops: List[str]) -> int:
    """
    过滤群体文件，只保留指定的群体
    Filter population file, keep only specified populations

    Args:
        pop_file: 原始群体文件|Original population file
        output_file: 输出文件路径|Output file path
        keep_pops: 要保留的群体列表|List of populations to keep

    Returns:
        保留的样本数|Number of samples kept
    """
    count = 0
    with open(pop_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.strip().startswith('#'):
                continue
            parts = line.strip().split('\t') if '\t' in line else line.strip().split(',')
            if len(parts) >= 2 and parts[1] in keep_pops:
                fout.write(line)
                count += 1
    return count
