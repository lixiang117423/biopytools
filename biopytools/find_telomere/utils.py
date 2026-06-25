"""
端粒识别工具模块|Telomere Finder Utility Module
"""

import os
import subprocess
import logging
import sys
import time
import shutil
import re
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'tidk' or '/path/to/tidk')

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        # 例如: /miniforge3/envs/tidk/bin/tidk
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    # 尝试找到conda基础目录|Try to find conda base directory
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        # 需要获取envs目录|Need to get envs directory
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令|Search for command in all environments
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    # 找到了|Found it
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    # 检查是否在conda环境中|Check if in conda environment
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用 conda run|Use conda run
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 直接调用|Direct call
        full_cmd = [command] + args

    return full_cmd


class TelomereLogger:
    """端粒识别日志管理器|Telomere Finder Logger Manager"""

    def __init__(self, output_dir: str, verbose: bool = False, log_file: Optional[str] = None):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            verbose: 详细输出模式|Verbose output mode
            log_file: 日志文件路径|Log file path
        """
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        self.log_file = log_file
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """设置日志系统|Setup logging system"""
        logger = logging.getLogger('TelomereFinder')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO及以下级别|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler (如果指定)|File handler (if specified)
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取logger对象|Get logger object"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
            output_dir: 输出目录|Output directory
        """
        self.logger = logger
        self.output_dir = Path(output_dir)

    def run_command(self, command: List[str], description: str = "") -> bool:
        """
        执行命令|Run command

        Args:
            command: 命令列表（应由build_conda_command构建）|Command list (should be built by build_conda_command)
            description: 命令描述|Command description

        Returns:
            bool: 是否成功|Whether successful
        """
        if description:
            self.logger.info(f"执行中|Executing: {description}")

        cmd_str = ' '.join(command)
        self.logger.debug(f"完整命令|Full command: {cmd_str}")

        start_time = time.time()

        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )

            elapsed_time = time.time() - start_time

            if result.stdout:
                # 输出到stdout (INFO级别)|Output to stdout (INFO level)
                for line in result.stdout.strip().split('\n'):
                    if line:
                        self.logger.info(line)

            if result.stderr:
                # 输出到stderr (WARNING级别)|Output to stderr (WARNING level)
                for line in result.stderr.strip().split('\n'):
                    if line and 'Error' not in line and 'error' not in line:
                        self.logger.warning(line)

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed with return code {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error message: {result.stderr}")
                return False

            self.logger.info(f"命令执行完成|Command completed in {elapsed_time:.2f}秒|seconds")
            return True

        except Exception as e:
            elapsed_time = time.time() - start_time
            self.logger.error(f"命令执行异常|Command execution exception: {str(e)}")
            self.logger.error(f"执行时间|Execution time: {elapsed_time:.2f}秒|seconds")
            return False


class OutputValidator:
    """输出验证器|Output Validator"""

    def __init__(self, logger: logging.Logger):
        """
        初始化验证器|Initialize validator

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def validate_output_files(self, expected_files: List[str]) -> bool:
        """
       验证输出文件|Validate output files

        Args:
            expected_files: 期望的输出文件列表|Expected output file list

        Returns:
            bool: 是否所有文件都存在|Whether all files exist
        """
        all_exist = True
        for file_path in expected_files:
            if Path(file_path).exists():
                self.logger.info(f"输出文件已生成|Output file generated: {file_path}")
            else:
                self.logger.warning(f"输出文件未生成|Output file not generated: {file_path}")
                all_exist = False

        return all_exist

    def get_file_stats(self, file_path: str) -> dict:
        """
        获取文件统计信息|Get file statistics

        Args:
            file_path: 文件路径|File path

        Returns:
            dict: 文件统计信息|File statistics
        """
        path = Path(file_path)
        if not path.exists():
            return {'exists': False}

        stats = {
            'exists': True,
            'size_mb': path.stat().st_size / (1024 * 1024),
            'size_bytes': path.stat().st_size
        }

        # 如果是文本文件，计算行数|If text file, count lines
        if path.suffix in ['.tsv', '.txt', '.bedgraph']:
            try:
                with open(file_path, 'r') as f:
                    stats['lines'] = sum(1 for _ in f)
            except Exception as e:
                self.logger.warning(f"无法读取文件|Cannot read file: {e}")

        return stats


def get_reverse_complement(sequence: str) -> str:
    """
    获取DNA序列的反向互补序列|Get reverse complement of DNA sequence

    Args:
        sequence: DNA序列|DNA sequence

    Returns:
        str: 反向互补序列|Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))


class CladeDatabase:
    """端粒重复序列数据库|Telomeric Repeat Database"""

    # 常见物种的端粒重复序列|Common telomeric repeats for various clades
    CLADES = {
        # 昆虫|Insects
        'Coleoptera': ['TTAGG', 'TTAGGG'],
        'Hymenoptera': ['TTAGG'],
        'Lepidoptera': ['TTAGG'],
        'Diptera': ['TTAGG'],
        'Hemiptera': ['TTAGG'],
        'Orthoptera': ['TTAGG', 'TTAGGG'],
        'Odonata': ['TTAGG'],
        'Plecoptera': ['TTAGG'],
        'Trichoptera': ['TTAGG'],
        'Symphypleona': ['TTAGG'],

        # 脊椎动物|Vertebrates
        'Mammalia': ['TTAGGG'],
        'Carnivora': ['TTAGGG'],
        'Rodentia': ['TTAGGG'],
        'Chiroptera': ['TTAGGG'],
        'Aves': ['TTAGGG'],
        'Accipitriformes': ['TTAGGG'],
        'Anura': ['TTAGGG'],
        'Caprimulgiformes': ['TTAGGG'],
        'Actinopterygii': ['TTAGGG'],
        'Perciformes': ['TTAGGG'],
        'Salmoniformes': ['TTAGGG'],
        'Cypriniformes': ['TTAGGG'],
        'Labriformes': ['TTAGGG'],
        'Syngnathiformes': ['TTAGGG'],
        'Carangiformes': ['TTAGGG'],
        'Pleuronectiformes': ['TTAGGG'],
        'Carcharhiniformes': ['TTAGGG'],

        # 植物|Plants
        'Arabidopsis': ['TTTAGGG'],
        'Poales': ['TTTAGGG'],
        'Rosales': ['TTTAGGG'],
        'Fabales': ['TTTAGGG'],
        'Malpighiales': ['TTTAGGG'],
        'Myrtales': ['TTTAGGG'],
        'Sapindales': ['TTTAGGG'],
        'Caryophyllales': ['TTTAGGG'],
        'Asterales': ['TTTAGGG'],
        'Lamiales': ['TTTAGGG'],
        'Solanales': ['TTTAGGG'],
        'Apiales': ['TTTAGGG'],
        'Fagales': ['TTTAGGG'],
        'Buxales': ['TTTAGGG'],
        'Hypnales': ['TTTAGGG'],
        'Chlamydomonadales': ['TTTAGGG'],

        # 其他生物|Other organisms
        'Nematoda': ['TTAGGC'],
        'Arachnida': ['TTAGGG'],
        'Crustacea': ['TTAGGG'],
        'Fungi': ['TTAGGG'],
    }

    @classmethod
    def get_clades(cls) -> List[str]:
        """获取所有支持的分类群|Get all supported clades"""
        return list(cls.CLADES.keys())

    @classmethod
    def get_repeats_for_clade(cls, clade: str) -> Optional[List[str]]:
        """获取指定分类群的端粒重复序列|Get telomeric repeats for specified clade"""
        return cls.CLADES.get(clade)

    @classmethod
    def print_clade_table(cls):
        """打印分类群表格|Print clade table"""
        print("支持的分类群及其端粒重复序列|Supported clades and their telomeric repeats:")
        print("=" * 80)
        print(f"{'分类群|Clade':<25} {'端粒重复序列|Telomeric Repeat'}")
        print("=" * 80)
        for clade, repeats in sorted(cls.CLADES.items()):
            print(f"{clade:<25} {', '.join(repeats)}")
        print("=" * 80)


def evaluate_telomere_search_results(windows_file: str,
                                     logger: Optional[logging.Logger] = None) -> dict:
    """
    评估端粒搜索结果的质量|Evaluate quality of telomere search results

    Args:
        windows_file: tidk search输出的windows文件路径|tidk search output windows file path
        logger: 日志对象(可选)|Logger object (optional)

    Returns:
        dict: 评估结果|Evaluation results
            - valid: 是否有效|Whether valid
            - total_windows: 端粒窗口总数|Total telomere windows
            - chromosomes: 涉及的染色体数量|Number of chromosomes involved
            - reason: 原因说明|Reason description
    """
    if not os.path.exists(windows_file):
        if logger:
            logger.warning(f"端粒窗口文件不存在|Telomere windows file not found: {windows_file}")
        return {
            'valid': False,
            'total_windows': 0,
            'chromosomes': 0,
            'reason': 'File not found'
        }

    try:
        with open(windows_file, 'r') as f:
            lines = f.readlines()

        # 跳过可能的注释行和header|Skip possible comment lines and header
        data_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#'):
                data_lines.append(line)

        total_windows = len(data_lines)

        # 提取涉及的染色体|Extract chromosomes involved
        chromosomes = set()
        for line in data_lines:
            parts = line.split('\t')
            if len(parts) >= 1:
                chromosomes.add(parts[0])

        num_chromosomes = len(chromosomes)

        # 评估规则|Evaluation rules
        min_windows = 10  # 最少需要10个端粒窗口|Minimum 10 telomere windows
        min_chromosomes = 1  # 至少1条染色体|At least 1 chromosome

        valid = total_windows >= min_windows and num_chromosomes >= min_chromosomes

        reason = ""
        if valid:
            reason = f"Found {total_windows} telomere windows across {num_chromosomes} chromosome(s)"
        else:
            if total_windows < min_windows:
                reason = f"Insufficient telomere windows ({total_windows} < {min_windows})"
            elif num_chromosomes < min_chromosomes:
                reason = f"Insufficient chromosome coverage ({num_chromosomes} < {min_chromosomes})"

        if logger:
            status = "VALID" if valid else "INVALID"
            logger.info(f"端粒搜索结果评估|Telomere search evaluation: {status}")
            logger.info(f"  端粒窗口数|Telomere windows: {total_windows}")
            logger.info(f"  涉及染色体|Chromosomes: {num_chromosomes}")
            logger.info(f"  原因|Reason: {reason}")

        return {
            'valid': valid,
            'total_windows': total_windows,
            'chromosomes': num_chromosomes,
            'reason': reason
        }

    except Exception as e:
        if logger:
            logger.error(f"评估端粒搜索结果时出错|Error evaluating telomere search results: {str(e)}")
        return {
            'valid': False,
            'total_windows': 0,
            'chromosomes': 0,
            'reason': f'Error: {str(e)}'
        }


def get_plant_telomere_priority_list() -> List[str]:
    """
    获取植物端粒序列的优先级列表|Get priority list of plant telomere repeats

    返回按优先级排序的植物端粒序列列表
    Returns plant telomere repeats sorted by priority

    Returns:
        List[str]: 植物端粒序列优先级列表|Plant telomere repeat priority list
    """
    # 第一优先级: 典型植物端粒|First priority: Typical plant telomere
    priority_1 = ['TTTAGGG', 'CCCTAAA']

    # 第二优先级: 其他常见植物端粒|Second priority: Other common plant telomeres
    plant_clades = {
        'Arabidopsis', 'Poales', 'Rosales', 'Fabales', 'Malpighiales',
        'Myrtales', 'Sapindales', 'Caryophyllales', 'Asterales', 'Lamiales',
        'Solanales', 'Apiales', 'Fagales', 'Buxales', 'Hypnales',
        'Chlamydomonadales'
    }

    priority_2 = []
    for clade in plant_clades:
        repeats = CladeDatabase.get_repeats_for_clade(clade)
        if repeats:
            priority_2.extend(repeats)

    # 去重并保持顺序|Deduplicate while maintaining order
    seen = set()
    unique_priority_2 = []
    for seq in priority_2:
        seq_upper = seq.upper()
        if seq_upper not in seen and seq_upper not in priority_1:
            seen.add(seq_upper)
            unique_priority_2.append(seq_upper)

    return priority_1 + unique_priority_2


def select_prioritized_telomere_repeat(candidates: List[str],
                                       logger: Optional[logging.Logger] = None) -> str:
    """
    从候选端粒序列中优先选择已知的植物端粒序列|Select prioritized telomere repeat from candidates

    优先级策略|Priority strategy:
    1. 植物端粒序列(TTTAGGG及其反向互补)|Plant telomere repeat (TTTAGGG and its reverse complement)
    2. 常见动物端粒序列(TTAGGG及其反向互补)|Common animal telomere repeat (TTAGGG and its reverse complement)
    3. 候选列表中的第一个序列|First sequence in the candidate list

    Args:
        candidates: 候选端粒序列列表|Candidate telomere repeat list
        logger: 日志对象(可选)|Logger object (optional)

    Returns:
        str: 选择的端粒序列|Selected telomere repeat
    """
    if not candidates:
        raise ValueError("候选序列列表不能为空|Candidate list cannot be empty")

    # 已知端粒序列优先级列表|Known telomere repeat priority list
    priority_repeats = [
        'TTTAGGG',  # 植物端粒序列|Plant telomere repeat
        'CCCTAAA',  # TTTAGGG的反向互补|Reverse complement of TTTAGGG
        'TTAGGG',   # 脊椎动物端粒序列|Vertebrate telomere repeat
        'CCCTAA',   # TTAGGG的反向互补|Reverse complement of TTAGGG
    ]

    # 标准化候选序列列表(去除空格,转大写)|Normalize candidate list (remove spaces, uppercase)
    normalized_candidates = [seq.strip().upper() for seq in candidates if seq.strip()]

    if logger:
        logger.debug(f"候选端粒序列|Candidate telomere repeats: {normalized_candidates}")

    # 在候选列表中查找优先级高的序列|Search for high-priority sequences in candidates
    for priority_seq in priority_repeats:
        for candidate in normalized_candidates:
            # 检查是否匹配或包含|Check if matches or contains
            if candidate == priority_seq or candidate in priority_seq or priority_seq in candidate:
                if logger:
                    logger.info(f"优先选择已知端粒序列|Prioritizing known telomere repeat: {candidate}")
                return candidate

    # 如果没有找到已知序列,返回第一个候选序列|If no known sequence found, return first candidate
    if logger:
        logger.info(f"未发现已知端粒序列,使用第一个候选序列|No known telomere repeat found, using first candidate: {normalized_candidates[0]}")

    return normalized_candidates[0]
