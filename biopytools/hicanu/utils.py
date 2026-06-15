"""
HiCanu组装工具函数模块|HiCanu Assembly Utility Functions Module
"""

import logging
import os
import sys
import subprocess
import shutil
import re
from datetime import datetime
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'canu' or '/path/to/canu')

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        # 例如: /miniforge3/envs/canu/bin/canu
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
        full_cmd = ['conda', 'run', '-n', conda_env, command] + args
    else:
        # 直接调用|Direct call
        full_cmd = [command] + args

    return full_cmd


class CanuLogger:
    """HiCanu组装日志管理器|HiCanu Assembly Logger Manager"""

    def __init__(self, output_dir: str, module_name: str = 'biopytools.hicanu'):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录路径|Output directory path
            module_name: 模块名称|Module name
        """
        self.output_dir = Path(output_dir)
        self.module_name = module_name

        # 创建日志文件路径|Create log file paths
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f'hicanu_{timestamp}.log'

        # 配置日志|Configure logging
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """配置日志系统|Setup logging system"""
        logger = logging.getLogger(self.module_name)
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO及以下级别|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志对象|Get logger object"""
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

    def run_command(self, cmd: list, dry_run: bool = False,
                   timeout: Optional[int] = None) -> bool:
        """
        运行命令|Run command

        Args:
            cmd: 命令列表|Command list
            dry_run: 模拟运行|Dry run mode
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            bool: 是否成功|Whether successful
        """
        cmd_str = ' '.join(cmd)

        if dry_run:
            self.logger.info(f"[DRY RUN] 待执行命令|Command to run: {cmd_str}")
            return True

        self.logger.info(f"执行命令|Executing command: {cmd_str}")
        self.logger.debug(f"完整命令|Full command: {cmd}")

        try:
            # 执行命令|Execute command
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )

            # 记录输出|Log output
            if result.stdout:
                self.logger.debug(f"命令输出|Command output:\n{result.stdout}")

            if result.stderr:
                self.logger.warning(f"命令错误输出|Command stderr:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(
                    f"命令执行失败|Command execution failed with return code {result.returncode}"
                )
                return False

            self.logger.info("命令执行成功|Command executed successfully")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command execution timeout after {timeout} seconds")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {str(e)}")
            return False

    def check_output_exists(self, filepath: str, error_on_missing: bool = False, silent: bool = False) -> bool:
        """
        检查输出文件是否存在|Check if output file exists

        Args:
            filepath: 文件路径|File path
            error_on_missing: 文件不存在时是否记录ERROR级别日志（否则WARNING或INFO）|Whether to log ERROR level when file is missing (otherwise WARNING or INFO)
            silent: 是否静默模式（不记录任何日志）|Whether to use silent mode (no logging)

        Returns:
            bool: 文件是否存在|Whether file exists
        """
        if not os.path.exists(filepath):
            if not silent:
                if error_on_missing:
                    self.logger.error(f"输出文件不存在|Output file does not exist: {filepath}")
                else:
                    # 对于可选文件使用INFO级别|Use INFO level for optional files
                    self.logger.info(f"输出文件不存在（可选文件）|Output file not found (optional): {filepath}")
            return False

        file_size = os.path.getsize(filepath)
        self.logger.info(f"输出文件已生成|Output file generated: {filepath} ({file_size:,} bytes)")
        return True


def parse_fasta_stats(fasta_file: str, logger: logging.Logger) -> dict:
    """
    解析FASTA文件统计信息|Parse FASTA file statistics

    Args:
        fasta_file: FASTA文件路径|FASTA file path
        logger: 日志对象|Logger object

    Returns:
        dict: 统计信息|Statistics
    """
    logger.info(f"解析FASTA文件|Parsing FASTA file: {fasta_file}")

    if not os.path.exists(fasta_file):
        logger.error(f"FASTA文件不存在|FASTA file does not exist: {fasta_file}")
        return {}

    stats = {
        'total_sequences': 0,
        'total_bp': 0,
        'min_length': float('inf'),
        'max_length': 0,
        'n50': 0,
        'l50': 0
    }

    lengths = []

    try:
        with open(fasta_file, 'r') as f:
            current_length = 0
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_length > 0:
                        lengths.append(current_length)
                    current_length = 0
                    stats['total_sequences'] += 1
                else:
                    current_length += len(line)

            # 添加最后一条序列|Add last sequence
            if current_length > 0:
                lengths.append(current_length)

        if lengths:
            stats['total_bp'] = sum(lengths)
            stats['min_length'] = min(lengths)
            stats['max_length'] = max(lengths)
            stats['avg_length'] = stats['total_bp'] / stats['total_sequences']

            # 计算N50和L50|Calculate N50 and L50
            lengths_sorted = sorted(lengths, reverse=True)
            cumulative = 0
            for i, length in enumerate(lengths_sorted, 1):
                cumulative += length
                if cumulative >= stats['total_bp'] / 2:
                    stats['n50'] = length
                    stats['l50'] = i
                    break

        logger.info(f"FASTA统计信息|FASTA statistics:")
        logger.info(f"  序列数量|Number of sequences: {stats['total_sequences']:,}")
        logger.info(f"  总长度|Total bp: {stats['total_bp']:,}")
        logger.info(f"  N50: {stats['n50']:,}")
        logger.info(f"  L50: {stats['l50']:,}")
        logger.info(f"  最短序列|Min length: {stats['min_length']:,}")
        logger.info(f"  最长序列|Max length: {stats['max_length']:,}")

        return stats

    except Exception as e:
        logger.error(f"解析FASTA文件失败|Failed to parse FASTA file: {str(e)}")
        return {}


def format_size(size_bp: int) -> str:
    """
    格式化基因组大小|Format genome size

    Args:
        size_bp: 碱基对数量|Number of base pairs

    Returns:
        str: 格式化的大小|Formatted size
    """
    if size_bp >= 1_000_000_000:
        return f"{size_bp / 1_000_000_000:.2f} Gb"
    elif size_bp >= 1_000_000:
        return f"{size_bp / 1_000_000:.2f} Mb"
    elif size_bp >= 1_000:
        return f"{size_bp / 1_000:.2f} Kb"
    else:
        return f"{size_bp} bp"


def generate_contig_reads_map(raw_dir: str, prefix: str, fasta_dir: str, logger: logging.Logger) -> bool:
    """
    生成contig到reads的映射文件|Generate contig to reads mapping file

    为每个contig生成其组成的reads列表文件，格式: contig_id read_name
    For each contig, generate a file containing its constituent reads, format: contig_id read_name

    Args:
        raw_dir: Canu原始输出目录|Canu raw output directory
        prefix: 输出文件前缀|Output file prefix
        fasta_dir: 输出目录|Output directory
        logger: 日志对象|Logger object

    Returns:
        bool: 是否成功|Whether successful
    """
    logger.info("=" * 60)
    logger.info("生成Contig-Reads映射文件|Generating Contig-Reads mapping file")
    logger.info("=" * 60)

    # 文件路径|File paths
    read_to_tig_file = os.path.join(raw_dir, f'{prefix}.contigs.layout.readToTig')
    read_names_file = os.path.join(raw_dir, f'{prefix}.seqStore', 'readNames.txt')
    output_file = os.path.join(fasta_dir, f'{prefix}.contig_reads.tsv')

    # 检查必需文件是否存在|Check if required files exist
    if not os.path.exists(read_to_tig_file):
        logger.warning(f"readToTig文件不存在，跳过|readToTig file not found, skipping: {read_to_tig_file}")
        return False

    if not os.path.exists(read_names_file):
        logger.warning(f"readNames.txt文件不存在，跳过|readNames.txt not found, skipping: {read_names_file}")
        return False

    try:
        # 第一步：读取readNames.txt，建立索引到read名称的映射
        # Step 1: Read readNames.txt, build index to read name mapping
        logger.info(f"读取read名称映射|Reading read name mapping: {read_names_file}")
        read_id_to_name = {}
        with open(read_names_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    read_idx = int(parts[0])
                    read_name = parts[1]
                    read_id_to_name[read_idx] = read_name

        logger.info(f"  共读取|Total reads loaded: {len(read_id_to_name):,}")

        # 第二步：读取readToTig文件，生成映射
        # Step 2: Read readToTig file, generate mapping
        logger.info(f"读取readToTig映射|Reading readToTig mapping: {read_to_tig_file}")

        # 使用字典按tigID分组|Use dict to group by tigID
        tig_to_reads = {}
        read_count = 0

        with open(read_to_tig_file, 'r') as f:
            # 跳过注释行|Skip comment lines
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    read_idx = int(parts[0])
                    tig_id = parts[1]

                    if read_idx in read_id_to_name:
                        read_name = read_id_to_name[read_idx]

                        if tig_id not in tig_to_reads:
                            tig_to_reads[tig_id] = []
                        tig_to_reads[tig_id].append(read_name)
                        read_count += 1

        logger.info(f"  共处理|Total reads processed: {read_count:,}")
        logger.info(f"  共发现|Total contigs found: {len(tig_to_reads):,}")

        # 第三步：写入输出文件（按tigID排序）
        # Step 3: Write to output file (sorted by tigID)
        logger.info(f"写入映射文件|Writing mapping file: {output_file}")

        with open(output_file, 'w') as f_out:
            # 写入header|Write header
            f_out.write("#contig_id\tread_name\n")

            # 按tigID排序后写入|Write sorted by tigID (numeric sort)
            # tigID格式化为tig00000005（8位数字，前面补0）|Format tigID as tig00000005 (8-digit with leading zeros)
            for tig_id in sorted(tig_to_reads.keys(), key=lambda x: int(x) if x.isdigit() else x):
                # 将数字ID格式化为tig00000005格式|Format numeric ID as tig00000005
                contig_id = f"tig{int(tig_id):08d}"
                for read_name in tig_to_reads[tig_id]:
                    f_out.write(f"{contig_id}\t{read_name}\n")

        # 统计信息|Statistics
        logger.info("映射统计|Mapping statistics:")
        logger.info(f"  总contig数|Total contigs: {len(tig_to_reads):,}")

        # 计算每个contig的reads数量|Calculate reads per contig
        reads_per_tig = [len(reads) for reads in tig_to_reads.values()]
        if reads_per_tig:
            logger.info(f"  平均reads数/contig|Avg reads/contig: {sum(reads_per_tig) / len(reads_per_tig):.1f}")
            logger.info(f"  最少reads数|Min reads: {min(reads_per_tig):,}")
            logger.info(f"  最多reads数|Max reads: {max(reads_per_tig):,}")

        logger.info(f"  输出文件|Output file: {output_file}")
        logger.info("=" * 60)
        logger.info("Contig-Reads映射文件生成完成|Contig-Reads mapping file generated successfully")
        logger.info("=" * 60)

        return True

    except Exception as e:
        logger.error(f"生成contig-reads映射失败|Failed to generate contig-reads mapping: {str(e)}")
        return False
