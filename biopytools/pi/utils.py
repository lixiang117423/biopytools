"""
Pi计算工具函数模块|Pi Calculation Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class PiRow:
    """Pi计算结果行|Pi calculation result row"""
    population: str
    chromosome: str
    window_start: Optional[int]  # None for genome-wide mode
    window_end: Optional[int]    # None for genome-wide mode
    pi_value: float
    n_sites: int = 0
    source: str = ""  # "vcftools" or "pixy"


class PiLogger:
    """Pi计算日志管理器|Pi Calculation Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "pi_analysis.log"):
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

        # stderr handler - WARNING及以上|Stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler, stderr_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


# conda包装统一走公共层(规范:模块内禁止复制实现, §13.2.3传完整路径)
# |conda wrapping via common layer (no local copies; full paths, §13.2.3)
from ..common.conda_runner import build_conda_command, get_conda_env


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "", timeout: int = 3600) -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
            timeout: 超时时间(秒)|Timeout in seconds
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
                timeout=timeout,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr[:1000]}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command execution timeout: {description}")
            return False


def parse_population_file(pop_file: str) -> Dict[str, List[str]]:
    """
    解析群体文件，返回群体到样本列表的映射
    Parse population file, return mapping of population to sample list

    Args:
        pop_file: 群体文件路径|Population file path

    Returns:
        {群体名: [样本列表]}|{population_name: [sample_list]}
    """
    pop_to_samples: Dict[str, List[str]] = {}

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
                if population not in pop_to_samples:
                    pop_to_samples[population] = []
                pop_to_samples[population].append(sample_id)

    return pop_to_samples


def get_software_version(tool_path: str, logger) -> str:
    """
    自动检测软件版本|Auto-detect software version

    Args:
        tool_path: 工具路径|Tool path
        logger: 日志器|Logger

    Returns:
        版本字符串|Version string
    """
    try:
        cmd = build_conda_command(tool_path, ['--version'])
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            return version
        else:
            return "unknown"
    except Exception as e:
        logger.warning(f"版本检测失败|Version detection failed: {e}")
        return "unknown"


def parse_fai_file(fai_file: str) -> Dict[str, int]:
    """
    解析samtools faidx生成的.fai文件，获取各染色体长度
    Parse .fai file generated by samtools faidx, get chromosome lengths

    Args:
        fai_file: .fai文件路径|.fai file path

    Returns:
        {染色体名: 碱基长度}|{chromosome_name: base_length}
    """
    chrom_lengths: Dict[str, int] = {}

    with open(fai_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])

    return chrom_lengths


def parse_fasta_lengths(fasta_file: str) -> Dict[str, int]:
    """
    纯Python解析FASTA文件获取各染色体长度（降级路径，不依赖samtools）
    Parse FASTA in pure Python to get chromosome lengths (fallback, no samtools)

    Args:
        fasta_file: FASTA文件路径|FASTA file path

    Returns:
        {染色体名: 碱基长度}|{chromosome_name: base_length}
    """
    chrom_lengths: Dict[str, int] = {}
    current_chrom = None

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_chrom = line[1:].split()[0]
                chrom_lengths[current_chrom] = 0
            elif current_chrom is not None:
                chrom_lengths[current_chrom] += len(line)

    return chrom_lengths


def format_number(num: float) -> str:
    """格式化数字|Format number"""
    return f"{num:.6f}"


def is_step_completed(output_file: str) -> bool:
    """检查步骤是否已完成（通过输出文件存在性判断）|Check if step is done"""
    return Path(output_file).exists()


def ensure_tabix_index(vcf_file: str, vcftools_path: str, logger) -> bool:
    """
    检查VCF文件的tabix索引(.tbi)是否存在，不存在则自动创建
    Check if VCF tabix index (.tbi) exists, create if missing

    Args:
        vcf_file: VCF文件路径|VCF file path
        vcftools_path: vcftools路径（用于定位conda环境中的tabix）|vcftools path (to locate tabix in conda env)
        logger: 日志器|Logger

    Returns:
        True if index exists or was created successfully, False otherwise
    """
    tbi_file = vcf_file + '.tbi'

    if os.path.exists(tbi_file):
        logger.info(f"VCF索引已存在|VCF index already exists: {tbi_file}")
        return True

    logger.info(f"VCF索引不存在，正在创建|VCF index not found, creating: {tbi_file}")

    # 使用vcftools所在conda环境中的tabix（同一环境通常都有htslib/tabix）
    # Use tabix from the same conda env as vcftools
    tabix_path = vcftools_path.replace('vcftools', 'tabix')
    if not os.path.exists(tabix_path):
        # 回退到直接用vcftools路径所在目录找tabix
        # Fallback: look for tabix in the same directory
        tabix_path = 'tabix'

    cmd = build_conda_command(tabix_path, ['-p', 'vcf', vcf_file])
    try:
        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            timeout=7200  # 大文件建索引可能需要较长时间
        )

        if os.path.exists(tbi_file):
            logger.info(f"VCF索引创建成功|VCF index created successfully: {tbi_file}")
            return True
        else:
            logger.error(f"tabix命令执行成功但索引文件未生成|tabix command succeeded but index file not created")
            return False

    except subprocess.CalledProcessError as e:
        logger.error(f"创建VCF索引失败|Failed to create VCF index")
        logger.error(f"错误代码|Error code: {e.returncode}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr[:1000]}")
        return False
    except subprocess.TimeoutExpired:
        logger.error(f"创建VCF索引超时|VCF index creation timeout")
        return False


def ensure_fai_index(fasta_path: str, samtools_path: str, logger) -> bool:
    """
    检查FASTA的.fai索引(.fai)是否存在，不存在则用samtools faidx自动创建
    Check if FASTA .fai index exists, auto-create with samtools faidx if missing

    Args:
        fasta_path: FASTA文件路径|FASTA file path
        samtools_path: samtools路径|samtools path
        logger: 日志器|Logger

    Returns:
        True if index exists or was created successfully, False otherwise
    """
    fai_file = fasta_path + '.fai'

    if os.path.exists(fai_file):
        logger.info(f"fai索引已存在|fai index already exists: {fai_file}")
        return True

    logger.info(f"fai索引不存在，正在创建|fai index not found, creating: {fai_file}")

    cmd = build_conda_command(samtools_path, ['faidx', fasta_path])
    # 记录完整命令（INFO级别）|Log complete command (INFO level)
    logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            timeout=3600  # 大基因组建索引可能较慢|Large genomes may take longer
        )

        if os.path.exists(fai_file):
            logger.info(f"fai索引创建成功|fai index created successfully: {fai_file}")
            return True
        else:
            logger.error(f"samtools faidx执行成功但索引文件未生成|samtools faidx succeeded but index file not created")
            return False

    except FileNotFoundError as e:
        logger.warning(f"samtools不存在，无法创建fai索引|samtools not found, cannot create fai index: {e.filename}")
        return False
    except subprocess.CalledProcessError as e:
        logger.warning(f"创建fai索引失败|Failed to create fai index")
        logger.warning(f"错误代码|Error code: {e.returncode}")
        if e.stderr:
            logger.warning(f"错误信息|Error message: {e.stderr[:1000]}")
        return False
    except subprocess.TimeoutExpired:
        logger.warning(f"创建fai索引超时|fai index creation timeout")
        return False
