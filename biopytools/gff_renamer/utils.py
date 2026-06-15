"""
GFF重命名工具函数模块|GFF Renamer Utility Functions Module
"""

import logging
import subprocess
import sys
import re
import shutil
import os
from pathlib import Path
from typing import Optional, Tuple, Dict, List


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name"""
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
    """构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment"""
    # AGAT是Perl工具，用basename检测conda环境
    cmd_name = os.path.basename(command)
    conda_env = get_conda_env(cmd_name)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class GFFRenamerLogger:
    """GFF重命名日志管理器|GFF Renamer Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 设置日志级别|Set log level
        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        # 配置logging|Configure logging
        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def parse_chromosome(seq_id: str, chr_mapping: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """
    从序列ID解析染色体信息|Parse chromosome information from sequence ID

    Args:
        seq_id: 序列ID|Sequence ID (e.g., OV12_RagTag, Chr1, chr01, scaffold_1)
        chr_mapping: 染色体映射字典|Chromosome mapping dictionary (optional)

    Returns:
        tuple: (标准化染色体ID, 染色体编号)|(Standardized chromosome ID, chromosome number)
            例如|e.g.: ("Chr12", "12")
    """
    # 如果提供了映射文件，优先使用映射|Use mapping file if provided
    if chr_mapping and seq_id in chr_mapping:
        mapped_id = chr_mapping[seq_id]
        # 从映射后的ID提取编号|Extract number from mapped ID
        match = re.search(r'(\d+)', mapped_id)
        if match:
            chr_num = match.group(1)
            return mapped_id, f"{int(chr_num):02d}"
        return mapped_id, "00"

    # 移除常见的版本信息|Remove common version info
    seq_id = re.sub(r'_RagTag$', '', seq_id)
    seq_id = re.sub(r'\.v\d+$', '', seq_id)
    seq_id = re.sub(r'_primary$', '', seq_id)
    seq_id = re.sub(r'_alternate$', '', seq_id)

    # 尝试提取染色体编号|Try to extract chromosome number
    match = re.search(r'(\d+)', seq_id)
    if match:
        chr_num = match.group(1)

        seq_id_lower = seq_id.lower()
        has_chr = 'chr' in seq_id_lower
        has_scaffold = 'scaffold' in seq_id_lower
        has_contig = 'contig' in seq_id_lower

        # 从关键字后提取编号，避免被前缀中的数字干扰
        # Extract number after keyword to avoid interference from numbers in prefix
        if has_chr:
            chr_match = re.search(r'[Cc]hr\D*(\d+)', seq_id)
            if chr_match:
                chr_num = chr_match.group(1)
        elif has_scaffold or has_contig:
            sc_match = re.search(r'(?:scaffold|contig)[_.]*(\d+)', seq_id, re.IGNORECASE)
            if sc_match:
                chr_num = sc_match.group(1)

        # 格式化为两位数|Format to two digits
        chr_num_formatted = f"{int(chr_num):02d}"

        if has_chr:
            return f"Chr{chr_num_formatted}", chr_num_formatted
        elif has_scaffold:
            # scaffold使用Scf前缀避免与染色体编号歧义
            # Use Scf prefix for scaffold to avoid ambiguity with chromosome numbers
            return f"Scf{chr_num_formatted}", f"Scf{chr_num_formatted}"
        elif has_contig:
            # contig使用Ctg前缀避免与染色体编号歧义
            # Use Ctg prefix for contig to avoid ambiguity with chromosome numbers
            return f"Ctg{chr_num_formatted}", f"Ctg{chr_num_formatted}"
        else:
            # 默认格式|Default format
            return f"Chr{chr_num_formatted}", chr_num_formatted

    # 如果无法提取编号，返回原始ID|If cannot extract number, return original ID
    return seq_id, "00"


def load_chr_mapping(mapping_file: str) -> Dict[str, str]:
    """
    加载染色体映射文件|Load chromosome mapping file

    Args:
        mapping_file: 映射文件路径|Mapping file path

    Returns:
        dict: {原始seqid: 标准化seqid}|Dictionary mapping {original_seqid: standardized_seqid}

    文件格式|File format:
        chr1    Chr1
        chr2    Chr2
        scaffold_1    Scf1
    """
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def generate_gene_id(prefix: str, species: str, chr_num: str, gene_index: int,
                     naming_format: str = "standard") -> str:
    """
    生成标准化的基因ID|Generate standardized gene ID

    Args:
        prefix: 前缀|Prefix (e.g., CDRT)
        species: 物种缩写|Species abbreviation (e.g., Ov)
        chr_num: 染色体编号|Chromosome number (e.g., 12, 01)
        gene_index: 基因索引|Gene index (1-based)
        naming_format: 命名格式|Naming format (standard/simple/compact)

    Returns:
        str: 基因ID|Gene ID
            - standard: CDRT_Ov12g000010
            - simple: CDRT12G000010
            - compact: CDRT12g000010
    """
    # 基因编号格式：000010, 000020, ... (递增10)|Gene number format: increments by 10
    gene_number = gene_index * 10

    if naming_format == "standard":
        # 标准格式|Standard format: {prefix}_{species}{chr_num}g{number:06d}
        gene_id = f"{prefix}_{species}{chr_num}g{gene_number:06d}"
    elif naming_format == "simple":
        # 简化格式|Simple format: {prefix}{chr_num}G{number:06d}
        gene_id = f"{prefix}{chr_num}G{gene_number:06d}"
    elif naming_format == "compact":
        # 紧凑格式|Compact format: {prefix}{chr_num}g{number:06d}
        gene_id = f"{prefix}{chr_num}g{gene_number:06d}"
    else:
        # 默认使用标准格式|Default to standard format
        gene_id = f"{prefix}_{species}{chr_num}g{gene_number:06d}"

    return gene_id
