"""
RagTag工具函数模块|RagTag Utility Functions Module
"""

import logging
import sys
import subprocess
from pathlib import Path


class RagTagLogger:
    """RagTag日志管理器|RagTag Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, output_dir):
        """初始化命令执行器|Initialize command runner"""
        self.logger = logger
        self.output_dir = Path(output_dir)

    def run_command(self, command, check=True):
        """
        执行命令|Execute command

        Args:
            command: 命令列表|Command list
            check: 是否检查返回码|Whether to check return code

        Returns:
            subprocess.CompletedProcess: 命令执行结果|Command execution result
        """
        self.logger.info(f"执行命令|Executing command: {' '.join(command)}")

        try:
            result = subprocess.run(
                command,
                check=check,
                capture_output=True,
                text=True
            )

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            if result.stderr:
                self.logger.debug(f"标准错误|Stderr: {result.stderr}")

            return result

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            raise

    def run_ragtag_scaffold(self, config):
        """
        运行RagTag scaffold命令|Run RagTag scaffold command

        Args:
            config: RagTagConfig对象|RagTagConfig object

        Returns:
            subprocess.CompletedProcess: 命令执行结果|Command execution result
        """
        # 构建基础命令|Build base command
        cmd = [config.ragtag_path, 'scaffold']

        # 添加输入文件|Add input files
        cmd.extend([config.reference, config.query])

        # 添加scaffolding选项|Add scaffolding options
        cmd.extend(['-t', str(config.threads)])
        cmd.extend(['-f', str(config.min_unique_length)])
        cmd.extend(['-q', str(config.min_mapq)])
        cmd.extend(['-d', str(config.max_merge_distance)])
        cmd.extend(['-i', str(config.min_grouping_confidence)])
        cmd.extend(['-a', str(config.min_location_confidence)])
        cmd.extend(['-s', str(config.min_orientation_confidence)])

        if config.concatenate_unplaced:
            cmd.append('-C')

        if config.infer_gaps:
            cmd.append('-r')
            cmd.extend(['-g', str(config.min_gap_size)])
            cmd.extend(['-m', str(config.max_gap_size)])

        if config.remove_small_alignments:
            cmd.append('--remove-small')

        # 添加比对器选项|Add aligner options
        cmd.extend(['--aligner', config.aligner])

        # 添加输出选项|Add output options
        cmd.extend(['-o', config.output_dir])

        return self.run_command(cmd)


def parse_fasta_header(header):
    """
    解析FASTA序列头|Parse FASTA sequence header

    Args:
        header: FASTA头（不包含>符号）|FASTA header (without > symbol)

    Returns:
        tuple: (序列ID, 描述信息)|Sequence ID and description
    """
    parts = header.split(None, 1)
    seq_id = parts[0]
    description = parts[1] if len(parts) > 1 else ''
    return seq_id, description


def rename_sequence_id(seq_id, sample_name):
    """
    重命名序列ID|Rename sequence ID

    转换规则|Conversion rules:
    - Chr12_RagTag -> Chr12
    - ptg000019l -> scaffold_1 (按顺序编号|numbered sequentially)

    Args:
        seq_id: 原始序列ID|Original sequence ID
        sample_name: 样品名称|Sample name

    Returns:
        str: 重命名后的序列ID|Renamed sequence ID
    """
    # 移除_RagTag后缀|Remove _RagTag suffix
    if seq_id.endswith('_RagTag'):
        return seq_id.replace('_RagTag', '')
    elif f"_{sample_name}" in seq_id:
        # 移除样品名称后缀（如果存在）|Remove sample name suffix (if exists)
        return seq_id.replace(f"_{sample_name}", '')
    else:
        # 保持原样或根据需要进行其他转换|Keep as is or perform other conversions as needed
        return seq_id


def is_scaffolded_sequence(seq_id):
    """
    判断序列是否已被scaffold|Check if sequence is scaffolded

    Args:
        seq_id: 序列ID|Sequence ID

    Returns:
        bool: 如果是scaffolded序列返回True，否则返回False
               Returns True if scaffolded, False otherwise

    判断规则|Judgment rules:
    1. 包含_RagTag后缀的序列是scaffolded的
    2. ptg/tig开头的序列是unscaffolded的
    3. chr/chrom开头的序列是scaffolded的
    """
    seq_id_lower = seq_id.lower()

    # 包含_RagTag后缀的序列是scaffolded的
    # Sequences with _RagTag suffix are scaffolded
    if '_ragtag' in seq_id_lower:
        return True

    # ptg/tig开头的序列是unscaffolded的 (RagTag默认的未定位序列前缀)
    # Sequences starting with ptg/tig are unscaffolded (RagTag default prefix)
    if seq_id_lower.startswith('ptg') or seq_id_lower.startswith('tig'):
        return False

    # chr/chrom开头的序列是scaffolded的
    # Sequences starting with chr/chrom are scaffolded
    if seq_id_lower.startswith('chr') or seq_id_lower.startswith('chrom'):
        return True

    # 其他情况：如果是scaffold_开头的（已重命名的），认为是unscaffolded
    # Otherwise: if starts with scaffold_ (already renamed), consider unscaffolded
    if seq_id_lower.startswith('scaffold_'):
        return False

    # 默认情况下认为是unscaffolded的 (更保守的策略)
    # By default, consider unscaffolded (more conservative strategy)
    return False
