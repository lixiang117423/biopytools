"""
BRAKER3基因组注释工具函数模块|BRAKER3 Genome Annotation Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import glob
import re
import shutil
import time
from pathlib import Path
from typing import List, Tuple, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 优先从传入的完整路径中检测|First check the passed path directly
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 再尝试从which结果检测|Then try from which result
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            command_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command_name)
                if os.path.exists(env_bin):
                    return env_name

    return None


class BrakerLogger:
    """BRAKER3注释日志管理器|BRAKER3 Annotation Logger Manager"""

    def __init__(self, log_file_path: str, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file_path: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = Path(log_file_path)
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # 创建格式化器|Create formatter
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 创建handlers列表|Create handlers list
        handlers = []

        # 添加文件handler|Add file handler
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

        # 配置stdout handler|Configure stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.setFormatter(formatter)
        handlers.append(stdout_handler)

        # 配置根日志记录器|Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)

        # 清除现有的handlers|Clear existing handlers
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

        # 添加新的handlers|Add new handlers
        for handler in handlers:
            root_logger.addHandler(handler)

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = None):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger object
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = working_dir if working_dir else "."

    def run_command(self, cmd: str, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description

        Returns:
            bool: 是否成功|Whether successful
        """
        # 自动包装conda环境的命令|Automatically wrap conda environment commands
        wrapped_cmd = self._wrap_conda_command(cmd)

        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {wrapped_cmd}")

        try:
            result = subprocess.run(
                wrapped_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description or 'Command'}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")  # 只显示前500字符|Show first 500 chars

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description or 'Command'}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                # 将完整stderr写入文件，避免截断丢失关键错误信息
                # Write full stderr to file to avoid truncation losing critical error info
                import tempfile
                stderr_file = os.path.join(self.working_dir, f"stderr_{int(time.time())}.log") if self.working_dir else None
                if stderr_file:
                    try:
                        with open(stderr_file, 'w', encoding='utf-8') as f:
                            f.write(e.stderr)
                        self.logger.error(f"完整错误信息已写入|Full stderr written to: {stderr_file}")
                        # 日志中显示最后2000字符，通常包含实际报错
                        # Show last 2000 chars in log, usually contains the actual error
                        tail = e.stderr[-2000:] if len(e.stderr) > 2000 else e.stderr
                        self.logger.error(f"错误信息(尾部)|Error message (tail): {tail}")
                    except Exception:
                        self.logger.error(f"错误信息|Error message: {e.stderr[:5000]}")
                else:
                    self.logger.error(f"错误信息|Error message: {e.stderr[:5000]}")
            return False

        except Exception as e:
            self.logger.error(f"执行异常|Execution exception: {e}")
            return False

    def _wrap_conda_command(self, cmd: str) -> str:
        """自动检测并包装conda环境中的命令|Automatically detect and wrap conda environment commands"""
        # 提取命令的第一个词（命令可执行文件路径）|Extract first word (command executable path)
        cmd = cmd.strip()
        match = re.match(r'^(\S+)(.*)', cmd)
        if not match:
            return cmd

        command_exe = match.group(1)
        rest_args = match.group(2)

        # 检查是否在conda环境中|Check if in conda environment
        conda_env = get_conda_env(command_exe)

        if conda_env:
            # 提取纯命令名，conda run会自动在环境PATH中找到|Extract command name, conda run finds it in env PATH
            command_name = os.path.basename(command_exe)
            # 使用--no-capture-output避免conda缓冲输出导致内存问题
            # Use --no-capture-output to avoid conda buffering output causing memory issues
            return f"conda run -n {conda_env} --no-capture-output {command_name}{rest_args}"
        else:
            # 直接返回原命令|Return original command directly
            return cmd

    def run_singularity_command(self, cmd: str, description: str = "") -> bool:
        """
        执行Singularity容器内的命令|Execute command inside Singularity container

        Args:
            cmd: 容器内命令|Command inside container
            description: 步骤描述|Step description

        Returns:
            bool: 是否成功|Whether successful
        """
        return self.run_command(cmd, description)


def format_number(num: int) -> str:
    """
    格式化数字为大单位|Format number to large units

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def check_file_exists(file_path: str, logger=None) -> bool:
    """
    检查文件是否存在|Check if file exists

    Args:
        file_path: 文件路径|File path
        logger: 日志器|Logger object

    Returns:
        bool: 是否存在|Whether exists
    """
    exists = os.path.exists(file_path)
    if not exists and logger:
        logger.warning(f"文件不存在|File not found: {file_path}")
    return exists


def check_step_completed(output_dir: str, marker_files: list, logger=None) -> bool:
    """
    检查步骤是否完成|Check if step is completed

    Args:
        output_dir: 输出目录|Output directory
        marker_files: 标记文件列表|List of marker files
        logger: 日志器|Logger object

    Returns:
        bool: 是否完成|Whether completed
    """
    all_exist = True
    for marker_file in marker_files:
        file_path = os.path.join(output_dir, marker_file)
        if not os.path.exists(file_path):
            all_exist = False
            if logger:
                logger.debug(f"步骤未完成，缺少文件|Step incomplete, missing file: {marker_file}")
            break

    if all_exist and logger:
        logger.info(f"检测到已完成步骤|Detected completed step in: {output_dir}")

    return all_exist


def find_long_reads_in_directory(directory: str, logger=None) -> Optional[str]:
    """
    在目录中查找三代全长转录本文件|Find long-read transcript files in directory

    Args:
        directory: 目录路径|Directory path
        logger: 日志器|Logger object

    Returns:
        str: 找到的文件路径，如果没有找到则返回 None|Found file path or None if not found
    """
    if not os.path.isdir(directory):
        if logger:
            logger.warning(f"三代转录组路径不是目录|Long-read directory is not a directory: {directory}")
        return None

    # 常见的三代转录组文件模式|Common long-read file patterns
    patterns = [
        "*.fastq", "*.fq", "*.fasta", "*.fa",
        "*.fastq.gz", "*.fq.gz", "*.fasta.gz", "*.fa.gz",
        "*.iso.seq.fa", "*.isoseq.fa", "*IsoSeq.fa",
        "*nanopore*.fastq", "*nanopore*.fq"
    ]

    for pattern in patterns:
        files = list(Path(directory).glob(pattern))
        if files:
            if logger:
                logger.info(f"找到三代转录组文件|Found long-read file: {files[0]}")
            return str(files[0])

    if logger:
        logger.warning(f"未找到三代转录组文件|No long-read files found in: {directory}")
    return None


def find_rnaseq_files_in_directory(directory: str, read1_pattern: str = "_1.fq.gz",
                                    read2_pattern: str = "_2.fq.gz", logger=None) -> List[Tuple[str, str]]:
    """
    在目录中查找二代RNA-seq配对文件|Find paired RNA-seq files in directory

    Args:
        directory: 目录路径|Directory path
        read1_pattern: R1文件模式|R1 file pattern
        read2_pattern: R2文件模式|R2 file pattern
        logger: 日志器|Logger object

    Returns:
        List[Tuple]: 配对文件列表 (R1, R2)|List of paired files (R1, R2)
    """
    if not os.path.isdir(directory):
        if logger:
            logger.warning(f"二代转录组路径不是目录|Short-read directory is not a directory: {directory}")
        return []

    paired_files = []

    # 查找所有可能的R1文件|Find all possible R1 files
    r1_files = list(Path(directory).glob(f"*{read1_pattern}"))

    if not r1_files:
        # 尝试其他常见模式|Try other common patterns
        alternative_patterns = ["_1.clean.fq.gz", "_R1.fastq.gz", "_R1.fq.gz", ".R1.fastq.gz"]
        for alt_pattern in alternative_patterns:
            r1_files = list(Path(directory).glob(f"*{alt_pattern}"))
            if r1_files:
                if logger:
                    logger.info(f"使用替代模式找到R1文件|Found R1 files with alternative pattern '{alt_pattern}'")
                read1_pattern = alt_pattern
                # 更新对应的R2模式|Update corresponding R2 pattern
                if alt_pattern == "_1.clean.fq.gz":
                    read2_pattern = "_2.clean.fq.gz"
                elif alt_pattern == "_R1.fastq.gz":
                    read2_pattern = "_R2.fastq.gz"
                elif alt_pattern == "_R1.fq.gz":
                    read2_pattern = "_R2.fq.gz"
                elif alt_pattern == ".R1.fastq.gz":
                    read2_pattern = ".R2.fastq.gz"
                break

    for r1_file in r1_files:
        # 获取样本名（去除R1模式）|Get sample name (remove R1 pattern)
        sample_name = str(r1_file.name).replace(read1_pattern, "")

        # 构建R2文件路径|Build R2 file path
        r2_file = r1_file.parent / f"{sample_name}{read2_pattern}"

        if r2_file.exists():
            paired_files.append((str(r1_file), str(r2_file)))
            if logger:
                logger.debug(f"找到配对文件|Found paired files: {sample_name}")
        else:
            if logger:
                logger.warning(f"未找到配对的R2文件|R2 file not found for: {r1_file.name}")

    if logger:
        logger.info(f"共找到 {len(paired_files)} 对RNA-seq文件|Found {len(paired_files)} pairs of RNA-seq files")

    return paired_files


def find_protein_files_in_directory(directory: str, logger=None) -> Optional[str]:
    """
    在目录中查找蛋白质序列文件|Find protein sequence files in directory

    Args:
        directory: 目录路径|Directory path
        logger: 日志器|Logger object

    Returns:
        str: 合并后的蛋白质文件路径，如果没有找到则返回 None|Path to merged protein file or None if not found
    """
    if not os.path.isdir(directory):
        # 如果是文件，直接返回|If it's a file, return directly
        if os.path.isfile(directory):
            if logger:
                logger.info(f"使用蛋白质文件|Using protein file: {directory}")
            return directory
        if logger:
            logger.warning(f"蛋白质序列路径不存在|Protein directory does not exist: {directory}")
        return None

    # 常见的蛋白质文件模式|Common protein file patterns
    patterns = [
        "*.fa", "*.faa", "*.fasta", "*.pep", "*.pep.fa", "*.protein.fa",
        "*.fa.gz", "*.faa.gz", "*.fasta.gz", "*.pep.gz"
    ]

    protein_files = []
    for pattern in patterns:
        files = list(Path(directory).glob(pattern))
        protein_files.extend(files)

    if not protein_files:
        if logger:
            logger.warning(f"未找到蛋白质序列文件|No protein files found in: {directory}")
        return None

    if len(protein_files) == 1:
        if logger:
            logger.info(f"找到蛋白质文件|Found protein file: {protein_files[0]}")
        # 即使是单个文件，也检查是否有重复ID
        # Even for single file, check for duplicate IDs
        protein_file = str(protein_files[0])
        return fix_duplicate_protein_ids(protein_file, logger)

    # 如果找到多个文件，需要合并|If multiple files found, need to merge
    if logger:
        logger.info(f"找到 {len(protein_files)} 个蛋白质文件，将进行合并|Found {len(protein_files)} protein files, will merge")

    import tempfile
    # 创建临时合并文件|Create temporary merged file
    temp_merged = tempfile.NamedTemporaryFile(mode='w', suffix='_merged_proteins.fa', delete=False, dir=str(directory))
    temp_merged_path = temp_merged.name
    temp_merged.close()

    try:
        # 合并所有蛋白质文件|Merge all protein files
        with open(temp_merged_path, 'w') as outfile:
            for protein_file in protein_files:
                if logger:
                    logger.debug(f"合并文件|Merging file: {protein_file.name}")

                # 读取并写入|Read and write
                with open(protein_file, 'r') as infile:
                    data = infile.read()
                    # 如果文件不以 > 开头，添加标题|If file doesn't start with >, add header
                    if not data.startswith('>'):
                        outfile.write(f">{protein_file.stem}\n")
                    outfile.write(data)

                    if not data.endswith('\n'):
                        outfile.write('\n')

        if logger:
            logger.info(f"蛋白质文件合并完成|Protein files merged: {temp_merged_path}")

        # 检查并修复重复ID
        # Check and fix duplicate IDs
        fixed_path = fix_duplicate_protein_ids(temp_merged_path, logger)

        return fixed_path

    except Exception as e:
        if logger:
            logger.error(f"合并蛋白质文件失败|Failed to merge protein files: {e}")
        # 清理临时文件|Clean up temp file
        if os.path.exists(temp_merged_path):
            os.remove(temp_merged_path)
        return None


def fix_duplicate_gtf_transcript_ids(gtf_file: str, logger=None) -> bool:
    """
    修复GTF文件中重复的transcript_id|Fix duplicate transcript_ids in GTF file

    GeneMark-ET在某些情况下会产生重复的transcript_id，导致TSEBRA合并时失败。
    GeneMark-ET may produce duplicate transcript_ids in some cases, causing TSEBRA to fail.

    Args:
        gtf_file: GTF文件路径|Path to GTF file
        logger: 日志器|Logger object

    Returns:
        bool: 是否进行了修复|Whether fixes were applied
    """
    from collections import defaultdict

    if not os.path.exists(gtf_file):
        return False

    # 空文件无需处理|Empty file needs no processing
    if os.path.getsize(gtf_file) == 0:
        return False

    # 读取所有行，按gene_id分组，检测重复transcript_id
    # Read all lines, group by gene_id, detect duplicate transcript_ids
    lines = []
    with open(gtf_file, 'r') as f:
        for line in f:
            lines.append(line)

    # 检查是否有重复的transcript_id
    # Check if there are duplicate transcript_ids
    tid_counts = defaultdict(int)
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        match = re.search(r'transcript_id "([^"]+)"', line)
        if match:
            tid_counts[match.group(1)] += 1

    duplicates = {tid: cnt for tid, cnt in tid_counts.items() if cnt > 1}
    if not duplicates:
        return False

    if logger:
        dup_total = sum(cnt - 1 for cnt in duplicates.values())
        logger.warning(
            f"GTF文件包含 {len(duplicates)} 个重复的transcript_id"
            f"|GTF file contains {len(duplicates)} duplicate transcript_ids: {os.path.basename(gtf_file)}"
        )
        logger.warning(
            f"共需修复 {dup_total} 处重复|Total {dup_total} duplicates to fix"
        )

    # 修复重复的transcript_id：为重复项添加序号后缀（全局去重）
    # TSEBRA要求GTF中所有transcript_id全局唯一
    # Fix duplicate transcript_ids: add sequential suffix (global dedup)
    # TSEBRA requires all transcript_ids in GTF to be globally unique
    tid_occurrence = defaultdict(int)
    fixed_lines = []
    # 当前transcript block的重命名映射
    current_tid_rename = None  # None=不重命名, str=新的transcript_id
    current_gid_rename = None

    for line in lines:
        if line.startswith('#') or not line.strip():
            fixed_lines.append(line)
            continue

        # 检测gene切换（新transcript block的开始）
        gid_match = re.search(r'gene_id "([^"]+)"', line)
        tid_match = re.search(r'transcript_id "([^"]+)"', line)

        if gid_match:
            gid = gid_match.group(1)
            if tid_match:
                tid = tid_match.group(1)
                tid_occurrence[tid] += 1

                if tid_occurrence[tid] > 1:
                    # 重复的transcript_id，添加序号后缀
                    # Duplicate transcript_id, add sequential suffix
                    suffix = tid_occurrence[tid]
                    current_tid_rename = f"{tid}.{suffix}"
                    current_gid_rename = re.sub(r'_g$', f'.{suffix}_g', gid)
                else:
                    current_tid_rename = None
                    current_gid_rename = None

        if tid_match and current_tid_rename:
            line = re.sub(r'transcript_id "[^"]+"', f'transcript_id "{current_tid_rename}"', line)
        if gid_match and current_gid_rename:
            line = re.sub(r'gene_id "[^"]+"', f'gene_id "{current_gid_rename}"', line)

        fixed_lines.append(line)

    # 写回文件|Write back to file
    import tempfile
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False, dir=os.path.dirname(gtf_file))
    try:
        for line in fixed_lines:
            temp_file.write(line)
        temp_file.close()
        os.replace(temp_file.name, gtf_file)
    except Exception:
        temp_file.close()
        if os.path.exists(temp_file.name):
            os.remove(temp_file.name)
        raise

    if logger:
        fixed_count = sum(1 for v in tid_occurrence.values() if v > 1)
        logger.info(
            f"GTF文件transcript_id去重完成|GTF transcript_id dedup completed: "
            f"{os.path.basename(gtf_file)}, 修复 {fixed_count} 个重复transcript"
            f"|fixed {fixed_count} duplicate transcripts"
        )

    return True


def fix_duplicate_protein_ids(fasta_file: str, logger=None) -> str:
    """
    检测并修复蛋白质FASTA文件中的重复序列ID
    Detect and fix duplicate sequence IDs in protein FASTA file

    Args:
        fasta_file: FASTA文件路径|Path to FASTA file
        logger: 日志器|Logger object

    Returns:
        str: 返回修复后的文件路径（如果无重复则返回原路径）|Return fixed file path (original if no duplicates)
    """
    from collections import defaultdict

    # 读取FASTA文件，记录每个ID出现的次数和位置
    # Read FASTA file, record occurrences and positions of each ID
    id_positions = []  # 存储所有ID及其位置|Store all IDs and their positions
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # 保存前一个序列
                # Save previous sequence
                if current_seq:
                    id_positions.append((seq_id, ''.join(current_seq)))
                    current_seq = []
                # 提取ID（去除>和可能的描述）
                # Extract ID (remove > and possible description)
                seq_id = line[1:].split()[0].strip()
            else:
                current_seq.append(line)

        # 保存最后一个序列
        # Save last sequence
        if current_seq:
            id_positions.append((seq_id, ''.join(current_seq)))

    # 检查重复ID
    # Check duplicate IDs
    id_counts = defaultdict(list)
    for idx, (seq_id, _) in enumerate(id_positions):
        id_counts[seq_id].append(idx)

    duplicates = {id_: positions for id_, positions in id_counts.items() if len(positions) > 1}

    if not duplicates:
        if logger:
            logger.debug("蛋白质文件无重复ID|No duplicate IDs in protein file")
        return fasta_file

    # 发现重复ID，需要修复
    # Found duplicate IDs, need to fix
    if logger:
        logger.warning(f"发现 {len(duplicates)} 个重复的蛋白质序列ID|Found {len(duplicates)} duplicate protein sequence IDs")
        for dup_id, positions in list(duplicates.items())[:5]:  # 只显示前5个|Show first 5 only
            logger.warning(f"  重复ID|Duplicate ID: {dup_id} (出现 {len(positions)} 次|occurs {len(positions)} times)")
        if len(duplicates) > 5:
            logger.warning(f"  ... 还有 {len(duplicates) - 5} 个重复ID|... and {len(duplicates) - 5} more duplicate IDs")

    # 创建修复后的文件
    # Create fixed file
    import tempfile
    temp_fixed = tempfile.NamedTemporaryFile(mode='w', suffix='_fixed_proteins.fa', delete=False, dir=os.path.dirname(fasta_file))
    fixed_path = temp_fixed.name

    # 用于跟踪每个ID的重复次数
    # Track duplicate count for each ID
    id_counters = defaultdict(int)

    with open(fixed_path, 'w') as out:
        for seq_id, seq in id_positions:
            id_counters[seq_id] += 1

            if id_counters[seq_id] > 1:
                # 重复ID，添加后缀使其唯一
                # Duplicate ID, add suffix to make it unique
                new_id = f"{seq_id}_dup{id_counters[seq_id]}"
                if logger and id_counters[seq_id] <= 5:  # 只记录前5个|Only log first 5
                    logger.debug(f"重命名序列ID|Renaming sequence ID: {seq_id} -> {new_id}")
                out.write(f">{new_id}\n{seq}")
            else:
                # 首次出现的ID，保持不变
                # First occurrence, keep unchanged
                out.write(f">{seq_id}\n{seq}")

    # 删除原文件，重命名修复后的文件
    # Remove original file, rename fixed file
    os.remove(fasta_file)
    os.rename(fixed_path, fasta_file)

    if logger:
        logger.info(f"蛋白质文件重复ID已修复|Duplicate protein IDs fixed: {fasta_file}")
        logger.info(f"共修复 {sum(len(positions) - 1 for positions in duplicates.values())} 个重复ID|Fixed {sum(len(positions) - 1 for positions in duplicates.values())} duplicate IDs in total")

    return fasta_file
