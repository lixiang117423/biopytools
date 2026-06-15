"""
Cactus工具函数|Cactus Utility Functions
"""

import os
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional


class CactusLogger:
    """Cactus日志类|Cactus Logger Class"""

    def __init__(self, log_file: Path, log_level: str = "INFO"):
        """
        初始化日志|Initialize logger

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR)
        """
        import sys

        self.log_file = Path(log_file)
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)

        # 确保日志目录存在|Ensure log directory exists
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 创建logger|Create logger
        self.logger = logging.getLogger("Cactus")
        self.logger.setLevel(logging.DEBUG)

        # 清除已有的处理器|Clear existing handlers
        self.logger.handlers = []
        self.logger.propagate = False

        # 格式化器|Formatter
        # 格式: YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 消息中文|Message English
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        # → 本地完整日志|→ Local complete log
        file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger对象|Get logger object"""
        return self.logger


def validate_singularity(singularity_path: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证Singularity是否可用|Validate Singularity availability

    Args:
        singularity_path: Singularity可执行文件路径|Singularity executable path
        logger: 日志对象|Logger object

    Returns:
        是否可用|Whether available
    """
    if logger:
        logger.info(" 检查Singularity|Checking Singularity")

    # 检查文件是否存在|Check if file exists
    singularity_path_expanded = Path(singularity_path).expanduser()
    if not singularity_path_expanded.exists():
        if logger:
            logger.error(f"   Singularity不存在|Singularity does not exist: {singularity_path_expanded}")
        return False

    # 测试Singularity是否可执行|Test if Singularity is executable
    try:
        result = subprocess.run(
            [str(singularity_path_expanded), "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            if logger:
                version_info = result.stdout.strip().split('\n')[0]
                logger.info(f"   Singularity可用|Singularity available: {version_info}")
            return True
        else:
            if logger:
                logger.error(f"   Singularity执行失败|Singularity execution failed: {result.stderr}")
            return False
    except Exception as e:
        if logger:
            logger.error(f"   Singularity检查失败|Singularity check failed: {e}")
        return False


def validate_cactus_sif(cactus_sif: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证Cactus SIF镜像是否存在|Validate Cactus SIF image existence

    Args:
        cactus_sif: Cactus SIF文件路径|Cactus SIF file path
        logger: 日志对象|Logger object

    Returns:
        是否存在|Whether exists
    """
    if logger:
        logger.info(" 检查Cactus SIF镜像|Checking Cactus SIF image")

    cactus_sif_expanded = Path(cactus_sif).expanduser()
    if not cactus_sif_expanded.exists():
        if logger:
            logger.error(f"   Cactus SIF不存在|Cactus SIF does not exist: {cactus_sif_expanded}")
        return False

    if logger:
        logger.info(f"   Cactus SIF找到|Cactus SIF found: {cactus_sif_expanded}")

    return True


def create_absolute_path_seqfile(seqfile: str, output_dir: str, logger: Optional[logging.Logger] = None) -> str:
    """
    创建使用绝对路径的seqfile|Create seqfile with absolute paths

    Args:
        seqfile: 原始seqfile路径|Original seqfile path
        output_dir: 输出目录|Output directory
        logger: 日志对象|Logger object

    Returns:
        新的seqfile路径（使用绝对路径）|New seqfile path (with absolute paths)
    """
    try:
        seqfile_path = Path(seqfile).expanduser()
        output_path = Path(output_dir)
        seqfile_dir = seqfile_path.parent

        # 读取原始seqfile|Read original seqfile
        with open(seqfile_path, 'r') as f:
            lines = f.readlines()

        # 创建新的seqfile文件名|Create new seqfile filename
        new_seqfile = output_path / f"{seqfile_path.stem}.abs{seqfile_path.suffix}"

        with open(new_seqfile, 'w') as f:
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sample_name = parts[0]
                    genome_path = ' '.join(parts[1:]) if len(parts) > 2 else parts[1]

                    # 转换为绝对路径|Convert to absolute path
                    path_obj = Path(genome_path)
                    if not path_obj.is_absolute():
                        path_obj = (seqfile_dir / genome_path).resolve()
                    else:
                        path_obj = path_obj.resolve()

                    f.write(f"{sample_name}\t{path_obj}\n")
                else:
                    f.write(line)

        if logger:
            logger.info(f" 创建绝对路径seqfile|Created absolute path seqfile: {new_seqfile}")

        return str(new_seqfile)

    except Exception as e:
        if logger:
            logger.error(f" 创建绝对路径seqfile失败|Failed to create absolute path seqfile: {e}")
        raise


def check_genome_files(seqfile: str, logger: Optional[logging.Logger] = None) -> Tuple[bool, list]:
    """
    检查序列文件中列出的所有基因组文件是否存在|Check all genome files listed in seqfile

    Args:
        seqfile: 序列文件路径|Sequence file path (两列格式：样本名 + 路径|Two columns: sample_name + path)
        logger: 日志对象|Logger object

    Returns:
        (是否全部存在, 基因组文件列表)|(Whether all exist, genome file list)
    """
    if logger:
        logger.info(" 检查基因组文件|Checking genome files")

    try:
        seqfile_path = Path(seqfile).expanduser()
        seqfile_dir = seqfile_path.parent

        with open(seqfile_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        # 解析两列格式：样本名 + 路径|Parse two-column format: sample_name + path
        genome_entries = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 2:
                sample_name = parts[0]
                # 路径可能是剩余部分拼接（如果路径中有空格）|Path may be rest of parts joined (if spaces in path)
                genome_path = ' '.join(parts[1:]) if len(parts) > 2 else parts[1]
                genome_entries.append((sample_name, genome_path))

        all_exist = True
        resolved_files = []

        for i, (sample_name, genome_file) in enumerate(genome_entries):
            # 尝试解析路径|Try to resolve path
            genome_path = Path(genome_file)

            # 如果是相对路径，尝试相对于seqfile目录|If relative path, try relative to seqfile dir
            if not genome_path.is_absolute():
                genome_path = seqfile_dir / genome_file

            # 展开用户目录|Expand user directory
            genome_path = genome_path.expanduser()

            if genome_path.exists():
                resolved_files.append(str(genome_path))
                if logger:
                    logger.info(f"   [{i+1}/{len(genome_entries)}] {sample_name}: 找到|Found: {genome_path}")
            else:
                all_exist = False
                if logger:
                    logger.warning(f"   [{i+1}/{len(genome_entries)}] {sample_name}: 未找到|Not found: {genome_file}")

        if all_exist:
            if logger:
                logger.info(f"   所有基因组文件都存在|All {len(genome_entries)} genome files exist")
        else:
            if logger:
                logger.error("   部分基因组文件缺失|Some genome files are missing")

        return all_exist, resolved_files

    except Exception as e:
        if logger:
            logger.error(f"   基因组文件检查失败|Genome file check failed: {e}")
        return False, []


def build_singularity_command(
    singularity_path: str,
    cactus_sif: str,
    jobstore: str,
    seqfile: str,
    output_dir: str,
    out_name: str,
    reference: str,
    output_formats: list,
    threads: int,
    max_memory: str,
    bind_args: list,
    work_dir: str,
    logger: Optional[logging.Logger] = None
) -> str:
    """
    构建Singularity调用cactus-pangenome的命令|Build Singularity command for cactus-pangenome

    Args:
        singularity_path: Singularity可执行文件路径|Singularity executable path
        cactus_sif: Cactus SIF镜像路径|Cactus SIF image path
        jobstore: Toil jobstore路径|Toil jobstore path
        seqfile: 序列文件路径|Sequence file path
        output_dir: 输出目录|Output directory
        out_name: 输出名称|Output name
        reference: 参考基因组名称|Reference genome name
        output_formats: 输出格式列表|Output format list
        threads: 线程数|Number of threads
        max_memory: 最大内存|Maximum memory
        bind_args: 绑定参数列表|Bind arguments list
        work_dir: 工作目录（用于临时文件）|Work directory (for temporary files)
        logger: 日志对象|Logger object

    Returns:
        完整的Singularity命令|Complete Singularity command
    """
    # 基础命令|Base command
    cmd_parts = []

    # 设置TMPDIR环境变量（修复cactus写入/tmp的bug）|Set TMPDIR env var (fix cactus /tmp bug)
    cmd_parts.append(f"TMPDIR={work_dir}")

    # Singularity执行|Singularity execution
    cmd_parts.append(f"{singularity_path} exec")

    # 绑定参数|Bind arguments
    if bind_args:
        cmd_parts.extend(bind_args)

    # Cactus SIF镜像|Cactus SIF image
    cmd_parts.append(f"{cactus_sif}")

    # cactus-pangenome命令|cactus-pangenome command
    cmd_parts.append("cactus-pangenome")

    # jobstore|jobstore
    cmd_parts.append(f"{jobstore}")

    # seqfile|seqfile
    cmd_parts.append(f"{seqfile}")

    # 输出目录|Output directory
    cmd_parts.append(f"--outDir {output_dir}")

    # 输出名称|Output name
    cmd_parts.append(f"--outName {out_name}")

    # 参考基因组|Reference genome
    cmd_parts.append(f"--reference {reference}")

    # 输出格式|Output formats
    # 注意：Cactus使用独立参数，不是--format|Note: Cactus uses separate flags, not --format
    # HAL是默认输出，不需要参数（且--hal在cactus中指输入HAL文件，不是输出格式）
    # HAL is default output, no flag needed (--hal means input HAL files in cactus, not output format)
    for fmt in output_formats:
        if fmt in ('psa', 'hal'):
            continue
        # 直接使用格式名作为参数|Use format name directly as flag
        cmd_parts.append(f"--{fmt}")

    # Toil参数|Toil parameters
    # 工作目录|Work directory (cactus会在此创建临时文件)
    cmd_parts.append(f"--workDir {work_dir}")

    # 工作线程|Worker threads (使用defaultCores而不是workers)
    cmd_parts.append(f"--defaultCores {threads}")

    # 最大内存|Maximum memory (Toil格式|Toil format: --defaultMemory)
    # 转换|Convert: 100G -> 100000000000
    memory_bytes = convert_memory_to_bytes(max_memory)
    cmd_parts.append(f"--defaultMemory {memory_bytes}")

    # 构建完整命令|Build complete command
    command = " ".join(cmd_parts)

    if logger:
        logger.debug(f"   构建的命令|Built command: {command}")

    return command


def convert_memory_to_bytes(memory_str: str) -> int:
    """
    将内存字符串转换为字节数|Convert memory string to bytes

    Args:
        memory_str: 内存字符串|Memory string (e.g., "100G", "50M")

    Returns:
        字节数|Number of bytes
    """
    memory_str = memory_str.upper().strip()

    if memory_str.endswith('G'):
        return int(memory_str[:-1]) * 1024 * 1024 * 1024
    elif memory_str.endswith('M'):
        return int(memory_str[:-1]) * 1024 * 1024
    elif memory_str.endswith('K'):
        return int(memory_str[:-1]) * 1024
    else:
        # 假设是字节|Assume bytes
        return int(memory_str)


def get_expected_output_files(out_name: str, output_formats: list) -> list:
    """
    根据输出格式获取cactus实际产生的文件名列表|Get actual output filenames produced by cactus based on formats

    cactus-pangenome (v3.1.4) 实际输出文件名规则|Actual output filename rules:
        - HAL: {out_name}.full.hal (不是|not {out_name}.hal)
        - GFA: {out_name}.gfa.gz (不是|not {out_name}.gfa)
        - GBZ: {out_name}.gbz
        - ODGI: {out_name}.full.og (不是|not {out_name}.odgi)

    Args:
        out_name: 输出名称|Output name
        output_formats: 输出格式列表|Output format list

    Returns:
        期望的文件名列表|List of expected filenames
    """
    # cactus实际输出文件名映射|Actual cactus output filename mapping
    format_to_filename = {
        "hal": f"{out_name}.full.hal",
        "gfa": f"{out_name}.gfa.gz",
        "gbz": f"{out_name}.gbz",
        "odgi": f"{out_name}.full.og",
    }

    expected_files = []
    for fmt in output_formats:
        filename = format_to_filename.get(fmt)
        if filename:
            expected_files.append(filename)

    return expected_files


def check_output_files(output_dir: str, out_name: str, output_formats: list) -> bool:
    """
    检查输出文件是否已存在|Check if output files already exist

    Args:
        output_dir: 输出目录|Output directory
        out_name: 输出名称|Output name
        output_formats: 输出格式列表|Output format list

    Returns:
        是否所有输出文件都存在|Whether all output files exist
    """
    output_path = Path(output_dir)
    expected_files = get_expected_output_files(out_name, output_formats)

    # 检查文件是否存在|Check if files exist
    for expected_file in expected_files:
        file_path = output_path / expected_file
        if not file_path.exists():
            return False

    return True


def get_genome_count(seqfile: str) -> int:
    """
    获取序列文件中的基因组数量|Get number of genomes in seqfile

    Args:
        seqfile: 序列文件路径|Sequence file path

    Returns:
        基因组数量|Number of genomes
    """
    try:
        with open(seqfile, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        return len(lines)
    except Exception:
        return 0


def parse_cactus_version(cactus_sif: str, singularity_path: str, logger: Optional[logging.Logger] = None) -> Optional[str]:
    """
    解析Cactus版本信息|Parse Cactus version information

    Args:
        cactus_sif: Cactus SIF镜像路径|Cactus SIF image path
        singularity_path: Singularity可执行文件路径|Singularity executable path
        logger: 日志对象|Logger object

    Returns:
        版本信息|Version information
    """
    try:
        result = subprocess.run(
            [singularity_path, "exec", cactus_sif, "cactus", "--version"],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            version = result.stdout.strip()
            if logger:
                logger.info(f"   Cactus版本|Cactus version: {version}")
            return version
        else:
            if logger:
                logger.warning("   无法获取Cactus版本|Cannot get Cactus version")
            return None

    except Exception as e:
        if logger:
            logger.warning(f"   Cactus版本解析失败|Cactus version parsing failed: {e}")
        return None
