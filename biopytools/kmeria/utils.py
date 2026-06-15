"""KMERIA工具函数模块|KMERIA Utility Functions Module"""

import logging
import sys
import os
import subprocess
import gzip
import shutil
import re
from pathlib import Path
from typing import List, Optional, Tuple
from datetime import datetime


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
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
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list

    Note:
        必须使用--no-capture-output避免conda缓冲输出导致内存溢出|Must use --no-capture-output to avoid OOM
        传递完整路径以正确识别conda环境，禁止使用os.path.basename()提取命令名|Pass full path to detect conda env, never use os.path.basename()
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class KMERIALogger:
    """KMERIA日志管理器|KMERIA Logger Manager

    按开发规范§2.3实现stdout/stderr/file三handler分离|Implements stdout/stderr/file separation per spec §2.3
    INFO→stdout→.out，WARNING+→stderr→.err，全部→本地日志文件|INFO→stdout→.out, WARNING+→stderr→.err, all→local log
    """

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.log_level = log_level
        self.logger = self._build_logger()

    def _build_logger(self):
        """构建日志器（避免basicConfig不可重入问题）|Build logger (avoids basicConfig non-reentrant issue)"""
        logger = logging.getLogger("biopytools.kmeria")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # stdout handler - INFO级别（被超算捕获到.out文件）|stdout handler - INFO level (captured to .out)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上（被超算捕获到.err文件）|stderr handler - WARNING+ (captured to .err)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        return logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, output_dir=None):
        self.logger = logger
        self.output_dir = output_dir
        self.commands_log = []

    def run_command(self, cmd: List[str], description: str = "",
                   check=True, capture_output=False, output_file=None) -> Tuple[bool, str]:
        """
        执行命令（自动检测conda环境）|Execute command (auto-detect conda environment)

        Args:
            cmd: 命令列表（首元素必须为完整路径）|Command list (first element must be full path)
            description: 命令描述（已包含中英文对照）|Command description (already bilingual)
            check: 是否检查返回码|Whether to check return code
            capture_output: 是否捕获输出|Whether to capture output
            output_file: 输出文件路径（将stdout重定向到文件）|Output file path (redirect stdout to file)

        Returns:
            (success, output): (是否成功|Success, 输出|Output)

        Note:
            传递完整路径给build_conda_command以正确识别conda环境|Pass full path to detect conda env
            命令完整内容记录在INFO级别便于追溯|Full command logged at INFO level for reproducibility
        """
        # 自动包装conda环境的命令（必须传完整路径）|Auto-wrap with conda env (must use full path)
        if cmd:
            wrapped_cmd = build_conda_command(cmd[0], cmd[1:])
        else:
            wrapped_cmd = cmd

        cmd_str = ' '.join(wrapped_cmd)

        if description:
            self.logger.info(description)
        self.logger.info(f"命令|Command: {cmd_str}")

        # 记录命令|Log command
        self.commands_log.append({
            'timestamp': datetime.now().isoformat(),
            'command': cmd_str,
            'description': description
        })

        try:
            # 处理输出重定向|Handle output redirection
            if output_file:
                # 打开输出文件|Open output file
                with open(output_file, 'w') as f:
                    result = subprocess.run(
                        wrapped_cmd,
                        stdout=f,
                        stderr=subprocess.PIPE,
                        check=check,
                        text=True
                    )
                return True, ""
            elif capture_output:
                result = subprocess.run(
                    wrapped_cmd,
                    check=check,
                    capture_output=True,
                    text=True
                )
                output = result.stdout
                if result.stderr:
                    output += result.stderr
                return True, output
            else:
                subprocess.run(wrapped_cmd, check=check)
                return True, ""

        except subprocess.CalledProcessError as e:
            error_msg = f"命令执行失败|Command execution failed: {e}"
            self.logger.error(error_msg)
            if capture_output and e.stderr:
                error_msg += f"\n{e.stderr}"
            return False, error_msg

        except FileNotFoundError:
            error_msg = f"命令未找到|Command not found: {wrapped_cmd[0]}"
            self.logger.error(error_msg)
            return False, error_msg

        except Exception as e:
            error_msg = f"未知错误|Unknown error: {e}"
            self.logger.error(error_msg)
            return False, error_msg

    def save_commands_log(self, output_file):
        """保存命令日志|Save commands log"""
        import json

        with open(output_file, 'w') as f:
            json.dump(self.commands_log, f, indent=2)

        self.logger.info(f"命令日志已保存|Commands log saved: {output_file}")


def format_number(num: float) -> str:
    """
    格式化数字|Format number

    大数字使用M(百万)单位显示|Large numbers displayed in M(million) units

    Args:
        num: 数字|Number

    Returns:
        格式化后的字符串|Formatted string
    """
    if isinstance(num, (int, float)):
        if num >= 1_000_000:
            return f"{num / 1_000_000:.2f}M"
        elif num >= 1_000:
            return f"{num / 1_000:.2f}K"
        else:
            return f"{num:.2f}"
    return str(num)


def read_sample_list(samples_file: str) -> List[str]:
    """
    读取样本列表|Read sample list

    Args:
        samples_file: 样本列表文件|Sample list file

    Returns:
        样本名列表|List of sample names
    """
    samples = []
    with open(samples_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                samples.append(line)

    return samples


def find_fastq_files(fastq_dir: str, samples: List[str],
                     pattern: str = '*.fq.gz') -> dict:
    """
    查找FASTQ文件|Find FASTQ files

    Args:
        fastq_dir: FASTQ文件目录|FASTQ directory
        samples: 样本列表|Sample list
        pattern: 文件模式|File pattern

    Returns:
        样本到FASTQ文件的映射|Sample to FASTQ files mapping
    """
    from glob import glob

    sample_files = {}

    for sample in samples:
        # 尝试多种命名模式|Try multiple naming patterns
        patterns = [
            os.path.join(fastq_dir, f"{sample}*.{pattern.replace('*', '')}"),
            os.path.join(fastq_dir, f"{sample}_R1*"),
            os.path.join(fastq_dir, f"{sample}_1*"),
            os.path.join(fastq_dir, f"{sample}*")
        ]

        files = []
        for p in patterns:
            found = glob(p)
            if found:
                files.extend(found)
                break

        if files:
            # 排序以确保R1在R2之前|Sort to ensure R1 comes before R2
            files.sort()
            sample_files[sample] = files
        else:
            print(f"警告|Warning: 未找到样本文件|Sample files not found: {sample}")

    return sample_files


def create_sample_list_from_fastq(fastq_dir: str,
                                   pattern: str = '*.fq.gz',
                                   output_file: str = 'samples.txt') -> List[str]:
    """
    从FASTQ文件创建样本列表|Create sample list from FASTQ files

    Args:
        fastq_dir: FASTQ文件目录|FASTQ directory
        pattern: 文件模式|File pattern
        output_file: 输出文件名|Output filename

    Returns:
        样本列表|Sample list
    """
    from glob import glob
    import re

    files = glob(os.path.join(fastq_dir, pattern))
    samples = set()

    for f in files:
        # 提取样本名|Extract sample name
        basename = os.path.basename(f)
        # 移除常见的后缀|Remove common suffixes
        sample = re.sub(r'[_\.][Rr]?[12]\.f[aq](\.gz)?$', '', basename)
        sample = re.sub(r'\.f[aq](\.gz)?$', '', sample)
        samples.add(sample)

    samples = sorted(list(samples))

    # 保存到文件|Save to file
    with open(output_file, 'w') as f:
        for sample in samples:
            f.write(f"{sample}\n")

    return samples


def calculate_sequencing_depth(bam_file: str) -> float:
    """
    计算测序深度|Calculate sequencing depth

    Args:
        bam_file: BAM文件|BAM file

    Returns:
        平均深度|Average depth
    """
    samtools_path = shutil.which('samtools')
    if not samtools_path:
        print("samtools未找到|samtools not found in PATH")
        return 0.0

    wrapped_cmd = build_conda_command(samtools_path, ['depth', bam_file])

    try:
        result = subprocess.run(wrapped_cmd, capture_output=True, text=True, check=True)

        total_depth = 0
        count = 0

        for line in result.stdout.split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    total_depth += int(parts[2])
                    count += 1

        if count > 0:
            return total_depth / count

    except Exception as e:
        print(f"计算深度失败|Failed to calculate depth: {e}")

    return 0.0


def parse_depth_file(depth_file: str) -> dict:
    """
    解析深度文件|Parse depth file

    格式|Format: sample_name depth [ploidy]

    Args:
        depth_file: 深度文件|Depth file

    Returns:
        样本到深度(和倍性)的映射|Sample to depth (and ploidy) mapping
    """
    depth_dict = {}

    with open(depth_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                sample = parts[0]
                depth = float(parts[1])
                ploidy = int(parts[2]) if len(parts) > 2 else None

                depth_dict[sample] = {
                    'depth': depth,
                    'ploidy': ploidy
                }

    return depth_dict


def check_file_exists(file_path: str, description: str = "文件|File") -> bool:
    """
    检查文件是否存在|Check if file exists

    Args:
        file_path: 文件路径|File path
        description: 文件描述|File description

    Returns:
        是否存在|Whether exists
    """
    if not os.path.exists(file_path):
        print(f"{description}不存在|{description} does not exist: {file_path}")
        return False
    return True


def check_kmeria_installation(kmeria_path: str) -> bool:
    """
    检查KMERIA安装|Check KMERIA installation

    Args:
        kmeria_path: KMERIA路径|KMERIA path

    Returns:
        是否正确安装|Whether correctly installed
    """
    kmeria_bin = os.path.join(kmeria_path, 'bin', 'kmeria')

    if not os.path.exists(kmeria_bin):
        print(f"KMERIA可执行文件不存在|KMERIA executable not found: {kmeria_bin}")
        return False

    # 测试运行（必须传完整路径以正确识别conda环境）|Test run (must use full path to detect conda env)
    wrapped_cmd = build_conda_command(kmeria_bin, [])

    try:
        result = subprocess.run(wrapped_cmd, capture_output=True, text=True, timeout=10)
        # kmeria不带参数会显示帮助|kmeria shows help without arguments
        return 'KMERIA' in result.stdout or 'kmeria' in result.stdout.lower()
    except Exception as e:
        print(f"KMERIA测试失败|KMERIA test failed: {e}")
        return False


def get_kmeria_command(kmeria_path: str, subcommand: str) -> List[str]:
    """
    获取kmeria命令|Get kmeria command

    Args:
        kmeria_path: KMERIA路径|KMERIA path
        subcommand: 子命令|Subcommand

    Returns:
        完整命令列表|Full command list
    """
    kmeria_bin = os.path.join(kmeria_path, 'bin', 'kmeria')
    return [kmeria_bin, subcommand]


def validate_phenotype_file(pheno_file: str, samples: List[str],
                             pheno_col: int = 1) -> Tuple[bool, List[str]]:
    """
    验证表型文件|Validate phenotype file

    Args:
        pheno_file: 表型文件|Phenotype file
        samples: 样本列表|Sample list
        pheno_col: 表型列(1-based)|Phenotype column (1-based)

    Returns:
        (是否有效|Valid, 错误信息|Error messages)
    """
    errors = []
    sample_set = set(samples)

    try:
        with open(pheno_file, 'r') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < pheno_col + 1:
                    errors.append(f"第{i}行：列数不足|Line {i}: Not enough columns")
                    continue

                sample = parts[0]
                if sample not in sample_set:
                    errors.append(f"第{i}行：样本{sample}不在样本列表中|Line {i}: Sample {sample} not in sample list")

                try:
                    phenotype = float(parts[pheno_col])
                except ValueError:
                    errors.append(f"第{i}行：表型值不是数字|Line {i}: Phenotype value is not numeric")

    except Exception as e:
        errors.append(f"文件读取错误|File reading error: {e}")

    return len(errors) == 0, errors


def generate_manhattan_data(asso_file: str, output_file: str):
    """
    生成曼哈顿图数据|Generate Manhattan plot data

    Args:
        asso_file: 关联分析结果文件|Association result file
        output_file: 输出文件|Output file
    """
    import pandas as pd
    import numpy as np

    try:
        # 读取关联结果|Read association results
        df = pd.read_csv(asso_file, sep='\t', header=0)

        # 提取需要的列|Extract needed columns
        # 假设格式：kmer, p_value, chr, pos
        if 'p_value' in df.columns or 'pval' in df.columns or 'P' in df.columns:
            pval_col = 'p_value' if 'p_value' in df.columns else ('pval' if 'pval' in df.columns else 'P')

            # 计算-log10(p)|Calculate -log10(p)
            df['neglog10p'] = -np.log10(df[pval_col])

            # 保存|Save
            df.to_csv(output_file, sep='\t', index=False)
            return True

    except Exception as e:
        print(f"生成曼哈顿图数据失败|Failed to generate Manhattan plot data: {e}")

    return False


def compress_file(input_file: str, output_file: str = None) -> str:
    """
    压缩文件|Compress file

    Args:
        input_file: 输入文件|Input file
        output_file: 输出文件|Output file

    Returns:
        输出文件路径|Output file path
    """
    if output_file is None:
        output_file = input_file + '.gz'

    with open(input_file, 'rb') as f_in:
        with gzip.open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)

    return output_file


def decompress_file(input_file: str, output_file: str = None) -> str:
    """
    解压文件|Decompress file

    Args:
        input_file: 输入文件|Input file
        output_file: 输出文件|Output file

    Returns:
        输出文件路径|Output file path
    """
    if output_file is None:
        output_file = input_file.rstrip('.gz')

    with gzip.open(input_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)

    return output_file
