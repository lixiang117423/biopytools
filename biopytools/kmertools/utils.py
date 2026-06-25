"""
K-mer工具工具函数模块|K-mer Tools Utility Functions Module
"""

import logging
import sys
import subprocess
import glob
import os
import shutil
import re
from pathlib import Path
from typing import List, Tuple, Optional


class KmerToolsLogger:
    """K-mer工具日志管理器|K-mer Tools Logger Manager"""

    def __init__(self, log_file: Path, log_name: str = "kmertools.log"):
        self.log_file = log_file
        self.log_name = log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir=None):
        self.logger = logger
        self.working_dir = working_dir or os.getcwd()

    def check_command_exists(self, command: str) -> bool:
        """检查命令是否存在|Check if command exists"""
        return shutil.which(command) is not None

    def run_command(self, cmd, description: str = "") -> bool:
        """执行命令（列表形式）|Execute command (list format)"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 如果是字符串，转换为列表
        if isinstance(cmd, str):
            cmd_list = cmd.split()
        else:
            cmd_list = cmd

        self.logger.info(f"命令|Command: {' '.join(cmd_list)}")

        try:
            result = subprocess.run(
                cmd_list,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"标准错误|Stderr: {result.stderr}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {e}")
            return False

    def run(self, cmd: str, description: str = "", check: bool = True) -> bool:
        """执行命令（字符串形式）|Execute command (string format)"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=check,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"标准错误|Stderr: {result.stderr}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            if check:
                raise
            return False


class KmtricksChecker:
    """kmtricks工具检查器|kmtricks Tool Checker"""

    def __init__(self, logger, kmtricks_path: str = 'kmtricks'):
        self.logger = logger
        self.kmtricks_path = kmtricks_path

    def check_kmtricks(self):
        """检查kmtricks是否已安装|Check if kmtricks is installed"""
        try:
            result = subprocess.run(
                [self.kmtricks_path, 'infos'],
                capture_output=True,
                text=True,
                check=True
            )
            # 从输出中提取版本信息|Extract version from output
            version = "unknown"
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if 'kmtricks' in line.lower():
                        version = line.strip()
                        break
            self.logger.info(f"找到kmtricks|Found kmtricks: {version}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"未找到kmtricks|kmtricks not found: {e}")
            self.logger.error("请先安装kmtricks|Please install kmtricks first")
            self.logger.error("安装方法|Installation: conda install -c bioconda kmtricks")
            return False


class BgzipChecker:
    """bgzip工具检查器|bgzip Tool Checker"""

    def __init__(self, logger, bgzip_path: str = 'bgzip'):
        self.logger = logger
        self.bgzip_path = bgzip_path

    def check_bgzip(self):
        """检查bgzip是否已安装|Check if bgzip is installed"""
        try:
            result = subprocess.run(
                [self.bgzip_path, '--version'],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"找到bgzip|Found bgzip")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"未找到bgzip|bgzip not found: {e}")
            self.logger.error("请先安装bgzip|Please install bgzip first")
            self.logger.error("安装方法|Installation: conda install -c bioconda htslib")
            return False


class RocksDBChecker:
    """RocksDB Python库检查器|RocksDB Python Library Checker"""

    def __init__(self, logger):
        self.logger = logger

    def check_rocksdb(self):
        """检查python-rocksdb是否已安装|Check if python-rocksdb is installed"""
        try:
            import rocksdb
            self.logger.info(f"找到python-rocksdb|Found python-rocksdb")
            return True
        except ImportError:
            self.logger.error(f"未找到python-rocksdb|python-rocksdb not found")
            self.logger.error("请先安装python-rocksdb|Please install python-rocksdb first")
            self.logger.error("安装方法|Installation: pip install python-rocksdb")
            return False


class KmindexChecker:
    """kmindex工具检查器|kmindex Tool Checker"""

    def __init__(self, logger, kmindex_path: str = 'kmindex'):
        self.logger = logger
        self.kmindex_path = kmindex_path

    def check_kmindex(self):
        """检查kmindex是否已安装|Check if kmindex is installed"""
        try:
            result = subprocess.run(
                [self.kmindex_path, '--version'],
                capture_output=True,
                text=True,
                check=True
            )
            # 从输出中提取版本信息|Extract version from output
            version = "unknown"
            if result.stdout:
                version = result.stdout.strip()
            elif result.stderr:
                version = result.stderr.strip()
            self.logger.info(f"找到kmindex|Found kmindex: {version}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"未找到kmindex|kmindex not found: {e}")
            self.logger.error("请先安装kmindex|Please install kmindex first")
            self.logger.error("安装方法|Installation: conda install -c bioconda kmindex")
            return False


def generate_fof_file(input_dir: Path, output_file: Path,
                      suffix_1: str = "_1.clean.fq.gz",
                      suffix_2: str = "_2.clean.fq.gz",
                      logger: logging.Logger = None) -> bool:
    """
    生成FOF文件|Generate FOF file

    Args:
        input_dir: 输入目录|Input directory
        output_file: 输出FOF文件|Output FOF file
        suffix_1: R1文件后缀|R1 file suffix
        suffix_2: R2文件后缀|R2 file suffix
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Success
    """
    from collections import defaultdict
    import re

    if logger:
        logger.info(f"正在生成FOF文件|Generating FOF file: {output_file}")

    # 关联字典：sample -> list of full paths
    samples = defaultdict(list)

    # 遍历所有匹配后缀的文件
    pattern1 = str(input_dir / f"*{suffix_1}")
    pattern2 = str(input_dir / f"*{suffix_2}")

    for pattern in [pattern1, pattern2]:
        for f in glob.glob(pattern):
            if os.path.isfile(f):
                base = os.path.basename(f)
                # 去掉对应后缀得到 sample
                if suffix_1 in base:
                    sample = base.replace(suffix_1, "")
                elif suffix_2 in base:
                    sample = base.replace(suffix_2, "")
                else:
                    sample = base
                samples[sample].append(f)

    # 生成输出文件
    try:
        with open(output_file, "w") as out:
            for sample, files in sorted(samples.items()):
                f1 = None
                f2 = None
                for file_path in files:
                    base_file = os.path.basename(file_path)
                    if base_file.endswith(suffix_1):
                        f1 = file_path
                    elif base_file.endswith(suffix_2):
                        f2 = file_path

                if f1 and f2:
                    # 解析软链接到真实文件路径，并转为绝对路径
                    f1_resolved = _resolve_symlink_to_real_path(Path(f1))
                    f2_resolved = _resolve_symlink_to_real_path(Path(f2))

                    # 清理样本名：移除所有特殊字符（包括下划线），只保留字母数字
                    # kmtricks要求样本名不能有任何符号
                    clean_sample = re.sub(r'[^a-zA-Z0-9]', '', sample)

                    # kmtricks期望格式：sample: /path/R1 ; /path/R2
                    # 样本名后用冒号和空格，两个路径间用空格分号空格
                    out.write(f"{clean_sample}: {f1_resolved} ; {f2_resolved}\n")
                else:
                    warning_msg = f"Warning: Sample {sample} missing {suffix_1} or {suffix_2} files"
                    if logger:
                        logger.warning(warning_msg)
                    else:
                        print(warning_msg, file=sys.stderr)

        if logger:
            logger.info(f"FOF文件生成成功|FOF file generated successfully: {output_file}")
            logger.info(f"共|Total {len(samples)} 个样本|samples")
        return True

    except Exception as e:
        error_msg = f"生成FOF文件失败|Failed to generate FOF file: {e}"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg, file=sys.stderr)
        return False


def _resolve_symlink_to_real_path(file_path: Path) -> str:
    """
    解析软链接到真实文件路径，并返回绝对路径|Resolve symlink to real file path and return absolute path

    Args:
        file_path: 文件路径（可能是软链接）|File path (may be symlink)

    Returns:
        str: 解析后的绝对路径|Resolved absolute path
    """
    if not file_path.is_symlink():
        return str(file_path.absolute())

    try:
        # 读取软链接目标
        link_target = os.readlink(file_path)

        # 如果是绝对路径，直接使用
        if os.path.isabs(link_target):
            resolved = Path(link_target)
            if resolved.exists():
                return str(resolved)
            else:
                # 如果目标不存在，返回原路径的绝对路径
                return str(file_path.absolute())

        # 相对路径，基于软链接所在目录解析
        symlink_dir = file_path.parent
        resolved = (symlink_dir / link_target).resolve()

        if resolved.exists():
            return str(resolved)
        else:
            # 如果解析后的路径不存在，尝试使用realpath命令
            try:
                result = subprocess.run(
                    ['realpath', str(file_path)],
                    capture_output=True,
                    text=True,
                    check=False
                )
                if result.returncode == 0:
                    realpath = Path(result.stdout.strip())
                    if realpath.exists():
                        return str(realpath)
            except Exception:
                pass

            # 最后尝试返回原路径的绝对路径
            return str(file_path.absolute())

    except Exception as e:
        # 解析失败，返回绝对路径
        return str(file_path.absolute())


def split_fasta_file(input_fasta: Path, output_dir: Path,
                     logger: logging.Logger = None) -> int:
    """
    分割FASTA文件|Split FASTA file

    Args:
        input_fasta: 输入FASTA文件|Input FASTA file
        output_dir: 输出目录|Output directory
        logger: 日志器|Logger

    Returns:
        int: 分割的序列数量|Number of sequences split
    """
    try:
        from Bio import SeqIO
    except ImportError:
        error_msg = "未安装Biopython|Biopython not installed. 请先安装|Please install: pip install biopython"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg, file=sys.stderr)
        return 0

    if logger:
        logger.info(f"正在分割FASTA文件|Splitting FASTA file: {input_fasta}")

    sequence_count = 0

    try:
        with open(input_fasta, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                sequence_count += 1

                # 使用 record.id 作为文件名
                output_filename = f"{record.id}.fa"
                output_filepath = output_dir / output_filename

                # 打印处理进度
                if logger and sequence_count % 100 == 0:
                    logger.debug(f"已处理序列|Processed sequences: {sequence_count}")

                # 将单个序列记录写入其对应的文件
                try:
                    with open(output_filepath, "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                except IOError as e:
                    error_msg = f"写入文件失败|Failed to write file {output_filepath}: {e}"
                    if logger:
                        logger.error(error_msg)
                    else:
                        print(error_msg, file=sys.stderr)

        if logger:
            logger.info(f"成功将|Successfully split {sequence_count} 个序列分割到目录|sequences to directory '{output_dir}'")

        return sequence_count

    except FileNotFoundError:
        error_msg = f"错误|Error: 未找到输入文件|Input file not found: {input_fasta}"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg, file=sys.stderr)
        return 0
    except Exception as e:
        error_msg = f"发生未知错误|Unknown error occurred: {e}"
        if logger:
            logger.error(error_msg)
        else:
            print(error_msg, file=sys.stderr)
        return 0


def run_command(cmd: List[str], logger: logging.Logger = None,
               check: bool = True, capture_output: bool = True) -> subprocess.CompletedProcess:
    """
    运行命令|Run command

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code
        capture_output: 是否捕获输出|Whether to capture output

    Returns:
        subprocess.CompletedProcess: 命令执行结果|Command execution result
    """
    if logger:
        logger.debug(f"执行命令|Running command: {' '.join(cmd)}")

    return subprocess.run(
        cmd,
        capture_output=capture_output,
        text=True,
        check=check
    )


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称
    Detect if command is in conda environment, return environment name

    策略|Strategy:
    1. 首先尝试从which命令路径检测（优先级高）
    2. 如果未找到，搜索所有conda环境（兜底方案）

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'kmtricks' or '/path/to/kmtricks')

    Returns:
        conda环境名称或None|conda environment name or None (e.g., 'kmtricks_v.1.5.1' or None)
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs
        # 例如: /miniforge3/envs/kmtricks_v.1.5.1/bin/kmtricks
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件
    Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表 (适用于subprocess.run(shell=False))
        Complete command list (for subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('kmtricks', ['--version'])
        ['conda', 'run', '-n', 'kmtricks_v.1.5.1', '--no-capture-output', 'kmtricks', '--version']

        >>> # 绝对路径且不在conda envs目录下时，直接调用
        >>> # Absolute path not under conda envs: called directly
        >>> build_conda_command('/usr/bin/tool', ['--help'])
        ['/usr/bin/tool', '--help']

    注意|Note:
        返回的列表应配合 subprocess.run(shell=False) 使用
        The returned list must be used with subprocess.run(shell=False)
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用
        # 如果command是命令名，conda run会自动找到环境中的版本
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用
        full_cmd = [command] + args

    return full_cmd
