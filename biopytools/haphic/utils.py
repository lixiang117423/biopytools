"""
HapHiC工具函数模块|HapHiC Utility Functions Module
"""

import os
import shutil
import re
import subprocess
import sys
import time
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any
import tempfile


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    策略|Strategy:
    1. 如果command是完整路径且包含/envs/，直接从路径提取环境名（优先级最高）
    2. 如果是命令名，使用which命令路径检测
    3. 如果未找到，搜索所有conda环境（兜底方案）

    Args:
        command: 命令名称或路径 (e.g., 'busco' or '/path/to/busco')

    Returns:
        conda环境名称或None (e.g., 'BUSCO_v.6.0.0' or None)
    """
    # 方法0: 如果是完整路径且包含envs，直接提取（优先级最高）
    # Method 0: If full path contains envs, extract directly (highest priority)
    if '/envs/' in command:
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 解析符号链接的真实路径
        # Resolve symbolic links to real path
        if os.path.islink(cmd_path):
            real_path = os.path.realpath(cmd_path)
            cmd_path = real_path

        # 检查路径中是否包含 envs
        # 例如: /miniforge3/envs/BUSCO_v.6.0.0/bin/busco
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
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
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表 (适用于subprocess.run(shell=False))|Complete command list (for subprocess.run(shell=False))

    Examples:
        >>> build_conda_command('busco', ['--version'])
        ['conda', 'run', '-n', 'BUSCO_v.6.0.0', '--no-capture-output', 'busco', '--version']

        >>> build_conda_command('/miniforge3/envs/BUSCO_v.6.0.0/bin/busco', ['--version'])
        ['conda', 'run', '-n', 'BUSCO_v.6.0.0', '--no-capture-output', 'busco', '--version']

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
        # 使用conda run调用（必须使用--no-capture-output避免内存溢出）
        # Use conda run invocation (must use --no-capture-output to avoid memory overflow)
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用
        full_cmd = [command] + args

    return full_cmd


class HapHiCLogger:
    """HapHiC日志管理器|HapHiC Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, verbose: bool = False):
        self.log_file = log_file
        self.verbose = verbose
        self.logger = None

        self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器|Setup logger"""
        self.logger = logging.getLogger("haphic")
        self.logger.setLevel(logging.DEBUG)

        # 清除现有处理器
        self.logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter('%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S')

        # 控制台处理器 - 输出到stdout|Console handler - output to stdout
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO if not self.verbose else logging.DEBUG)
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件处理器
        if self.log_file:
            try:
                log_dir = os.path.dirname(os.path.abspath(self.log_file))
                Path(log_dir).mkdir(parents=True, exist_ok=True)

                file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
                file_handler.setLevel(logging.DEBUG)
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)
            except Exception as e:
                # 使用stderr作为fallback来显示错误信息
                print(f"无法创建日志文件|Cannot create log file: {e}", file=sys.stderr)

        # 确保日志不重复
        self.logger.propagate = False

    def get_logger(self):
        """获取日志记录器|Get logger"""
        return self.logger


class CommandRunner:
    """命令运行器|Command Runner"""

    def __init__(self, logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run

    def run_command(self, cmd: List[str], description: str = "") -> bool:
        """
        运行命令（自动检测conda环境）|Run command (auto-detect conda environment)

        Args:
            cmd: 命令列表|Command list
            description: 命令描述|Command description

        Returns:
            是否成功|Success status
        """
        if description:
            self.logger.info(f"{description}")

        # 自动包装conda环境的命令|Auto-wrap conda environment commands
        # 传递完整路径，让 build_conda_command 能够正确检测conda环境
        # Pass full path, so build_conda_command can detect conda environment correctly
        if cmd:
            wrapped_cmd = build_conda_command(cmd[0], cmd[1:])
        else:
            wrapped_cmd = cmd

        cmd_str = " ".join(wrapped_cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

        if self.dry_run:
            self.logger.info(f"[DRY RUN] {cmd_str}")
            return True

        try:
            # 设置环境变量
            env = os.environ.copy()
            env["LC_ALL"] = "C.UTF-8"
            env["LANG"] = "C.UTF-8"

            # 运行命令
            start_time = time.time()
            result = subprocess.run(
                wrapped_cmd,
                capture_output=True,
                text=True,
                env=env,
                timeout=None,  # 无时间限制
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )
            elapsed_time = time.time() - start_time

            # 统一处理stdout和stderr，很多工具的标准输出实际是stderr
            all_output = []
            if result.stdout:
                all_output.append(result.stdout)
            if result.stderr:
                all_output.append(result.stderr)

            combined_output = "\n".join(all_output).strip() if all_output else ""

            if result.returncode == 0:
                self.logger.info(f"命令执行成功|Command executed successfully (耗时: {elapsed_time:.2f}秒)")
                # 记录所有输出信息
                if combined_output:
                    # 检查是否是重要信息
                    if any(keyword in combined_output for keyword in ["Progress", "Completed", "Warning", "ERROR", "error"]):
                        self.logger.info(f"命令输出|Command Output:\n{combined_output}")
                    else:
                        self.logger.debug(f"命令输出|Command Output:\n{combined_output}")
                return True
            else:
                self.logger.error(f"命令执行失败|Command failed (返回码: {result.returncode})")
                # 记录所有输出信息
                if combined_output:
                    self.logger.error(f"命令输出|Command Output:\n{combined_output}")
                return False

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command timeout after 1 hour")
            return False

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {e}")
            return False

    def check_tool(self, tool_path: str) -> bool:
        """检查工具是否可用|Check if tool is available"""
        try:
            result = subprocess.run(
                [tool_path, "--help"],
                capture_output=True,
                text=True,
                timeout=10,
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )
            if result.returncode == 0:
                self.logger.debug(f"工具检查通过|Tool check passed: {tool_path}")
            else:
                self.logger.error(f"工具检查失败|Tool check failed: {tool_path} (返回码: {result.returncode})")
                if result.stderr:
                    self.logger.debug(f"工具错误信息|Tool error: {result.stderr}")
            return result.returncode == 0
        except Exception as e:
            self.logger.error(f"工具检查异常|Tool check exception: {tool_path} - {e}")
            return False


class FileManager:
    """文件管理器|File Manager"""

    @staticmethod
    def ensure_directory(dir_path: str) -> bool:
        """确保目录存在|Ensure directory exists"""
        try:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            return True
        except Exception as e:
            logging.error(f"无法创建目录|Cannot create directory: {dir_path}")
            logging.error(f"错误|Error: {e}")
            return False

    @staticmethod
    def check_file_exists(file_path: str) -> bool:
        """检查文件是否存在|Check if file exists"""
        return os.path.isfile(file_path)

    @staticmethod
    def get_file_size(file_path: str) -> int:
        """获取文件大小|Get file size"""
        try:
            return os.path.getsize(file_path)
        except OSError:
            return 0

    @staticmethod
    def get_file_extension(file_path: str) -> str:
        """获取文件扩展名|Get file extension"""
        return os.path.splitext(file_path)[1].lower()

    @staticmethod
    def validate_fasta(file_path: str) -> bool:
        """验证FASTA文件格式|Validate FASTA file format"""
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except Exception:
            return False

    @staticmethod
    def validate_bam(file_path: str) -> bool:
        """验证BAM文件格式|Validate BAM file format"""
        try:
            # 使用samtools检查BAM文件头
            result = subprocess.run(
                ["samtools", "view", "-H", file_path],
                capture_output=True,
                text=True,
                timeout=30,
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )
            return result.returncode == 0
        except Exception:
            return False

    @staticmethod
    def count_sequences_in_fasta(fasta_path: str) -> int:
        """统计FASTA文件中的序列数量|Count sequences in FASTA file"""
        try:
            count = 0
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count
        except Exception:
            return 0

    @staticmethod
    def create_temp_dir(prefix: str = "haphic_", base_dir: Optional[str] = None) -> Optional[str]:
        """
        创建临时目录|Create temporary directory

        Args:
            prefix: 临时目录前缀|Temp dir prefix
            base_dir: 父目录;传入则在该目录下创建(超算场景应传 output/tmp 避开系统 /tmp);
                      为 None 时回退系统默认临时目录(向后兼容)|Parent dir; when given,
                      creates under it (on HPC pass output/tmp to avoid filling system /tmp);
                      None falls back to system default (backward-compatible)

        Returns:
            临时目录路径或 None|Temp dir path or None
        """
        try:
            # 先确保 base_dir 存在,否则 mkdtemp(dir=...) 在父目录缺失时会失败
            # Ensure base_dir exists; mkdtemp(dir=...) fails if parent missing
            if base_dir:
                os.makedirs(base_dir, exist_ok=True)
                return tempfile.mkdtemp(prefix=prefix, dir=base_dir)
            # 老 caller 不传 base_dir → 回退系统默认临时目录(向后兼容)
            # Old callers omit base_dir → fall back to system default (backward-compatible)
            return tempfile.mkdtemp(prefix=prefix)
        except Exception as e:
            logging.error(f"无法创建临时目录|Cannot create temporary directory: {e}")
            return None


class QualityController:
    """质量控制器|Quality Controller"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def validate_inputs(self) -> bool:
        """验证输入文件|Validate input files"""
        self.logger.info("验证输入文件|Validating input files")

        # 验证基因组文件
        if not FileManager.validate_fasta(self.config.asm_file):
            self.logger.error(f"基因组文件格式错误|Invalid assembly file format: {self.config.asm_file}")
            return False

        seq_count = FileManager.count_sequences_in_fasta(self.config.asm_file)
        self.logger.info(f"基因组序列数量|Assembly sequence count: {seq_count}")

        # 根据输入类型验证Hi-C文件|Validate Hi-C file based on input type
        if self.config.hic_file_type == "bam":
            # 验证BAM文件|Validate BAM file
            if not self.config.bam_file:
                self.logger.error("BAM文件路径为空|BAM file path is empty")
                return False

            if not FileManager.validate_bam(self.config.bam_file):
                self.logger.error(f"BAM文件格式错误|Invalid BAM file format: {self.config.bam_file}")
                return False

            bam_size = FileManager.get_file_size(self.config.bam_file)
            self.logger.info(f"BAM文件大小|BAM file size: {bam_size / (1024**3):.2f} GB")

        elif self.config.hic_file_type == "fastq":
            # 验证FASTQ文件|Validate FASTQ files
            if not self.config.hic_file:
                self.logger.error("Hi-C FASTQ文件路径为空|Hi-C FASTQ file path is empty")
                return False

            # 验证第一个FASTQ文件|Validate first FASTQ file
            if not FileManager.check_file_exists(self.config.hic_file):
                self.logger.error(f"Hi-C FASTQ文件不存在|Hi-C FASTQ file not found: {self.config.hic_file}")
                return False

            # 验证第二个FASTQ文件 (如果存在)|Validate second FASTQ file (if exists)
            if hasattr(self.config, 'hic2_file') and self.config.hic2_file:
                if not FileManager.check_file_exists(self.config.hic2_file):
                    self.logger.error(f"Hi-C FASTQ第二个文件不存在|Hi-C FASTQ second file not found: {self.config.hic2_file}")
                    return False

                fastq1_size = FileManager.get_file_size(self.config.hic_file)
                fastq2_size = FileManager.get_file_size(self.config.hic2_file)
                self.logger.info(f"FASTQ文件大小|FASTQ file sizes: R1 {fastq1_size / (1024**3):.2f} GB, R2 {fastq2_size / (1024**3):.2f} GB")
            else:
                fastq_size = FileManager.get_file_size(self.config.hic_file)
                self.logger.info(f"FASTQ文件大小|FASTQ file size: {fastq_size / (1024**3):.2f} GB")

        # 验证HapHiC工具
        if not os.path.exists(self.config.haphic_bin):
            self.logger.error(f"HapHiC工具不存在|HapHiC tool not found: {self.config.haphic_bin}")
            return False

        self.logger.info(f"输入文件验证通过|Input files validated")
        return True

    def check_system_resources(self) -> bool:
        """检查系统资源|Check system resources"""
        self.logger.info("检查系统资源|Checking system resources")

        try:
            # 检查内存
            import psutil
            memory = psutil.virtual_memory()
            available_gb = memory.available / (1024**3)
            self.logger.info(f"可用内存|Available memory: {available_gb:.1f} GB")

            # 检查磁盘空间
            disk_usage = psutil.disk_usage(self.config.output_dir)
            free_gb = disk_usage.free / (1024**3)
            self.logger.info(f"可用磁盘空间|Available disk space: {free_gb:.1f} GB")

            # 检查CPU核心数
            cpu_count = psutil.cpu_count()
            self.logger.info(f"CPU核心数|CPU cores: {cpu_count}")

            # 基本资源检查
            if available_gb < 4:
                self.logger.warning("可用内存不足4GB，可能影响运行|Available memory less than 4GB")

            if free_gb < 10:
                self.logger.warning("可用磁盘空间不足10GB，可能影响输出|Available disk space less than 10GB")

            if cpu_count < self.config.threads:
                self.logger.warning(f"请求的线程数({self.config.threads})超过CPU核心数({cpu_count})|Requested threads exceed CPU cores")

            return True

        except ImportError:
            self.logger.info("psutil未安装，跳过系统资源检查|psutil not installed, skipping system resource check")
            return True

        except Exception as e:
            self.logger.warning(f"系统资源检查失败|System resource check failed: {e}")
            return True  # 继续执行，不因为资源检查失败而终止


class ResultValidator:
    """结果验证器|Result Validator"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def validate_output_files(self) -> bool:
        """验证输出文件|Validate output files"""
        self.logger.info("验证输出文件|Validating output files")

        output_files = self.config.get_output_files()
        success = True

        # 验证主要输出文件
        main_files = {
            "scaffolds_fa": output_files["scaffolds_fasta"],
            "scaffolds_agp": output_files["scaffolds_agp"]
        }

        for file_type, file_path in main_files.items():
            if FileManager.check_file_exists(file_path):
                size = FileManager.get_file_size(file_path)
                self.logger.info(f"{file_type}: {file_path} ({size:,} bytes)")
            else:
                self.logger.warning(f"{file_type} 未生成: {file_path}")
                success = False

        # 验证可选输出文件
        optional_files = {
            "scaffolds_raw_agp": output_files["scaffolds_raw_agp"],
            "juicebox_script": output_files["juicebox_script"],
            "contact_map": output_files["contact_map"]
        }

        for file_type, file_path in optional_files.items():
            if FileManager.check_file_exists(file_path):
                size = FileManager.get_file_size(file_path)
                self.logger.info(f"{file_type}: {file_path} ({size:,} bytes)")

        return success

    def get_statistics(self) -> Dict[str, Any]:
        """获取输出统计信息|Get output statistics"""
        output_files = self.config.get_output_files()
        stats = {}

        # 统计scaffolds信息
        if FileManager.check_file_exists(output_files["scaffolds_fasta"]):
            seq_count = FileManager.count_sequences_in_fasta(output_files["scaffolds_fasta"])
            stats["scaffold_count"] = seq_count

        # 统计文件大小
        for file_type, file_path in output_files.items():
            if FileManager.check_file_exists(file_path):
                stats[f"{file_type}_size"] = FileManager.get_file_size(file_path)

        return stats


def format_time(seconds: float) -> str:
    """格式化时间|Format time"""
    if seconds < 60:
        return f"{seconds:.1f}秒"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}分{secs:.1f}秒"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}小时{minutes}分{secs:.1f}秒"


def format_size(bytes_size: int) -> str:
    """格式化文件大小|Format file size"""
    if bytes_size < 1024:
        return f"{bytes_size} B"
    elif bytes_size < 1024**2:
        return f"{bytes_size / 1024:.1f} KB"
    elif bytes_size < 1024**3:
        return f"{bytes_size / (1024**2):.1f} MB"
    else:
        return f"{bytes_size / (1024**3):.1f} GB"


def generate_software_versions_yml(
    config,
    output_file: str,
    logger: logging.Logger
) -> None:
    """
    生成software_versions.yml文件|Generate software_versions.yml file

    Args:
        config: 配置对象(需实现get_software_info)|Configuration object (must implement get_software_info)
        output_file: 输出文件路径|Output file path
        logger: 日志器|Logger
    """
    import yaml
    from datetime import datetime

    try:
        info = config.get_software_info()
        info['execution'] = {
            'timestamp': datetime.now().isoformat()
        }

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, sort_keys=False)

        logger.info(f"版本信息已保存|Version info saved: {output_file}")

    except Exception as e:
        logger.warning(f"生成版本信息文件失败|Failed to generate version info: {e}")


def create_step_directory(output_dir: str, step_name: str) -> str:
    """创建步骤目录|Create step directory"""
    step_dir = os.path.join(output_dir, step_name)
    FileManager.ensure_directory(step_dir)
    return step_dir


def get_step_progress_log(step_name: str, total: int, current: int = 0) -> str:
    """获取步骤进度日志|Get step progress log"""
    progress = (current / total * 100) if total > 0 else 0
    return f"[{step_name}] 进度: {current}/{total} ({progress:.1f}%)"


# 导入aligner模块|Import aligner module
from .aligner import BWAAligner
