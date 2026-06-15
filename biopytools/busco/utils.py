"""
BUSCO分析工具函数模块|BUSCO Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import glob
import os
import re
import shutil
from pathlib import Path
from typing import List, Tuple, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path (e.g., 'busco' or '/path/to/busco')

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先检查传入的command是否本身就是完整路径
    # First check if command itself is a full path
    if os.path.isabs(command):
        # 直接从完整路径中提取环境名|Extract env name directly from full path
        # 例如|e.g.: /miniforge3/envs/BUSCO_v.6.0.0/bin/busco
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 如果不是完整路径，尝试用which查找|If not full path, try which
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        # 例如|e.g.: /miniforge3/envs/BUSCO_v.6.0.0/bin/busco
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令|Search for command in all environments
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
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用，添加--no-capture-output避免缓冲输出导致内存问题
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, direct call
        full_cmd = [command] + args

    return full_cmd


class BUSCOLogger:
    """BUSCO分析日志管理器|BUSCO Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "busco_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / "99_logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logger = logging.getLogger("busco")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - 仅INFO|stdout handler - INFO only
        # 规范§2.3: INFO→stdout→.out, WARNING+→stderr→.err
        # Python setLevel是>=语义，需加filter限制上限|Python setLevel is >=, need filter to cap upper bound
        class _MaxLevelFilter(logging.Filter):
            def __init__(self, max_level):
                super().__init__()
                self.max_level = max_level

            def filter(self, record):
                return record.levelno <= self.max_level

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(_MaxLevelFilter(logging.INFO))
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()

    def run(self, cmd: list, description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False, e.stdout, e.stderr


class FileManager:
    """文件管理器|File Manager"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def get_input_files(self) -> List[Tuple[str, str]]:
        """获取输入文件列表|Get input file list

        Returns:
            List[Tuple[str, str]]: [(file_path, sample_name), ...]
        """
        input_path = Path(self.config.input_path)
        files = []

        if input_path.is_file():
            # 单个文件|Single file
            sample_name = self.extract_sample_name(input_path.name)
            files.append((str(input_path), sample_name))
            self.logger.info(f"发现单个输入文件|Found single input file: {input_path.name}")

        elif input_path.is_dir():
            # 目录批处理|Directory batch processing
            pattern = self.config.sample_suffix.replace('*', '*')
            search_pattern = input_path / pattern

            matching_files = glob.glob(str(search_pattern))
            if not matching_files:
                # 尝试常见序列文件扩展名，合并所有匹配结果|Try common extensions, merge all matches
                fallback_patterns = [
                    '*.fa', '*.faa', '*.pep', '*.pep.fa', '*.aa',
                    '*.fasta', '*.protein', '*.protein.fa',
                    '*.fna', '*.cds', '*.cds.fa', '*.gene.fa',
                    '*.gff', '*.gff3',
                ]
                for fb in fallback_patterns:
                    fb_files = glob.glob(str(input_path / fb))
                    if fb_files:
                        matching_files.extend(fb_files)
                        self.logger.info(f"匹配到文件|Matched files: {len(fb_files)} 个文件|files ({fb})")

            # 去重|Deduplicate
            matching_files = sorted(set(matching_files))

            if not matching_files:
                raise ValueError(f"目录中未找到匹配模式的文件|No files matching pattern '{pattern}' found in directory")

            for file_path in sorted(matching_files):
                file_name = Path(file_path).name
                sample_name = self.extract_sample_name(file_name)
                files.append((file_path, sample_name))

            self.logger.info(f"发现批处理文件|Found batch files: {len(files)} 个文件|files")

        else:
            raise ValueError(f"输入路径无效|Invalid input path: {input_path}")

        return files

    def extract_sample_name(self, filename: str) -> str:
        """从文件名提取样本名|Extract sample name from filename"""
        suffix_pattern = self.config.sample_suffix

        # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
        escaped_pattern = re.escape(suffix_pattern)
        regex_pattern = escaped_pattern.replace(r'\*', r'(.*)')

        match = re.match(regex_pattern, filename)
        if match:
            return match.group(1)
        else:
            # 如果模式不匹配，使用文件名去除扩展名|If pattern doesn't match, use filename without extension
            return Path(filename).stem


def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        cmd = build_conda_command(config.busco_path, ['--version'])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            version_info = result.stdout.strip()
            logger.info(f"BUSCO 可用|BUSCO available: {version_info}")
        else:
            raise RuntimeError("BUSCO版本检查失败|BUSCO version check failed")
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"BUSCO不可用|BUSCO not available: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True
