"""
叶绿体基因组组装工具模块|Plastome Assembly Utility Module
"""

import os
import subprocess
import logging
import sys
import time
import glob
import re
import shutil
import fnmatch
from pathlib import Path
from typing import Optional, Tuple, List


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
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    # 检查是否在conda环境中|Check if in conda environment
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用 conda run|Use conda run
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 直接调用|Direct call
        full_cmd = [command] + args

    return full_cmd


class PlastomeLogger:
    """叶绿体组装日志管理器|Plastome Assembly Logger Manager"""

    def __init__(self, output_dir: str, verbose: bool = False, log_file: Optional[str] = None):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            verbose: 详细输出模式|Verbose output mode
            log_file: 日志文件路径|Log file path
        """
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        self.log_file = log_file
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """设置日志系统|Setup logging system"""
        logger = logging.getLogger('PlastomeAssembler')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO 及以下级别|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件 handler (如果指定)|File handler (if specified)
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取 logger 对象|Get logger object"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str, conda_env_path: Optional[str] = None,
                 working_dir: Optional[str] = None):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
            output_dir: 输出目录|Output directory
            conda_env_path: conda环境路径|Conda environment path
            working_dir: 工作目录|Working directory (for avoiding Chinese characters in path)
        """
        self.logger = logger
        self.output_dir = Path(output_dir)
        self.conda_env_path = conda_env_path
        self.working_dir = working_dir

    def _get_env_with_conda(self) -> dict:
        """
        获取包含conda环境的环境变量|Get environment variables with conda environment

        Returns:
            dict: 环境变量字典|Environment variables dictionary
        """
        import copy
        env = copy.deepcopy(os.environ)

        if self.conda_env_path:
            # 将conda环境的bin目录添加到PATH|Add conda environment bin directory to PATH
            conda_bin = os.path.join(self.conda_env_path, 'bin')
            if 'PATH' in env:
                env['PATH'] = f"{conda_bin}:{env['PATH']}"
            else:
                env['PATH'] = conda_bin

            self.logger.debug(f"设置PATH环境变量|Setting PATH: {env['PATH'][:200]}...")

        return env

    def run_command(self, command: List[str], description: str = "") -> bool:
        """
        执行命令|Run command

        Args:
            command: 命令列表|Command list
            description: 命令描述|Command description

        Returns:
            bool: 是否成功|Whether successful
        """
        if description:
            self.logger.info(f"执行中|Executing: {description}")

        # 自动包装conda环境的命令|Automatically wrap conda environment commands
        if len(command) > 0:
            cmd_name = os.path.basename(command[0])
            wrapped_cmd = build_conda_command(cmd_name, command[1:])
        else:
            wrapped_cmd = command

        cmd_str = ' '.join(wrapped_cmd)
        self.logger.debug(f"完整命令|Full command: {cmd_str}")

        # 确定工作目录|Determine working directory
        cwd = self.working_dir if self.working_dir else None
        if cwd:
            self.logger.debug(f"工作目录|Working directory: {cwd}")

        start_time = time.time()

        try:
            result = subprocess.run(
                wrapped_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                cwd=cwd  # 设置工作目录|Set working directory
            )

            elapsed_time = time.time() - start_time

            if result.stdout:
                # 输出到 stdout (INFO 级别)|Output to stdout (INFO level)
                for line in result.stdout.strip().split('\n'):
                    if line:
                        self.logger.info(line)

            if result.stderr:
                # 输出到 stderr (WARNING 级别)|Output to stderr (WARNING level)
                for line in result.stderr.strip().split('\n'):
                    if line:
                        self.logger.warning(line)

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed with return code {result.returncode}")
                return False

            self.logger.info(f"命令执行完成|Command completed in {elapsed_time:.2f} 秒|seconds")
            return True

        except Exception as e:
            elapsed_time = time.time() - start_time
            self.logger.error(f"命令执行异常|Command execution exception: {str(e)}")
            self.logger.error(f"执行时间|Execution time: {elapsed_time:.2f} 秒|seconds")
            return False


def extract_sample_name(filename: str, pattern: str) -> str:
    """
    从文件名中提取样本名，支持通配符模式|Extract sample name from filename, supports wildcard patterns

    Args:
        filename: 文件名|Filename
        pattern: 文件名模式（可包含通配符）|Filename pattern (may contain wildcards)

    Returns:
        提取的样本名|Extracted sample name
    """
    # 如果模式不包含通配符，直接替换|If pattern doesn't contain wildcard, replace directly
    if '*' not in pattern and '?' not in pattern:
        return filename.replace(pattern, "")

    # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
    # 将 * 替换为 (.+) 捕获组以提取样本名|Replace * with (.+) capture group to extract sample name
    regex_pattern = fnmatch.translate(pattern)
    # 修改正则表达式，将通配符部分设为捕获组|Modify regex to make wildcard part a capture group
    regex_pattern = regex_pattern.replace(r'\.\*.*', r'(.+)\..*')
    regex_pattern = regex_pattern.replace(r'\Z', r'\Z')

    match = re.match(regex_pattern, filename)
    if match:
        return match.group(1)

    # 如果正则匹配失败，回退到简单替换|If regex match fails, fallback to simple replacement
    # 去除已知的后缀|Remove known suffix
    name = filename
    # 尝试去除模式中通配符后面的部分|Try to remove part after wildcard in pattern
    if '*' in pattern:
        suffix = pattern.split('*', 1)[1]
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


class SampleFinder:
    """样本查找器|Sample Finder"""

    def __init__(self, input_dir: str, read1_suffix: str, read2_suffix: str,
                 logger: Optional[logging.Logger] = None):
        """
        初始化样本查找器|Initialize sample finder

        Args:
            input_dir: 输入目录路径|Input directory path
            read1_suffix: R1文件后缀模式|R1 file suffix pattern
            read2_suffix: R2文件后缀模式|R2 file suffix pattern
            logger: 日志对象|Logger object
        """
        # 保持相对路径|Keep relative path to avoid Chinese characters in absolute path
        self.input_dir = Path(input_dir)
        self.read1_suffix = read1_suffix
        self.read2_suffix = read2_suffix
        self.logger = logger

    def find_sample_pairs(self) -> List[Tuple[str, Path, Path]]:
        """
        查找样本配对文件|Find sample paired files

        Returns:
            样本配对列表，每个元素为 (样本名, read1文件, read2文件)|List of sample pairs (sample_name, read1_file, read2_file)
        """
        sample_pairs = []

        # 查找所有read1文件|Find all read1 files
        # 如果后缀已经包含通配符，直接使用；否则在前面加*|If suffix contains wildcard, use directly; otherwise add *
        if '*' in self.read1_suffix:
            read1_pattern = self.read1_suffix
        else:
            read1_pattern = f"*{self.read1_suffix}"
        read1_files = list(self.input_dir.glob(read1_pattern))

        if not read1_files:
            if self.logger:
                self.logger.warning(
                    f"在输入目录中未找到匹配模式的文件|"
                    f"No files matching pattern found in input directory: {read1_pattern}"
                )
            return sample_pairs

        if self.logger:
            self.logger.info(f"找到 {len(read1_files)} 个read1文件|Found {len(read1_files)} read1 files")

        for read1_file in read1_files:
            # 提取样本名|Extract sample name
            sample_name = extract_sample_name(read1_file.name, self.read1_suffix)

            # 构建read2文件路径|Construct read2 file path
            read2_file = self.input_dir / f"{sample_name}{self.read2_suffix}"

            if read2_file.exists():
                sample_pairs.append((sample_name, read1_file, read2_file))
                if self.logger:
                    self.logger.debug(f"找到配对样本|Found paired sample: {sample_name}")
            else:
                # 尝试自动推导配对文件名|Try to auto-derive paired filename
                read2_file = self._find_paired_file(read1_file)
                if read2_file:
                    sample_pairs.append((sample_name, read1_file, read2_file))
                    if self.logger:
                        self.logger.debug(f"通过自动推导找到配对样本|Found paired sample via auto-detection: {sample_name}")
                else:
                    if self.logger:
                        self.logger.warning(
                            f"找不到样本 {sample_name} 的配对文件|"
                            f"Cannot find paired file for sample {sample_name}"
                        )

        return sample_pairs

    def _find_paired_file(self, read1_file: Path) -> Optional[Path]:
        """
        自动查找配对的Read2文件|Auto-find paired Read2 file

        Args:
            read1_file: Read1文件路径|Read1 file path

        Returns:
            Read2文件路径（如果存在）|Read2 file path (if exists), None otherwise
        """
        # 获取文件所在的目录|Get the directory of the file
        file_dir = read1_file.parent
        file_name = read1_file.name

        # 尝试多种常见的配对模式|Try multiple common pairing patterns
        patterns_to_try = [
            # 替换 _1 为 _2 | Replace _1 with _2
            file_name.replace('_1.', '_2.'),
            # 替换 .R1 为 .R2 | Replace .R1 with .R2
            file_name.replace('.R1.', '.R2.'),
            # 替换 _R1 为 _R2 | Replace _R1 with _R2
            file_name.replace('_R1.', '_R2.'),
            # 替换 .r1 为 .r2 | Replace .r1 with .r2
            file_name.replace('.r1.', '.r2.'),
            # 替换 _r1 为 _r2 | Replace _r1 with _r2
            file_name.replace('_r1.', '_r2.'),
            # 替换 /read1 为 /read2 | Replace /read1 with /read2
            file_name.replace('read1', 'read2'),
        ]

        for pattern in patterns_to_try:
            if pattern != file_name:  # 确保模式确实改变了文件名|Ensure pattern actually changed the filename
                paired_file = file_dir / pattern
                if paired_file.exists():
                    if self.logger:
                        self.logger.debug(f"使用模式找到配对文件|Found paired file using pattern: {pattern}")
                    return paired_file

        return None


def detect_samples_and_reads(input_dir: str, logger: Optional[logging.Logger] = None,
                             read1_suffix: str = '_1.clean.fq.gz',
                             read2_suffix: str = '_2.clean.fq.gz') -> dict:
    """
    自动检测输入目录中的多个样品及其reads文件|Auto-detect multiple samples and their reads files

    Args:
        input_dir: 输入目录路径|Input directory path
        logger: 日志对象|Logger object
        read1_suffix: R1文件后缀模式|R1 file suffix pattern
        read2_suffix: R2文件后缀模式|R2 file suffix pattern

    Returns:
        dict: 样品信息字典|Sample information dictionary
              {
                  'sample_name': {
                      'r1': 'path/to/R1.fq.gz',
                      'r2': 'path/to/R2.fq.gz',
                  }
              }
    """
    if logger:
        logger.info(f"自动检测样品和reads文件|Auto-detecting samples and reads files in: {input_dir}")

    # 使用 SampleFinder 查找样品配对|Use SampleFinder to find sample pairs
    finder = SampleFinder(input_dir, read1_suffix, read2_suffix, logger)
    sample_pairs = finder.find_sample_pairs()

    if not sample_pairs:
        if logger:
            logger.error(f"未找到任何有效的样品配对|No valid sample pairs found in {input_dir}")
        return {}

    # 转换为字典格式|Convert to dictionary format
    samples = {}
    for sample_name, r1_file, r2_file in sample_pairs:
        samples[sample_name] = {
            'r1': str(r1_file),
            'r2': str(r2_file) if r2_file else None,
        }

    # 输出检测结果|Output detection results
    if logger:
        logger.info(f"检测到 {len(samples)} 个样品|Detected {len(samples)} samples:")
        for sample_name, files in sorted(samples.items()):
            logger.info(f"  样品|Sample: {sample_name}")
            if files['r1']:
                logger.info(f"    R1: {os.path.basename(files['r1'])}")
            if files['r2']:
                logger.info(f"    R2: {os.path.basename(files['r2'])}")

    return samples


def detect_reads_files(input_dir: str, logger: Optional[logging.Logger] = None,
                       read1_suffix: str = '_1.clean.fq.gz',
                       read2_suffix: str = '_2.clean.fq.gz') -> Tuple[Optional[str], Optional[str], None]:
    """
    自动检测输入目录中的reads文件（单个样品）| Auto-detect reads files in input directory (single sample)

    Args:
        input_dir: 输入目录路径|Input directory path
        logger: 日志对象|Logger object
        read1_suffix: R1文件后缀模式|R1 file suffix pattern
        read2_suffix: R2文件后缀模式|R2 file suffix pattern

    Returns:
        Tuple[r1_file, r2_file, None]: R1文件, R2文件, None (保留第三个参数以兼容旧接口)
    """
    samples = detect_samples_and_reads(input_dir, logger, read1_suffix, read2_suffix)

    if not samples:
        return None, None, None

    # 如果只有一个样品，返回其文件信息|If only one sample, return its file info
    if len(samples) == 1:
        sample_name = list(samples.keys())[0]
        files = samples[sample_name]
        return files['r1'], files['r2'], None

    # 如果有多个样品，使用第一个|If multiple samples, use the first one
    sample_name = sorted(samples.keys())[0]
    files = samples[sample_name]

    if logger:
        logger.warning(f"检测到多个样品，仅处理第一个: {sample_name}|Detected multiple samples, only processing the first: {sample_name}")

    return files['r1'], files['r2'], None
