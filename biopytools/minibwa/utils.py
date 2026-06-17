"""
Minibwa工具函数模块|Minibwa Utility Functions Module
"""

import glob
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


class MinibwaLogger:
    """Minibwa日志管理器|Minibwa Logger Manager

    按规范2.3分离stdout/stderr/file：
    - stdout: INFO及以上 → 超算.out
    - stderr: WARNING及以上 → 超算.err
    - file:   DEBUG及以上 → 本地完整日志
    """

    def __init__(self, log_file: Path, log_name: str = "minibwa_pipeline.log"):
        self.log_dir = log_file if isinstance(log_file, Path) and log_file.is_dir() else Path(log_file).parent
        self.log_file = self.log_dir / log_name if (self.log_dir / log_name).parent == self.log_dir else log_file
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 清理旧日志|Clear old log
        if isinstance(self.log_file, Path) and self.log_file.exists():
            self.log_file.unlink()

        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 配置命名logger，避免污染root|Use named logger, don't touch root
        logger = logging.getLogger('biopytools.minibwa')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        # 文件handler|File handler (DEBUG+)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # stdout handler|Stdout handler (INFO+)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler|Stderr handler (WARNING+)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner

    支持两种执行方式：
    - run(cmd_list): 单条命令，shell=False
    - run_pipeline(cmd_str): 含管道或重定向的命令，shell=True
    """

    def __init__(self, logger, working_dir: Optional[Path] = None):
        self.logger = logger
        self.working_dir = str(working_dir) if working_dir else None

    def run(self, cmd: List[str], description: str = "") -> bool:
        """
        执行单条命令（无管道）|Execute single command (no pipe)

        Args:
            cmd: 由build_conda_command()构建的命令列表|Command list
            description: 步骤描述|Step description

        Returns:
            bool: 执行是否成功|Whether succeeded
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        # 规范2.2.1: 完整命令必须记录到INFO级别|Full command at INFO level
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir,
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed: {description}")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误输出|Stderr: {result.stderr.strip()}")
                return False

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout.strip()}")
            self.logger.info(f"命令执行成功|Command succeeded: {description}")
            return True

        except Exception as e:
            self.logger.error(f"命令执行异常|Command exception: {description}")
            self.logger.error(f"异常信息|Exception: {e}")
            return False

    def run_pipeline(self, cmd_str: str, description: str = "") -> bool:
        """
        执行含管道的命令|Execute pipelined command

        用于minibwa map | samtools sort这类必须用shell管道的场景。

        Args:
            cmd_str: 完整shell命令字符串|Full shell command string
            description: 步骤描述|Step description

        Returns:
            bool: 执行是否成功|Whether succeeded
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd_str,
                shell=True,
                executable='/bin/bash',
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir,
            )

            if result.returncode != 0:
                self.logger.error(f"管道执行失败|Pipeline failed: {description}")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误输出|Stderr: {result.stderr.strip()}")
                return False

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout.strip()}")
            self.logger.info(f"管道执行成功|Pipeline succeeded: {description}")
            return True

        except Exception as e:
            self.logger.error(f"管道执行异常|Pipeline exception: {description}")
            self.logger.error(f"异常信息|Exception: {e}")
            return False


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中|Detect conda environment from command path

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名或None|Conda env name or None
    """
    # 绝对路径：只检查 /envs/ 模式，不匹配则返回None（尊重用户显式指定的路径）
    # |Absolute path: only check /envs/ pattern, return None if no match
    # (respect user's explicit path choice)
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        return match.group(1) if match else None

    # 命令名：先用shutil.which解析实际路径|Command name: resolve via shutil.which
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 命令名未在PATH中找到，才搜索conda envs兜底
    # |Only fallback to envs search if name not in PATH
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run包装命令|Build conda run wrapped command

    严格按规范13：传递完整路径，自动检测环境，添加--no-capture-output
    |Per spec 13: pass full path, auto-detect env, add --no-capture-output

    Args:
        command: 完整命令路径|Full command path
        args: 参数列表|Argument list

    Returns:
        包装后的命令列表|Wrapped command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        # conda run使用命令名即可，因为已经通过-n指定了环境
        # |conda run uses command name, env is specified via -n
        cmd_name = os.path.basename(command)
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', cmd_name] + args
    # 非conda环境，直接用完整路径|Non-conda, use full path directly
    return [command] + args


def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    deps = [
        (config.minibwa_path, 'minibwa', ['version']),
        (config.samtools_path, 'samtools', ['--version']),
    ]

    missing = []
    for path, name, version_args in deps:
        try:
            # 用build_conda_command包装，支持conda env内的samtools
            # |Wrap with build_conda_command to support conda env samtools
            cmd = build_conda_command(path, version_args)
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=30
            )
            if result.returncode == 0:
                first_line = (result.stdout or result.stderr).strip().split('\n')[0]
                logger.info(f"{name} 可用|{name} available: {first_line}")
            else:
                logger.warning(f"{name}版本检查返回非零|{name} version check nonzero: rc={result.returncode}")
                missing.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"{name}不可用|{name} not available: {e}")
            missing.append(name)

    if missing:
        raise RuntimeError(f"缺少依赖软件|Missing dependencies: {', '.join(missing)}")

    return True


def find_fastq_pairs(input_dir: str, pattern: str, logger) -> List[Tuple[str, str, str]]:
    """
    查找FASTQ配对文件|Find paired FASTQ files

    Args:
        input_dir: FASTQ目录|FASTQ directory
        pattern: R1匹配模式（如_1.fq.gz或_1.clean.fq.gz）|R1 pattern
        logger: 日志器|Logger

    Returns:
        List of (sample_name, read1_path, read2_path)
    """
    logger.info(f"查找FASTQ文件|Finding FASTQ files")
    logger.info(f"输入目录|Input directory: {input_dir}")
    logger.info(f"匹配模式|Pattern: {pattern}")

    search = os.path.join(input_dir, f"*{pattern}")
    r1_files = sorted(glob.glob(search))

    if not r1_files:
        raise ValueError(f"未找到匹配的FASTQ文件|No FASTQ found matching: {search}")

    # 推断R2模式：将pattern中的_1替换为_2|Infer R2 pattern
    if '_1.' in pattern:
        r2_pattern = pattern.replace('_1.', '_2.', 1)
    elif '_1' in pattern:
        r2_pattern = pattern.replace('_1', '_2', 1)
    elif '1.fq' in pattern:
        r2_pattern = pattern.replace('1.fq', '2.fq', 1)
    elif 'R1' in pattern:
        r2_pattern = pattern.replace('R1', 'R2', 1)
    else:
        r2_pattern = None

    pairs = []
    skipped = []
    for r1 in r1_files:
        basename = os.path.basename(r1)
        sample_name = basename.replace(pattern, '', 1)
        if not sample_name:
            sample_name = Path(basename).stem

        if r2_pattern is None:
            r2 = None
        else:
            r2 = r1.replace(pattern, r2_pattern, 1)
            if not os.path.exists(r2):
                skipped.append((sample_name, r2))
                r2 = None

        pairs.append((sample_name, r1, r2))

    for sample_name, r1, r2 in pairs:
        if r2:
            logger.info(f"找到样品|Found sample (PE): {sample_name}")
        else:
            logger.info(f"找到样品|Found sample (SE): {sample_name}")

    for sample_name, missing_r2 in skipped:
        logger.warning(f"样品|Sample {sample_name}: 未找到R2|missing R2: {missing_r2}, 按单端处理|treating as single-end")

    logger.info(f"共找到|Total: {len(pairs)} 个样品|samples")
    return pairs


def format_number(n: int) -> str:
    """格式化大数字（M单位）|Format large number with M unit"""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n / 1_000:.2f}K"
    return str(n)
