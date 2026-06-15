"""
OrthoFinder泛基因组分析工具函数模块|OrthoFinder Pangenome Analysis Utility Functions Module
"""

import logging
import re
import shutil
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Optional, Tuple


class PangenomeLogger:
    """泛基因组分析日志管理器|Pangenome Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "pangenome_analysis.log"):
        self.output_dir = output_dir
        self.log_dir = output_dir / "99_logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.log_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        logger = logging.getLogger("orthofinder_pangenome")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
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
        self.tmp_dir = self.working_dir / "tmp"
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"临时目录|Temp directory: {self.tmp_dir}")

    def run(self, cmd: list, description: str = "", long_running: bool = False) -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
            long_running: 是否为长时间运行命令|Whether the command runs for a long time
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        try:
            if long_running:
                return self._run_long_command(cmd, description)
            else:
                return self._run_short_command(cmd, description)

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {description}: {e}")
            return False

    def _run_short_command(self, cmd: list, description: str = "") -> bool:
        """执行短命令（捕获输出）|Execute short command (capture output)"""
        env = os.environ.copy()
        env['TMPDIR'] = str(self.tmp_dir)

        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            check=True,
            cwd=self.working_dir,
            env=env
        )

        self.logger.info(f"命令执行成功|Command executed successfully: {description}")

        if result.stdout:
            self.logger.debug(f"标准输出|Stdout: {result.stdout}")

        return True

    def _run_long_command(self, cmd: list, description: str = "") -> bool:
        """执行长时间命令（流式输出）|Execute long-running command (streaming output)"""
        import selectors

        env = os.environ.copy()
        env['TMPDIR'] = str(self.tmp_dir)

        proc = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=self.working_dir,
            env=env
        )

        sel = selectors.DefaultSelector()
        sel.register(proc.stdout, selectors.EVENT_READ)
        sel.register(proc.stderr, selectors.EVENT_READ)

        stdout_lines = []
        stderr_lines = []

        while proc.poll() is None:
            events = sel.select(timeout=1)
            for key, _ in events:
                line = key.fileobj.readline()
                if not line:
                    continue
                decoded = line.decode('utf-8', errors='replace').rstrip()
                if not decoded:
                    continue
                if key.fileobj is proc.stdout:
                    stdout_lines.append(decoded)
                    self.logger.debug(f"[stdout] {decoded}")
                else:
                    stderr_lines.append(decoded)
                    self.logger.info(f"[stderr] {decoded}")

        sel.close()

        if proc.returncode != 0:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {proc.returncode}")
            if stderr_lines:
                self.logger.error(f"错误信息|Error message:\n" + "\n".join(stderr_lines))
            return False

        self.logger.info(f"命令执行成功|Command executed successfully: {description}")
        return True


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda env, return env name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|Conda env name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        real_path = os.path.realpath(cmd_path)
        match = re.search(r'/envs/([^/]+)', real_path)
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
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to invoke software in conda env

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    cmd = build_conda_command(config.orthofinder_path, ['-h'])
    logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        output = result.stdout + result.stderr
        if "OrthoFinder" in output or "orthofinder" in output:
            logger.info(f"OrthoFinder 可用|OrthoFinder available")
        else:
            raise RuntimeError(f"程序响应异常|Program response abnormal: {config.orthofinder_path}")

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"程序响应超时|Program response timeout: {config.orthofinder_path}")
    except FileNotFoundError:
        raise RuntimeError(f"未找到OrthoFinder|OrthoFinder not found: {config.orthofinder_path}")

    return True


def count_fasta_files(input_dir: str) -> Tuple[int, List[str]]:
    """统计FASTA文件数量和列表|Count FASTA files and get file list"""
    input_path = Path(input_dir)

    fasta_extensions = ['*.fa', '*.faa', '*.fas', '*.fasta', '*.pep', '*.protein']

    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(input_path.glob(ext))

    sample_names = [f.stem for f in fasta_files]

    return len(fasta_files), sample_names


def validate_fasta_files(input_dir: str, logger) -> bool:
    """验证FASTA文件格式|Validate FASTA file format"""
    logger.info("验证输入文件格式|Validating input file format")

    count, sample_names = count_fasta_files(input_dir)

    if count == 0:
        logger.error("未找到FASTA文件|No FASTA files found")
        return False

    logger.info(f"找到 {count} 个FASTA文件|Found {count} FASTA files")
    logger.info(f"样本名称|Sample names: {', '.join(sample_names[:5])}{'...' if count > 5 else ''}")

    input_path = Path(input_dir)
    empty_files = []
    invalid_files = []

    for sample_name in sample_names:
        possible_files = list(input_path.glob(f"{sample_name}.*"))
        if possible_files:
            file_path = possible_files[0]
            if file_path.stat().st_size == 0:
                empty_files.append(sample_name)
                continue
            for line in open(file_path):
                if line.startswith('>'):
                    continue
                invalid_chars = set(line.strip()) - set('ACDEFGHIKLMNPQRSTVWYXOUBJZ*-')
                if invalid_chars:
                    invalid_files.append((sample_name, invalid_chars))
                    break

    if empty_files:
        logger.warning(f"发现空文件|Found empty files: {', '.join(empty_files)}")

    if invalid_files:
        logger.error(f"发现含非法字符的FASTA文件|Found FASTA files with invalid characters:")
        for name, chars in invalid_files:
            logger.error(f"  {name}: 非法字符|invalid chars: {chars}")
        return False

    return True


def clean_fasta_sequences(input_dir: str, logger) -> int:
    """清洗FASTA序列中的非法字符（仅保留标准氨基酸字符）|Clean invalid characters from FASTA sequences"""
    valid_chars = set('ACDEFGHIKLMNPQRSTVWYXOUBJZ*-')
    cleaned_count = 0
    input_path = Path(input_dir)

    for fasta_file in input_path.glob('*.fa*'):
        modified = False
        lines = []
        for line in open(fasta_file):
            if line.startswith('>'):
                lines.append(line)
            else:
                cleaned = ''.join(c for c in line.strip() if c in valid_chars)
                if len(cleaned) != len(line.strip()):
                    modified = True
                if cleaned:
                    lines.append(cleaned + '\n')

        if modified:
            with open(fasta_file, 'w') as f:
                f.writelines(lines)
            cleaned_count += 1
            logger.info(f"已清洗|Cleaned: {fasta_file.name}")

    if cleaned_count:
        logger.info(f"共清洗 {cleaned_count} 个FASTA文件|Cleaned {cleaned_count} FASTA files")
    else:
        logger.info("无需清洗|No cleaning needed")

    return cleaned_count


def format_number(num: int) -> str:
    """格式化大数字|Format large number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return f"{num:,}"
