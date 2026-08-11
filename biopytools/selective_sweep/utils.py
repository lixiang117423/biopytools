"""
选择性扫荡检测工具函数模块|Selective Sweep Detection Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中|Detect if command is in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 仅当command是裸命令名(不含路径分隔符)时,才搜索所有conda环境。
    # 绝对路径(如静态二进制 ~/software/sweed/SweeD-P)必须跳过,否则
    # os.path.join(envs_dir, env, 'bin', command) 会被绝对路径塌缩为原路径
    # (存在),误判为位于第一个枚举到的env而静默包装另一个同名二进制。
    # |Search all conda envs ONLY for bare command names. Absolute paths MUST
    # skip, else os.path.join collapses to that path (which exists) and falsely
    # reports the first env iterated, silently wrapping a different binary.
    if os.path.sep not in command:
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
    """构建conda run命令(含--no-capture-output)|Build conda run command (with --no-capture-output)

    Args:
        command: 命令名称或完整路径(禁用os.path.basename提取)|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表(配合subprocess.run(shell=False))|Full command list (for shell=False)
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class SweepModuleLogger:
    """选择性扫荡检测日志管理器|Selective Sweep Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = 'selective_sweep.log',
                 log_level: str = 'INFO'):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.log_level = log_level  # 必须在setup_logging之前赋值|assign before setup_logging
        self.setup_logging()

    def setup_logging(self):
        """设置日志(stdout INFO / stderr WARNING+ / file DEBUG+)|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logger = logging.getLogger('selective_sweep')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False  # 避免重复|avoid duplicates

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - log_level级别|stdout handler - log_level
        console_level = getattr(logging, self.log_level.upper(), logging.INFO)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(console_level)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING+|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file handler - DEBUG+|file handler - DEBUG and above
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

    def run(self, cmd: list, description: str = '', cwd: Optional[Path] = None) -> bool:
        """执行命令(记录完整命令到INFO)|Execute command (log full command at INFO)

        Args:
            cmd: 命令列表(由build_conda_command构建)|Command list (built by build_conda_command)
            description: 步骤描述|Step description
            cwd: 可选工作目录(默认self.working_dir)|Optional working dir

        Returns:
            bool: 是否成功|Success flag
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        work_dir = Path(cwd).resolve() if cwd else self.working_dir
        self.logger.info(f"工作目录|Working directory: {work_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=str(work_dir)
            )
            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout:\n{result.stdout}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message:\n{e.stderr}")
            return False

        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return False


def _probe_version(cmd: list, timeout: int = 15) -> str:
    """尽力探测工具版本(失败返回'unknown'),绝不抛异常。
    |Best-effort version probe (returns 'unknown' on failure); never raises.

    版本行选择|Version line selection:
    - 优先取版本行(带点号的version X.Y或RAiSD vX.Y);RAiSD -v输出的是按时间
      排序的changelog横幅,最新版本在最后,故取最后一条版本行而非第一条
      |prefer version lines (dotted "version X.Y" or "RAiSD vX.Y"); RAiSD -v
      prints a chronological changelog banner with the current version last,
      so the LAST version line wins
    - 无版本行(bcftools/vcftools)回退首个非空行
      |no version line (bcftools/vcftools) -> fall back to first non-empty line
    - 版本号必须带点号,避免误匹配bcftools license行"GNU GPL version 3 or later"
      |dotted version required to avoid matching bcftools' license line
    """
    try:
        # capture_output与stderr=STDOUT互斥,改用显式PIPE合并两流|capture_output conflicts
        # with stderr=STDOUT; use explicit PIPE merging both streams
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=timeout,
        )
        if result.returncode != 0:
            return 'unknown'
        combined = (result.stdout or '') + (result.stderr or '')
        lines = [line.strip() for line in combined.splitlines() if line.strip()]
        if not lines:
            return 'unknown'
        version_pat = re.compile(r'version\s*[0-9]+\.[0-9]+|RAiSD\s+v[0-9]+\.[0-9]+',
                                 re.IGNORECASE)
        version_lines = [ln for ln in lines if version_pat.search(ln)]
        if version_lines:
            return version_lines[-1]
        return lines[0]
    except Exception:
        return 'unknown'


def generate_software_versions_yml(config, logger) -> None:
    """生成00_pipeline_info/software_versions.yml(纯文本,不依赖yaml包)。
    |Generate 00_pipeline_info/software_versions.yml as plain text (no yaml dependency)."""
    try:
        versions = {}
        for tool_name, tool_path, ver_args in [
            ('bcftools', config.bcftools_path, ['--version']),
            ('vcftools', config.vcftools_path, ['--version']),
            ('raisd', config.raisd_path, ['-v']),
            ('sweed', config.sweed_path, ['-v']),
        ]:
            probe_cmd = build_conda_command(tool_path, ver_args)
            versions[tool_name] = _probe_version(probe_cmd)

        params = [
            ('win', config.win),
            ('step', config.step),
            ('top_quantile', config.top_quantile),
            ('merge_gap', config.merge_gap),
            ('min_maf', config.min_maf),
            ('max_missing', config.max_missing),
            ('raisd_window', config.raisd_window),
            ('raisd_min_samples', config.raisd_min_samples),
            ('include_mu_low_n', config.include_mu_low_n),
            ('sweed_grid', config.sweed_grid),
            ('sweed_min_samples', config.sweed_min_samples),
            ('include_sweed_low_n', config.include_sweed_low_n),
            ('sweed_folded', config.sweed_folded),
            ('threads', config.threads),
        ]

        lines = [
            '# 00_pipeline_info/software_versions.yml',
            'pipeline:',
            '  name: "biopytools selective_sweep"',
            '  version: "1.0.0"',
            f'  generated: "{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}"',
            'tools:',
        ]
        for name, ver in versions.items():
            lines.append(f'  {name}:')
            lines.append(f'    version: "{ver}"')
            lines.append(f'    path: "{getattr(config, name + "_path")}"')
        lines.append('parameters:')
        for k, v in params:
            lines.append(f'  {k}: {v}')

        out_file = Path(config.info_dir) / 'software_versions.yml'
        out_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        logger.info(f"版本信息已保存|Software versions saved: {out_file}")
    except Exception as e:
        # 版本记录为辅助信息,失败不影响主流程|version record is auxiliary; never break the run
        logger.warning(f"生成software_versions.yml失败|Failed to write software_versions.yml: {e}")


def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")
    all_ok = True
    for tool_name, tool_path, ver_args in [
        ('bcftools', config.bcftools_path, ['--version']),
        ('vcftools', config.vcftools_path, ['--version']),
        ('RAiSD', config.raisd_path, ['-v']),
        ('SweeD-P', config.sweed_path, ['-v']),
    ]:
        try:
            cmd = build_conda_command(tool_path, ver_args)
            logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                logger.info(f"{tool_name} 可用|{tool_name} available")
            else:
                logger.error(f"{tool_name} 不可用|{tool_name} not available (rc={result.returncode})")
                all_ok = False
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"{tool_name} 不可用|{tool_name} not available: {e}")
            all_ok = False
    return all_ok
