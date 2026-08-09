"""
VCF转树工具函数模块|VCF to Tree Utility Functions Module
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
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
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
    """构建conda run命令|Build conda run command

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表(适用于shell=False)|Complete command list (for shell=False)
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class Vcf2TreeLogger:
    """VCF转树日志管理器|VCF to Tree Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "vcf2tree.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logger = logging.getLogger('vcf2tree')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING+
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file handler - DEBUG+
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

    def run(self, cmd: list, description: str = "") -> bool:
        """执行命令|Execute command

        Args:
            cmd: 命令列表(由build_conda_command构建)|Command list (built by build_conda_command)
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

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
                self.logger.debug(f"标准输出|Stdout:\n{result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message:\n{e.stderr}")
            return False


def _probe_version(cmd: list, timeout: int = 15) -> str:
    """尽力探测工具版本(失败返回'unknown'), 绝不抛异常。
    |Best-effort version probe (returns 'unknown' on failure); never raises."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            stderr=subprocess.STDOUT,
        )
        combined = (result.stdout or '') + (result.stderr or '')
        # 取首个非空行作为版本信息|first non-empty line as version info
        for line in combined.splitlines():
            line = line.strip()
            if line:
                return line
        return 'unknown'
    except Exception:
        return 'unknown'


def generate_software_versions_yml(config, logger) -> None:
    """生成00_pipeline_info/software_versions.yml(纯文本, 不依赖yaml包)。
    |Generate 00_pipeline_info/software_versions.yml as plain text (no yaml dependency).

    §12.5合规: 文件名不含版本号, 版本信息集中记录到yml。
    |§12.5 compliance: no version in filenames; record versions centrally in yml.
    """
    try:
        # 探测建树工具版本|probe tree-tool version
        if config.method == 'fasttree':
            tool_name = 'fasttree'
            tool_path = config.fasttree_path
            probe_cmd = build_conda_command(tool_path, [])
        else:
            tool_name = 'iqtree'
            tool_path = config.iqtree_path
            probe_cmd = build_conda_command(tool_path, ['--version'])
        version = _probe_version(probe_cmd)

        params = [
            ('method', config.method),
            ('threads', config.threads),
            ('min_samples_locus', config.min_samples_locus),
            ('outgroup', config.outgroup or '(none)'),
        ]
        if config.method == 'fasttree':
            params.append(('fasttree_params', config.fasttree_params or '(none)'))
        else:
            params.append(('iqtree_model', config.iqtree_model or 'MFP(auto)'))
            params.append(('iqtree_bootstrap', config.iqtree_bootstrap))
            params.append(('iqtree_asc', config.iqtree_asc))

        lines = [
            '# 00_pipeline_info/software_versions.yml',
            'pipeline:',
            '  name: "biopytools vcf2tree"',
            '  version: "1.0.0"',
            f'  generated: "{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}"',
            'tools:',
            f'  {tool_name}:',
            f'    version: "{version}"',
            f'    path: "{tool_path}"',
            'parameters:',
        ]
        for k, v in params:
            lines.append(f'  {k}: {v}')

        out_file = Path(config.info_dir) / 'software_versions.yml'
        out_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        logger.info(f"版本信息已保存|Software versions saved: {out_file}")
    except Exception as e:
        # 版本记录为辅助信息, 失败不影响主流程|version record is auxiliary; never break the run
        logger.warning(f"生成software_versions.yml失败|Failed to write software_versions.yml: {e}")


def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    if config.method == 'fasttree':
        # FastTree: 经build_conda_command包装检查(与IQ-TREE一致)
        # |FastTree: check via build_conda_command (same as IQ-TREE)
        try:
            cmd = build_conda_command(config.fasttree_path, [])
            result = subprocess.run(
                cmd,
                capture_output=True, text=True, timeout=10
            )
            combined = (result.stdout or '') + (result.stderr or '')
            if 'FastTree' in combined:
                version_line = combined.strip().split('\n')[0]
                logger.info(f"FastTree 可用|FastTree available: {version_line}")
                return True
            else:
                raise RuntimeError("FastTree 未正确安装|FastTree not properly installed")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"FastTree 不可用|FastTree not available: {e}")
            raise RuntimeError(f"缺少依赖软件|Missing dependency: FastTree")

    elif config.method == 'iqtree':
        try:
            cmd = build_conda_command(config.iqtree_path, ['--version'])
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                version_line = (result.stdout or '').split('\n')[0]
                logger.info(f"IQ-TREE 可用|IQ-TREE available: {version_line}")
                return True
            else:
                raise RuntimeError("IQ-TREE 未正确安装|IQ-TREE not properly installed")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"IQ-TREE 不可用|IQ-TREE not available: {e}")
            raise RuntimeError(f"缺少依赖软件|Missing dependency: IQ-TREE")

    return True
