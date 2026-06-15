"""
Samplot Utilities
Samplot工具类
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from typing import List, Optional, Tuple


class SamplotLogger:
    """Samplot日志管理器|Samplot Logger Manager"""

    def __init__(self, output_dir: str, log_prefix: str = "samplot"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            log_prefix: 日志文件前缀|Log file prefix
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(output_dir, f"{log_prefix}_{timestamp}.log")

        self._setup_logging(log_file)

    def _setup_logging(self, log_file: str):
        """
        设置日志|Setup logging

        Args:
            log_file: 日志文件路径|Log file path
        """
        self.logger = logging.getLogger("Samplot")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO及以上|stdout handler - INFO and above
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # file handler - 所有级别|file handler - all levels
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中|Detect if command is in conda environment

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

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
    """
    构建conda run命令|Build conda run command

    Args:
        command: 命令路径|Command path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


class SamplotRunner:
    """Samplot命令运行器|Samplot Command Runner"""

    def __init__(self, logger):
        """
        初始化运行器|Initialize runner

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def run_command(
        self,
        cmd: List[str],
        description: str = "",
    ) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description

        Returns:
            (成功与否, stdout, stderr)|(success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        cmd_str = ' '.join(cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
            )
            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"退出码|Exit code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False, e.stdout or "", e.stderr or ""
        except FileNotFoundError as e:
            self.logger.error(f"命令不存在|Command not found: {e}")
            return False, "", str(e)

    def build_plot_command(self, config) -> List[str]:
        """
        构建samplot plot命令|Build samplot plot command

        Args:
            config: SamplotPlotConfig配置对象|SamplotPlotConfig object

        Returns:
            命令列表|Command list
        """
        args = ['plot']

        # 必需参数|Required parameters
        args.extend(['-b'] + config.bams)
        args.extend(['-c', config.chrom])
        args.extend(['-s', str(config.start)])
        args.extend(['-e', str(config.end)])

        # 可选参数|Optional parameters
        if config.sv_type:
            args.extend(['-t', config.sv_type])

        if config.output_file:
            args.extend(['-o', config.output_file])

        if config.output_dir and config.output_dir != ".":
            args.extend(['--output_dir', config.output_dir])

        if config.reference:
            args.extend(['-r', config.reference])

        if config.max_depth != 1:
            args.extend(['-d', str(config.max_depth)])

        if config.window is not None:
            args.extend(['-w', str(config.window)])

        if config.z != 4:
            args.extend(['-z', str(config.z)])

        if config.plot_height is not None:
            args.extend(['-H', str(config.plot_height)])

        if config.plot_width != 8:
            args.extend(['-W', str(config.plot_width)])

        if config.dpi != 300:
            args.extend(['--dpi', str(config.dpi)])

        if config.long_read != 1000:
            args.extend(['--long_read', str(config.long_read)])

        if config.coverage_only:
            args.append('--coverage_only')

        if config.same_yaxis_scales:
            args.append('--same_yaxis_scales')

        if config.titles:
            args.extend(['-n'] + config.titles)

        return build_conda_command(config.samplot_path, args)

    def build_vcf_command(self, config) -> List[str]:
        """
        构建samplot vcf命令|Build samplot vcf command

        Args:
            config: SamplotVcfConfig配置对象|SamplotVcfConfig object

        Returns:
            命令列表|Command list
        """
        args = ['vcf']

        # 必需参数|Required parameters
        args.extend(['--vcf', config.vcf])
        args.extend(['-b'] + config.bams)

        # 输出参数|Output parameters
        if config.output_dir != "samplot-out":
            args.extend(['-d', config.output_dir])

        if config.output_type != "png":
            args.extend(['-O', config.output_type])

        # 线程|Threads
        if config.threads != 1:
            args.extend(['-t', str(config.threads)])

        # 采样参数|Sampling parameters
        if config.downsample != 1:
            args.extend(['--downsample', str(config.downsample)])

        # 过滤参数|Filter parameters
        if config.min_bp != 20:
            args.extend(['--min_bp', str(config.min_bp)])

        if config.max_mb is not None:
            args.extend(['--max_mb', str(config.max_mb)])

        if config.sample_ids:
            args.extend(['--sample_ids'] + config.sample_ids)

        if config.plot_all:
            args.append('--plot_all')

        if config.min_call_rate is not None:
            args.extend(['--min_call_rate', str(config.min_call_rate)])

        if config.max_hets is not None:
            args.extend(['--max_hets', str(config.max_hets)])

        if config.min_entries != 6:
            args.extend(['--min_entries', str(config.min_entries)])

        if config.max_entries != 10:
            args.extend(['--max_entries', str(config.max_entries)])

        if config.gff3:
            args.extend(['--gff3', config.gff3])

        return build_conda_command(config.samplot_path, args)
