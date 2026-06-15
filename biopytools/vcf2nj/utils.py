"""
VCF系统发育分析工具函数模块|VCF Phylogenetic Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
import shutil
import re
from pathlib import Path
from typing import Optional, List


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
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, command] + args
    else:
        return [command] + args

class PhyloLogger:
    """系统发育分析日志管理器|Phylogenetic Analysis Logger Manager"""

    def __init__(self, output_dir: Path, output_prefix: str, log_name: str = "phylo_analysis.log"):
        self.output_dir = Path(output_dir)
        self.output_prefix = output_prefix
        self.log_file = self.output_dir / f"{output_prefix}.log"
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 如果日志文件存在则删除|Remove log file if exists
        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令（自动检测conda环境）|Execute command (auto-detect conda environment)"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        # 提取命令名称用于检测conda环境|Extract command name for conda environment detection
        cmd_parts = cmd.strip().split()
        if cmd_parts:
            cmd_name = os.path.basename(cmd_parts[0])

            # 自动检测conda环境|Auto-detect conda environment
            conda_env = get_conda_env(cmd_name)

            if conda_env:
                # 使用conda run包装命令|Use conda run to wrap command
                full_cmd = f"conda run -n {conda_env} {cmd}"
                self.logger.debug(f"检测到conda环境|Detected conda environment: {conda_env}")
            else:
                # 直接执行命令|Execute command directly
                full_cmd = cmd
                self.logger.debug(f"未检测到conda环境，直接执行|No conda environment detected, executing directly")
        else:
            full_cmd = cmd

        try:
            result = subprocess.run(
                full_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = [
        (config.vcf2dis_path, "VCF2Dis")
    ]

    # 如果需要重根化，检查nw_reroot|Check nw_reroot if rerouting is needed
    if config.outgroup_list:
        dependencies.append((config.nw_reroot_path, "nw_reroot"))

    missing_deps = []

    for cmd, name in dependencies:
        try:
            # nw_reroot使用-h参数，其他工具使用--help|nw_reroot uses -h, others use --help
            help_param = "-h" if name == "nw_reroot" else "--help"
            cmd_name = os.path.basename(cmd)

            # 自动包装conda环境的命令|Auto-wrap conda environment commands
            wrapped_cmd = build_conda_command(cmd_name, [help_param])

            result = subprocess.run(wrapped_cmd,
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 or "usage" in result.stdout.lower() or "usage" in result.stderr.lower() or "synopsis" in result.stdout.lower():
                logger.info(f"{name} 可用|available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)

    # 检查Python依赖|Check Python dependencies
    python_deps = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("scikit-bio", "skbio")  # 注意：包名是scikit-bio，但导入名是skbio
    ]

    python_missing = []
    for dep_name, import_name in python_deps:
        try:
            __import__(import_name)
            logger.info(f"Python包 {dep_name} 可用|Python package {dep_name} available")
        except ImportError:
            python_missing.append(dep_name)

    if missing_deps:
        error_msg = f"缺少外部依赖软件|Missing external dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    if python_missing:
        error_msg = f"缺少Python依赖包|Missing Python dependencies: {', '.join(python_missing)}"
        logger.error(error_msg)

        # 提供安装建议|Provide installation suggestions
        logger.error("安装建议|Installation suggestions:")
        for dep in python_missing:
            if dep == "scikit-bio":
                logger.error(f"  pip install {dep}")
                logger.error(f"  或者|Or: conda install -c conda-forge {dep}")
                logger.error(f"  或者|Or: conda install -c bioconda {dep}")
            else:
                logger.error(f"  pip install {dep}")

        raise RuntimeError(error_msg)

    return True
