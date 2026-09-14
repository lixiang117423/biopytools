"""
VCF系统发育分析工具函数模块|VCF Phylogenetic Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command, conda_run_prefix


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

            # 自动检测conda环境(同源conda绝对路径+run -p前缀,§13)|Auto-detect conda env
            prefix = conda_run_prefix(cmd_name)

            if prefix:
                # 使用conda run前缀包装命令|Wrap command with the conda run prefix
                full_cmd = f"{prefix} {cmd}"
                self.logger.debug(f"检测到conda环境|Detected conda environment: {prefix}")
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
