"""
ADMIXTURE分析工具函数模块|ADMIXTURE Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
import os
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(规范:模块内禁止复制实现)|conda wrapping via common layer (no local copies)
from ..common.conda_runner import build_conda_command, get_conda_env



class AdmixtureLogger:
    """ADMIXTURE分析日志管理器|ADMIXTURE Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "admixture_analysis.log",
                 log_level: str = "INFO", quiet: bool = False):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录路径|Output directory path
            log_name: 日志文件名|Log file name
            log_level: 日志级别 (DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level
            quiet: 静默模式（只输出ERROR）| Quiet mode (ERROR only)
        """
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.log_level = log_level
        self.quiet = quiet
        self.setup_logging()

    def setup_logging(self):
        """设置日志系统|Setup logging system"""
        # 确保输出目录存在|Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 使用命名日志器,隔离root避免污染全局(§2.3)|Named logger, isolate from root to avoid global pollution
        logger = logging.getLogger("biopytools.admixture")
        logger.handlers.clear()
        logger.propagate = False

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - 输出INFO及以下级别|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        if self.quiet:
            # 静默模式：只输出ERROR到stderr，stdout不输出|Quiet mode: only ERROR to stderr
            stdout_handler.setLevel(logging.CRITICAL + 1)  # 禁用stdout
        else:
            stdout_handler.setLevel(logging.INFO)
            stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - 输出WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # 设置日志级别|Set logger level
        logger.setLevel(logging.DEBUG)
        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path, dry_run: bool = False):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
            working_dir: 工作目录|Working directory
            dry_run: 模拟模式,只记录命令不执行|Dry run: log command without executing
        """
        self.logger = logger
        self.working_dir = working_dir
        self.dry_run = dry_run

    def run(self, cmd_list: list, description: str = "") -> str:
        """
        执行命令|Execute command

        Args:
            cmd_list: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description

        Returns:
            命令输出|Command output
        """
        if description:
            self.logger.info(f"开始|Starting: {description}")

        self.logger.info(f"执行命令|Executing command: {' '.join(cmd_list)}")

        # 模拟模式:记录命令后跳过执行|Dry run: log command, skip execution
        if self.dry_run:
            self.logger.info(f"模拟运行,跳过执行|Dry run, skipping execution: {description}")
            return ""

        try:
            result = subprocess.run(
                cmd_list,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                check=True,
                capture_output=True,
                text=True,
                cwd=self.working_dir
            )

            if result.stdout.strip():
                self.logger.info(f"命令输出|Command output:\n{result.stdout.strip()}")

            if result.stderr.strip():
                self.logger.warning(f"命令警告|Command warnings:\n{result.stderr.strip()}")

            self.logger.info(f"完成|Completed: {description}")
            return result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败|Command execution failed: {' '.join(cmd_list)}\n"
                f"   - 返回码|Return code: {e.returncode}\n"
                f"   - 标准输出|Stdout: {e.stdout.strip()}\n"
                f"   - 标准错误|Stderr: {e.stderr.strip()}"
            )
            raise


class SoftwareChecker:
    """软件环境检查器|Software Environment Checker"""

    def __init__(self, logger, config=None):
        """
        初始化软件检查器|Initialize software checker

        Args:
            logger: 日志对象|Logger object
            config: 配置对象(提供工具完整路径)|Config object (provides full tool paths)
        """
        self.logger = logger
        self.config = config

    def _probe(self, name: str, path: str, version_args: Optional[List[str]] = None) -> bool:
        """
        通过build_conda_command探活工具(非阻断)|Probe tool via build_conda_command (non-blocking)

        优先尝试运行版本参数;失败则回退到路径存在性检查。
        |Try running version args first; fallback to path existence on failure.
        """
        version_args = version_args if version_args is not None else ['--version']
        # 1) 尝试运行版本参数|Try running version flag
        try:
            cmd = build_conda_command(path, version_args)
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            if result.returncode == 0 or result.stdout or result.stderr:
                ver = (result.stdout or result.stderr or '').strip().split('\n')[0]
                self.logger.info(f"{name} 可用|{name} available: {ver or path}")
                return True
        except Exception as e:
            self.logger.debug(f"{name} 版本探活异常,回退路径检查|{name} version probe failed, fallback to path check: {e}")
        # 2) 回退:路径存在即视为可用|Fallback: path exists = available
        if os.path.exists(path) or shutil.which(path):
            self.logger.info(f"{name} 路径存在|{name} path exists: {path}")
            return True
        self.logger.warning(
            f"{name} 未就绪(路径|path={path});实际调用时若仍失败会报错|"
            f"{name} not ready; will error at call if still missing"
        )
        return False

    def check_software(self, software: str) -> bool:
        """
        检查软件是否可用(向后兼容)|Check software availability (backward compatible)

        Args:
            software: 工具完整路径或名称|Full tool path or name

        Returns:
            是否可用|Whether available
        """
        return self._probe(software, software)

    def check_dependencies(self) -> bool:
        """
        检查所有依赖软件(非阻断,仅warning)|Check all dependencies (non-blocking, warning only)

        必需工具缺失只发warning不退出——真实环境常通过conda env就绪,登录节点PATH检测
        会误判。实际命令失败时由对应步骤报错。
        |Missing required tools only warn, don't exit — tools are often ready via
        conda env on the compute node; login-node PATH check is unreliable. Real
        command failure is reported by the step that runs it.

        Returns:
            始终True(不阻断)|Always True (non-blocking)
        """
        self.logger.info("=== 软件环境检查|Software Environment Check ===")

        if self.config is not None:
            method = getattr(self.config, 'method', 'admixture')
            required = {
                'plink': (getattr(self.config, 'plink_path', 'plink'), ['--version']),
                'bcftools': (getattr(self.config, 'bcftools_path', 'bcftools'), ['--version']),
            }
            # 按method选核心工具|Select core tool by method
            if method == 'adamixture':
                required['adamixture'] = (getattr(self.config, 'adamixture_path', 'adamixture'), ['--help'])
            else:
                required['admixture'] = (getattr(self.config, 'admixture_path', 'admixture'), ['--version'])
            for name, (path, vargs) in required.items():
                self._probe(name, path, vargs)
        else:
            # 无config时回退到裸名探活|Fallback to bare-name probe when no config
            for software in ["plink", "admixture", "bcftools"]:
                self._probe(software, software)

        # 可选工具|Optional tools
        for software in ["Rscript"]:
            self._probe(software, software, ['--version'])

        self.logger.info("软件环境检查完成(缺失仅warning,不阻断)|Dependency check done (missing=warning, non-blocking)")
        return True


def check_file_format(file_path: str, expected_format: str) -> bool:
    """
    检查文件格式|Check file format

    Args:
        file_path: 文件路径|File path
        expected_format: 期望格式 (vcf/bed)|Expected format

    Returns:
        格式是否匹配|Whether format matches
    """
    if expected_format == "vcf":
        return file_path.endswith(('.vcf', '.vcf.gz'))
    elif expected_format == "bed":
        return file_path.endswith('.bed')
    return False
