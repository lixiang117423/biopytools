"""
工具函数模块|Utility Functions Module
"""

import subprocess
import logging
import os
from pathlib import Path
from typing import List, Optional, Union

from ..common.conda_runner import build_conda_command, check_tools

class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run

    def run(self, cmd: Union[List[str], str], description: str = "",
            check_output: bool = False, output_file: Optional[str] = None) -> Optional[str]:
        """执行命令|Execute command

        列表命令自动conda包装(shell=False), 字符串命令shell=True
        |List commands auto conda-wrapped (shell=False), strings shell=True
        """
        if description:
            self.logger.info(f"{description}")

        if isinstance(cmd, (list, tuple)):
            full_cmd = build_conda_command(cmd[0], cmd[1:])
            use_shell = False
            display = ' '.join(full_cmd)
        else:
            full_cmd = cmd
            use_shell = True
            display = cmd

        self.logger.info(f"命令|Command: {display}")

        if self.dry_run:
            self.logger.warning("干运行模式 - 跳过实际执行|Dry-run mode - Skipping actual execution")
            return None

        out_fh = None
        try:
            if output_file is not None:
                out_fh = open(output_file, 'w')
                result = subprocess.run(
                    full_cmd, shell=use_shell, stdout=out_fh,
                    stderr=subprocess.PIPE, text=True, check=True
                )
                out_fh.close()
                out_fh = None
                self.logger.info(f"执行成功|Execution successful")
                return True
            if check_output:
                result = subprocess.run(
                    full_cmd, shell=use_shell, capture_output=True,
                    text=True, check=True
                )
                return result.stdout.strip()
            else:
                result = subprocess.run(
                    full_cmd, shell=use_shell, check=True
                )
                self.logger.info(f"执行成功|Execution successful")
                return True

        except subprocess.CalledProcessError as e:
            if out_fh is not None:
                out_fh.close()
            self.logger.error(f"命令执行失败|Command execution failed")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if hasattr(e, 'stderr') and e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            raise

    def run_with_output(self, cmd: Union[List[str], str], description: str = "") -> str:
        """执行命令并返回输出|Execute command and return output"""
        return self.run(cmd, description, check_output=True)

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    missing_deps = check_tools([
        (config.samtools_path, "SAMtools", ["--version"]),
        (config.gatk_path, "GATK", ["--version"]),
    ], logger)

    for name in missing_deps:
        logger.error(f"  {name}: 未找到|Not found")

    if missing_deps:
        error_msg = f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    logger.info("所有依赖软件检查完成|All dependencies checked")
    return True

def check_file_exists(filepath: Path, logger) -> bool:
    """检查文件是否存在|Check if file exists"""
    if filepath.exists():
        logger.info(f"  文件已存在|File exists: {filepath.name}")
        return True
    return False
