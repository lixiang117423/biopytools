"""全长转录本分析工具函数|Full-length RNA analysis utilities

日志(§2.3三handler分离)、conda命令执行(§13)、版本信息yml(§12.5)
|Logging (§2.3 3-handler), conda command execution (§13), versions yml (§12.5)
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

from ..common.conda_runner import build_conda_command


class RnaIsoLogger:
    """rna_iso日志管理器|rna_iso Logger Manager

    §2.3超算日志分离:INFO→stdout→.out,WARNING+→stderr→.err,DEBUG+→file
    """

    def __init__(self, log_dir, log_name: str = "rna_iso.log"):
        self.log_dir = Path(log_dir)
        self.log_file = self.log_dir / log_name
        self._setup_logging()

    def _setup_logging(self):
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 独立logger不污染root(§参考phylo_trim教训)|Independent logger, no root pollution
        logger = logging.getLogger("rna_iso")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        file_handler = logging.FileHandler(self.log_file, encoding='utf-8', mode='a')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        return self.logger


def run_command(cmd: List[str], logger, work_dir: Optional[str] = None,
                capture_output: bool = False) -> bool:
    """执行单条命令(自动conda包装,shell=False)|Run single command (auto conda wrap)

    Args:
        cmd: 命令列表,cmd[0]为工具完整路径|Command list, cmd[0] is tool full path
        logger: 日志器|Logger
        work_dir: 工作目录|Working directory
        capture_output: 是否捕获输出|Whether to capture output
    """
    if not cmd:
        logger.error("空命令|Empty command")
        return False
    wrapped = build_conda_command(cmd[0], cmd[1:])
    logger.info(f"命令|Command: {' '.join(wrapped)}")
    try:
        result = subprocess.run(
            wrapped, shell=False, cwd=work_dir,
            capture_output=capture_output, text=True, check=True
        )
        if capture_output and result.stdout:
            logger.debug(f"输出|Output: {result.stdout[:500]}")
        logger.info("命令执行成功|Command executed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command failed (rc={e.returncode}): {' '.join(wrapped)}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False
    except FileNotFoundError as e:
        logger.error(f"命令不存在|Command not found: {e}")
        return False


def generate_software_versions_yml(config, output_file: str, logger) -> None:
    """生成软件版本信息yml(§12.5)|Generate software versions yml (§12.5)"""
    from . import __version__

    tools = {
        'ccs': config.ccs_path,
        'isoseq3': config.isoseq3_path,
        'isoquant': config.isoquant_path,
        'samtools': config.samtools_path,
    }
    lines = [
        "pipeline:",
        f"  name: biopytools rna_iso",
        f"  version: '{__version__}'",
        "tools:",
    ]
    for name, path in tools.items():
        version = "unknown"
        try:
            # isoquant等Python入口冷启动import慢,30s宽限|Python entry tools import slowly on cold start
            result = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=30
            )
            out = (result.stdout.strip() or result.stderr.strip()).split("\n")[0]
            if out:
                version = out
        except Exception as e:
            logger.warning(f"版本探测失败|Version probe failed ({name}): {e}")
        lines.append(f"  {name}:")
        lines.append(f"    version: '{version}'")
        lines.append(f"    path: {path}")
    lines += [
        "parameters:",
        f"  engine: {config.engine}",
        f"  input_kind: {config.input_kind}",
        f"  data_type: {config.data_type}",
        f"  isoquant_data_type: {config.isoquant_data_type}",
        f"  threads: {config.threads}",
        f"  prefix: {config.prefix}",
    ]
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info(f"版本信息已写入|Versions written: {output_file}")
