"""MGA工具函数模块|MGA utility functions module"""

import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


def build_mga_command(mga_path: str, conda_env: str, reads: str,
                      output_dir: str, threads: int) -> list:
    """
    构建MGA命令(显式conda run,因MGA二进制不在env)|Build MGA command (explicit conda run)

    MGA二进制不在conda env中,但运行时需要env内的minimap2/samtools/python,
    因此不能用build_conda_command(get_conda_env返回None),必须显式包装。
    """
    args = ["--reads", reads, "--output", output_dir, "--threads", str(threads)]
    return ["conda", "run", "-n", conda_env, "--no-capture-output", mga_path] + args


class MGALogger:
    """MGA日志管理器|MGA logger manager(日志分离:stdout INFO / stderr WARNING+ / file DEBUG+)"""

    def __init__(self, logs_path, log_name: str = "mga.log"):
        self.logs_path = Path(str(logs_path))
        self.logs_path.mkdir(parents=True, exist_ok=True)
        self.log_file = self.logs_path / log_name
        self._setup_logging()

    def _setup_logging(self):
        self.logger = logging.getLogger("MGA")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        stdout_h = logging.StreamHandler(sys.stdout)
        stdout_h.setLevel(logging.INFO)
        stdout_h.setFormatter(formatter)
        stderr_h = logging.StreamHandler(sys.stderr)
        stderr_h.setLevel(logging.WARNING)
        stderr_h.setFormatter(formatter)
        file_h = logging.FileHandler(self.log_file, encoding="utf-8")
        file_h.setLevel(logging.DEBUG)
        file_h.setFormatter(formatter)

        self.logger.addHandler(stdout_h)
        self.logger.addHandler(stderr_h)
        self.logger.addHandler(file_h)

    def get_logger(self):
        """获取logger|Get logger"""
        return self.logger


def check_dependencies(config, logger) -> Optional[dict]:
    """
    检查conda env存在+MGA可执行+env内工具/python依赖|Check env, MGA, tools, python deps

    Returns:
        版本dict(成功)或None(失败)|version dict (ok) or None (fail)
    """
    env = config.conda_env
    logger.info(f"检查依赖|Checking dependencies (env={env})")
    versions = {}

    def _run(cmd, timeout=120):
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)

    try:
        # 1. MGA --help(验证二进制可执行+env可激活)|MGA --help
        r = _run(["conda", "run", "-n", env, "--no-capture-output", config.mga_path, "--help"])
        if r.returncode != 0:
            logger.error(f"MGA --help失败(检查env/mga_path)|MGA --help failed: {(r.stderr or '').strip()}")
            return None
        versions["MGA"] = {"version": "consensusLJA (no --version flag)", "path": config.mga_path}

        # 2. env内工具|minimap2 / samtools in env
        for tool in ["minimap2", "samtools"]:
            r = _run(["conda", "run", "-n", env, "--no-capture-output", tool, "--version"], timeout=60)
            out = (r.stdout or r.stderr or "").strip()
            if r.returncode != 0 or not out:
                logger.error(f"{tool}不可用(env={env})|{tool} not available: {out}")
                return None
            versions[tool] = {"version": out.split()[0]}

        # 3. python依赖|python deps (pysam/Bio/numpy)
        r = _run(["conda", "run", "-n", env, "python", "-c",
                  "import pysam,Bio,numpy;print(pysam.__version__,Bio.__version__,numpy.__version__)"],
                 timeout=60)
        if r.returncode != 0:
            logger.error(f"python依赖缺失(env={env})|python deps missing: {(r.stderr or '').strip()}")
            return None
        parts = r.stdout.strip().split()
        versions["python_deps"] = {
            "pysam": parts[0] if len(parts) > 0 else "unknown",
            "biopython": parts[1] if len(parts) > 1 else "unknown",
            "numpy": parts[2] if len(parts) > 2 else "unknown",
        }
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"依赖检查异常|Dependency check error: {e}")
        return None

    logger.info(f"依赖检查通过|Dependencies OK: {versions}")
    return versions


def write_software_versions(config, logger, start_time: datetime,
                            end_time: datetime, versions: dict) -> None:
    """生成00_pipeline_info/software_versions.yml|Generate software_versions.yml"""
    try:
        import yaml
    except ImportError:
        logger.warning("缺少pyyaml,跳过software_versions.yml|Missing pyyaml, skip")
        return

    info = {
        "pipeline": {"name": "biopytools mga", "version": "1.0.0"},
        "tools": versions,
        "parameters": {
            "reads": config.reads,
            "output_dir": config.output_dir,
            "threads": config.threads,
            "conda_env": config.conda_env,
            "dry_run": config.dry_run,
        },
        "execution": {
            "start_time": start_time.strftime("%Y-%m-%d %H:%M:%S"),
            "end_time": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "runtime_seconds": int((end_time - start_time).total_seconds()),
        },
    }
    out_file = config.info_path / "software_versions.yml"
    with open(out_file, "w") as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
    logger.info(f"版本信息已写入|Version info written: {out_file}")
