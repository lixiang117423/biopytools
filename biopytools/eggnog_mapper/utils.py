"""
eggnog-mapper 工具函数|eggnog-mapper utilities
conda 包装、日志管理器、命令运行器|conda wrapping, logger, runner
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令所属 conda 环境|Detect conda env of a command.

    优先从命令完整路径的 /envs/<name> 段提取|Extract from /envs/<name> in full path first.
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r"/envs/([^/]+)", cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get("CONDA_EXE")
    if conda_base:
        base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(base_dir, "envs")
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, "bin", command)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建 conda run 命令|Build conda run command.

    传完整路径(command 参数)以正确检测 env;必须含 --no-capture-output(规范§13.2.0)。
    Pass full path so env is detected; must include --no-capture-output.
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ["conda", "run", "-n", conda_env, "--no-capture-output", command] + args
    return [command] + args


class EggnogMapperLogger:
    """eggnog-mapper 日志管理器|logger manager (stdout INFO + stderr WARNING + file DEBUG)."""

    def __init__(self, output_dir, log_name: str = "emapper.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / "99_logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self._setup()

    def _setup(self):
        """配置日志 handler|Configure handlers."""
        if self.log_file.exists():
            self.log_file.unlink()
        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(self.log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上(规范§2.3.1, 超算 .err 捕获)|stderr handler - WARNING+
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        self.logger = logging.getLogger("eggnog_mapper")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.addHandler(file_handler)
        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.propagate = False

    def get_logger(self):
        """获取日志器|Get logger."""
        return self.logger


def build_emapper_args(config) -> List[str]:
    """
    构建 emapper 命令参数(纯函数,便于测试)|Build emapper args (pure).

    注意|Note: emapper 用 --output_dir(目录) + --output(前缀)。
    """
    args = [
        "-i", config.input_file,
        "--output_dir", str(Path(config.output_dir) / "01_emapper"),
        "--output", config.prefix,
        "--data_dir", config.data_dir,
        "--cpu", str(config.cpu),
        "-m", config.mode,
        "--itype", config.itype,
        "--sensmode", config.sensmode,
        "--seed_ortholog_evalue", str(config.seed_ortholog_evalue),
    ]
    if config.translate:
        args.append("--translate")
    if config.resume:
        args.append("--resume")
    if config.override:
        args.append("--override")
    return args


class EggnogMapperRunner:
    """eggnog-mapper 命令运行器|eggnog-mapper command runner."""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def check_engine(self) -> bool:
        """检查所选模式引擎可用(diamond 二进制等)|Check engine availability."""
        if self.config.mode != "diamond":
            return True
        try:
            cmd = build_conda_command("diamond", ["--version"])
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                self.logger.error(
                    "diamond不可用|diamond not available。安装|install: "
                    "conda install -n eggnog-mapper_v.2.1.15 -c bioconda diamond"
                )
                return False
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.logger.error(
                "diamond未安装|diamond not installed。安装|install: "
                "conda install -c bioconda diamond"
            )
            return False
        return True

    def run(self) -> bool:
        """运行注释流程|Run annotation pipeline."""
        self.logger.info("=" * 70)
        self.logger.info("开始eggnog-mapper功能注释|Starting eggnog-mapper annotation")
        self.logger.info("=" * 70)
        self.logger.info(f"输入|Input: {self.config.input_file}")
        self.logger.info(f"输出|Output: {self.config.output_dir}")
        self.logger.info(
            f"模式|Mode: {self.config.mode}  itype: {self.config.itype}  "
            f"cpu: {self.config.cpu}"
        )

        if self.config.translate and self.config.itype not in ("CDS", "genome", "metagenome"):
            self.logger.warning(
                "--translate对proteins无意义,忽略|translate ignored for proteins"
            )

        if not self.check_engine():
            return False

        annotations = str(
            Path(self.config.output_dir) / "01_emapper"
            / f"{self.config.prefix}.emapper.annotations"
        )

        # 断点续传:输出存在且非 override 则跳过 search|Checkpoint resume
        if Path(annotations).exists() and not self.config.override:
            self.logger.info(f"跳过已完成|Skipping (output exists): {annotations}")
        else:
            wrapped = build_conda_command(
                self.config.emapper_path, build_emapper_args(self.config)
            )
            self.logger.info(f"命令|Command: {' '.join(wrapped)}")
            try:
                result = subprocess.run(wrapped, shell=False)
            except FileNotFoundError as e:
                self.logger.error(f"emapper未找到|emapper not found: {e}")
                return False
            if result.returncode != 0:
                self.logger.error(
                    f"emapper执行失败|emapper failed (exit {result.returncode})"
                )
                return False
            self.logger.info("emapper运行完成|emapper finished")

        # 重排版(可跳过)|Reformat (skippable)
        if Path(annotations).exists() and not self.config.no_format:
            from .formatter import AnnotationsReformatter
            reformatter = AnnotationsReformatter(
                annotations,
                str(Path(self.config.output_dir) / "01_emapper"),
                self.logger,
            )
            reformatter.format()
        else:
            self.logger.info("跳过重排版|Skipping reformat")

        self._write_versions()
        return True

    def _write_versions(self):
        """写 software_versions.yml(失败不阻断)|Write versions yml."""
        versions_path = (
            Path(self.config.output_dir) / "00_pipeline_info" / "software_versions.yml"
        )
        version_str = self._detect_emapper_version()
        content = (
            "pipeline:\n"
            "  name: biopytools eggnog_mapper\n"
            "  version: '1.0.0'\n"
            "tools:\n"
            "  emapper:\n"
            f"    version: '{version_str}'\n"
            f"    path: '{self.config.emapper_path}'\n"
            "parameters:\n"
            f"  mode: {self.config.mode}\n"
            f"  itype: {self.config.itype}\n"
            f"  cpu: {self.config.cpu}\n"
            f"  sensmode: {self.config.sensmode}\n"
            f"  seed_ortholog_evalue: {self.config.seed_ortholog_evalue}\n"
        )
        try:
            versions_path.write_text(content, encoding="utf-8")
            self.logger.info(f"版本信息已保存|Versions saved: {versions_path}")
        except Exception as e:
            self.logger.warning(f"版本信息写入失败|Failed to write versions: {e}")

    def _detect_emapper_version(self) -> str:
        """探测 emapper 版本(失败返回 unknown)|Detect emapper version."""
        try:
            cmd = build_conda_command(self.config.emapper_path, ["--version"])
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            out = (result.stdout or result.stderr or "").strip().splitlines()
            return out[0] if out else "unknown"
        except Exception:
            return "unknown"
