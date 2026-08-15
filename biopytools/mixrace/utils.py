"""mixrace 工具函数|mixrace utilities

日志管理器、conda 包装、命令执行器、断点续传、版本信息。
组合自 genome2sv(ModuleLogger/get_conda_env/build_conda_command/write_software_versions)、
cim(CommandRunner SIGTERM 进程组管理)、rnaseq2vcf(CheckpointManager + 路径串优先的 env 探测)。
|Logger, conda wrapping, command runner, checkpoints, version info.
Assembled from genome2sv/cim/rnaseq2vcf sibling modules.
"""
import glob
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

LOG_FORMAT = "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"


class ModuleLogger:
    """模块日志管理器(三 handler: stdout/stderr/file)|Module logger (3 handlers)."""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("mixrace")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(LOG_FORMAT, LOG_DATEFMT)
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out|stdout
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)   # WARNING+ → 超算 .err|stderr
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:                              # 全级别 → 文件|file
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger."""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在 conda 环境|Detect conda env of a command.

    优先从路径串的 /envs/<name>/ 解析(确定性,不依赖 which),否则 shutil.which 兜底,
    最后搜索 CONDA_EXE 下所有 env。|Parse /envs/<name>/ from path string first
    (deterministic), else shutil.which, else search all envs under CONDA_EXE.
    """
    m = re.search(r"/envs/([^/]+)/", command)
    if m:
        return m.group(1)
    resolved = shutil.which(command)
    if resolved:
        m2 = re.search(r"/envs/([^/]+)/", resolved)
        if m2:
            return m2.group(1)
    # 不做"搜索所有 env"兜底:对绝对路径会误判(os.path.join 塌缩到该路径),
    # 且 mixrace 总是传完整工具路径,/envs/<name>/ 正则已能确定性命中。
    # |No "search all envs" fallback: it misfires on absolute paths (os.path.join
    # collapses to the path) and mixrace always passes full tool paths, so the
    # /envs/<name>/ regex already matches deterministically.
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建 conda run 命令(必带 --no-capture-output)|Build conda run command.

    传完整 tool 路径(勿 basename);非 conda 工具原样直调。|Pass full tool path
    (never basename); non-conda tools called directly.
    """
    env = get_conda_env(command)
    if env:
        return ["conda", "run", "-n", env, "--no-capture-output", command] + args
    return [command] + args


class CommandRunner:
    """命令执行器|Command runner.

    run(): str→shell=True(管道/重定向), list→shell=False;均捕获 stdout/stderr,
    返回 (success, stdout, stderr)。run_conda(): 传 tool 路径+args,自动 conda 包装。
    子进程以独立进程组启动,SIGTERM/超时/Ctrl-C 时整组杀(无孤儿 conda/samtools)。
    |run(): str→shell=True (pipes/redirects), list→shell=False; captures stdout/stderr,
    returns (success, stdout, stderr). run_conda(): tool path + args, auto conda-wrapped.
    Child in own process group; group-killed on SIGTERM/timeout/Ctrl-C (no orphans).
    """

    def __init__(self, logger: logging.Logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run
        self._child_proc: Optional[subprocess.Popen] = None
        try:
            signal.signal(signal.SIGTERM, self._on_terminate_signal)
        except (ValueError, OSError):
            pass  # 非主线程无法注册(测试环境)|cannot register outside main thread

    def _on_terminate_signal(self, signum, frame):
        self._kill_child_group()
        sys.exit(128 + signum)

    def _kill_child_group(self):
        proc = self._child_proc
        if proc is None:
            return
        try:
            pgid = os.getpgid(proc.pid)
        except ProcessLookupError:
            return
        self.logger.info(f"终止子进程组|Terminating child process group (pgid={pgid})")
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=5)
            return
        except subprocess.TimeoutExpired:
            pass
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    def run(self, command, description: str = "",
            timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """执行命令(str 或 list)|Execute command (str or list)."""
        is_str = isinstance(command, str)
        cmd_str = command if is_str else " ".join(shlex.quote(str(x)) for x in command)
        if description:
            self.logger.info(f"执行|Running: {description}")
        self.logger.info(f"命令|Command: {cmd_str}")
        if self.dry_run:
            return True, "", ""
        try:
            proc = subprocess.Popen(
                command, shell=is_str,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, start_new_session=True,
            )
        except Exception as e:
            self.logger.error(f"命令启动失败|Command launch error: {e}")
            return False, "", str(e)
        self._child_proc = proc
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command timed out ({timeout}s): {cmd_str}")
            self._kill_child_group()
            try:
                stdout, stderr = proc.communicate(timeout=10)
            except Exception:
                stdout, stderr = "", ""
            self._child_proc = None
            return False, stdout or "", stderr or "Timeout"
        except KeyboardInterrupt:
            self.logger.error("用户中断,终止子进程组|User interrupted, killing child group")
            self._kill_child_group()
            try:
                proc.communicate(timeout=10)
            except Exception:
                pass
            self._child_proc = None
            raise
        self._child_proc = None
        if proc.returncode != 0:
            self.logger.error(f"命令执行失败|Command failed (exit {proc.returncode})")
            if stderr:
                for line in stderr.strip().splitlines()[-20:]:
                    self.logger.error(f"  {line}")
            return False, stdout, stderr
        return True, stdout, stderr

    def run_conda(self, tool_path: str, args: List[str], description: str = "",
                  timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """conda 包装执行单工具|Run a single tool via conda wrapping."""
        cmd = build_conda_command(tool_path, args)
        return self.run(cmd, description, timeout)


class CheckpointManager:
    """断点续传管理(基于 .done 文件存在性)|Checkpoint manager (file-existence based)."""

    def __init__(self, checkpoint_dir: str, logger: logging.Logger):
        self.checkpoint_dir = checkpoint_dir
        self.logger = logger
        os.makedirs(checkpoint_dir, exist_ok=True)

    def _path(self, step: str) -> str:
        return os.path.join(self.checkpoint_dir, f"{step}.done")

    def exists(self, step: str) -> bool:
        return os.path.exists(self._path(step))

    def create(self, step: str):
        try:
            open(self._path(step), "w").close()
        except Exception as e:
            self.logger.warning(f"写断点失败|Failed to write checkpoint {step}: {e}")

    def remove(self, step: str):
        if self.exists(step):
            os.remove(self._path(step))

    def list_completed(self) -> List[str]:
        return [os.path.basename(p)[:-5]
                for p in glob.glob(os.path.join(self.checkpoint_dir, "*.done"))]


def format_number(num) -> str:
    """大数字(>=1M)用 M 单位保留2位小数|Format big numbers (>=1M) as M with 2 decimals."""
    try:
        n = int(num)
    except (TypeError, ValueError):
        return str(num)
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f}M"
    return str(n)


def write_software_versions(config, logger: logging.Logger, output_path: str,
                            start_time: Optional[datetime] = None) -> None:
    """生成 software_versions.yml|Generate software_versions.yml.

    探测 mixrace 工具链版本(samtools/bcftools/bedtools/Rscript)+ 记录参数与运行时间。
    |Probe mixrace tool versions + record parameters/runtime.
    """
    import yaml
    tools = {
        "samtools": (config.samtools_path, ["--version"]),
        "bcftools": (config.bcftools_path, ["--version"]),
        "bedtools": (config.bedtools_path, ["--version"]),
        "Rscript": (config.rscript_path, ["--version"]),
    }
    versions = {}
    for name, (path, args) in tools.items():
        try:
            cmd = build_conda_command(path, args)
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            raw = (r.stdout or r.stderr).strip()
            ver = raw.splitlines()[0] if raw else "unknown"
            versions[name] = {"version": ver, "path": path}
        except Exception as e:
            logger.warning(f"版本探测失败|Version probe failed [{name}]: {e}")
            versions[name] = {"version": "unknown", "path": path}
    param_keys = ["threads", "kmer_size", "read_length", "min_qual", "min_dp",
                  "min_alt_reads", "freebayes_min_coverage",
                  "freebayes_min_alternate_fraction", "het_pure", "het_suspicious",
                  "het_impure", "min_depth", "repeat_bed", "skip_tree"]
    info = {
        "pipeline": {"name": "biopytools mixrace", "version": "0.1.0"},
        "tools": versions,
        "parameters": {k: getattr(config, k, None) for k in param_keys},
    }
    end_time = datetime.now()
    if start_time is not None:
        info["execution"] = {
            "start_time": start_time.strftime("%Y-%m-%d %H:%M:%S"),
            "end_time": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "runtime_seconds": int((end_time - start_time).total_seconds()),
        }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
