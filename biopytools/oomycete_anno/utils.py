"""疫霉菌注释工具函数|Oomycete annotation utilities

conda 包装 / 日志管理器 / 命令运行器
|conda wrapping, logger manager, command runner
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# ============================================================
# conda 环境检测与命令构建|Conda env detection & command building
# ============================================================

def get_conda_env(command: str) -> Optional[str]:
    """检测命令所属 conda 环境|Detect conda env of a command.

    优先从命令完整路径的 /envs/<name> 段提取; 找不到则扫描所有 conda 环境。
    |Extract from /envs/<name> in full path first; else scan all envs.
    """
    # 必须传完整路径(规范§13.6), 命令名无法可靠检测 env|Full path required (spec §13.6)
    cmd_path = shutil.which(command) or command
    if cmd_path:
        match = re.search(r"/envs/([^/]+)", cmd_path)
        if match:
            return match.group(1)

    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(base_dir, "envs")
        if os.path.isdir(envs_dir):
            cmd_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, "bin", cmd_name)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建 conda run 命令(传完整路径)|Build conda run command (full path).

    传完整路径(command)以正确检测 env; 必须含 --no-capture-output(规范§13.2.0)。
    |Pass full path so env is detected; must include --no-capture-output.
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ["conda", "run", "-n", conda_env, "--no-capture-output", command] + args
    return [command] + args


def build_genemark_command(
    gmes_petap_path: str, perl_env: str, args: List[str]
) -> List[str]:
    """构建 GeneMark-ES/ET/EP+ 命令(特殊)|Build GeneMark command (special).

    GeneMark 安装路径无 /envs/, get_conda_env 检测不到; 且其 .pl shebang 为
    #!/usr/bin/env perl, 必须跑在带 7 个 CPAN 模块的环境里(已验证 braker_v.3.0.8)。
    所以显式用 conda run -n <perl_env> perl <gmes_petap_path> 调用。
    |GeneMark lives outside /envs/ so env auto-detection fails; its shebang is
    env-perl, so it MUST run inside an env with the 7 CPAN modules (braker_v.3.0.8).
    Hence explicit: conda run -n <perl_env> perl <gmes_petap_path>.
    """
    return (
        ["conda", "run", "-n", perl_env, "--no-capture-output", "perl", gmes_petap_path]
        + args
    )


# ============================================================
# 日志管理器|Logger manager (stdout INFO + stderr WARNING + file DEBUG)
# ============================================================

class OomyceteAnnoLogger:
    """疫霉菌注释日志管理器|Logger manager.

    三路 handler 分离(规范§2.3): file=DEBUG, stdout=INFO, stderr=WARNING+。
    超算上 INFO -> .out, WARNING/ERROR -> .err。
    |Three-way split (spec §2.3): file DEBUG, stdout INFO, stderr WARNING+.
    On the cluster INFO -> .out, WARNING/ERROR -> .err.
    """

    def __init__(self, log_dir, log_name: str = "oomycete_anno.log"):
        self.log_file = Path(log_dir) / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self._setup()

    def _setup(self):
        """配置 handler|Configure handlers."""
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

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        self.logger = logging.getLogger("oomycete_anno")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.addHandler(file_handler)
        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.propagate = False

    def get_logger(self):
        """获取日志器|Get logger."""
        return self.logger


# ============================================================
# 命令执行器|Command runner
# ============================================================

class CommandRunner:
    """命令执行器|Command runner.

    重型步骤(hisat2/augustus/genemark 等)默认不捕获输出, 让其流式写到超算
    .out/.err, 避免内存膨胀且保留进度。仅 --version 等小输出用 run_capture。
    |Heavy steps stream (no capture) to cluster .out/.err to avoid OOM and keep
    progress; only small outputs (e.g. --version) use run_capture.
    """

    def __init__(self, logger: logging.Logger, working_dir: str = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(
        self,
        cmd: List[str],
        description: str = "",
        cwd: str = None,
        env: dict = None,
        stdout=None,
        stdin=None,
        timeout: int = 86400,
    ) -> bool:
        """执行命令(流式, 不捕获)|Execute command (streaming, no capture).

        Args:
            cmd: 命令列表(build_conda_command 构建)|command list
            description: 步骤描述|step description
            cwd: 工作目录|working directory
            env: 环境变量(如 AUGUSTUS_CONFIG_PATH)|env vars
            stdout: 输出文件句柄(augustus 重定向用)|output file handle
            stdin: 输入文件句柄(gtf2gff.pl 喂 GTF 用)|input file handle
            timeout: 超时秒|timeout seconds

        Returns:
            bool: 是否成功|success
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        # 规范§2.2.1: 执行前必须记录完整命令到 INFO|log full command at INFO (spec §2.2.1)
        self.logger.info(f"命令|Command: {' '.join(str(c) for c in cmd)}")

        try:
            subprocess.run(
                cmd,
                shell=False,
                check=True,
                cwd=cwd or self.working_dir,
                env=env,
                stdout=stdout,
                stdin=stdin,
                timeout=timeout,
            )
            if description:
                self.logger.info(f"{description} 完成|completed")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败|Command failed (exit {e.returncode}): {description}"
            )
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(
                f"命令超时|Command timed out ({timeout}s): {description}"
            )
            return False
        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False

    def run_capture(
        self, cmd: List[str], description: str = "", timeout: int = 60
    ) -> Tuple[bool, str, str]:
        """捕获输出执行(仅小输出如 --version)|Run with capture (small outputs only).

        Returns:
            (success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(str(c) for c in cmd)}")
        try:
            r = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=self.working_dir,
            )
            return (r.returncode == 0), r.stdout, r.stderr
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令超时|Command timed out: {description}")
            return False, "", "timeout"
        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False, "", f"Command not found: {cmd[0]}"


# ============================================================
# 辅助函数|Helpers
# ============================================================

def format_number(num: int) -> str:
    """格式化大数字(M/K)|Format large numbers (M/K)."""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    if num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def safe_ascii_path(path: str, cache_dir: str = None) -> str:
    """非 ASCII 路径 -> ASCII 软链; 已是 ASCII 则原样返回
    |Non-ASCII path -> ASCII symlink; return as-is if already ASCII.

    GeneMark probuild / Augustus 等 C 二进制会把 UTF-8 中文路径字节解成乱码导致
    fopen 失败。用 ~/tmp(或 cache_dir) 下的 ASCII 名软链指向真实文件,C 二进制只看到
    ASCII 路径即可正常读写。
    |C binaries (GeneMark probuild, Augustus) mangle UTF-8 CJK path bytes, breaking
    fopen. Symlink the real file to an ASCII name in ~/tmp (or cache_dir) so the C
    binary only sees an ASCII path.
    """
    if not path:
        return path
    try:
        path.encode("ascii")
        return path
    except UnicodeEncodeError:
        pass
    import hashlib
    base = Path(cache_dir) if cache_dir else Path.home() / "tmp"
    base.mkdir(parents=True, exist_ok=True)
    target = os.path.abspath(path)
    link = base / f"oo_{hashlib.md5(target.encode('utf-8')).hexdigest()[:12]}"
    if link.is_symlink():
        if os.readlink(str(link)) != target:
            link.unlink()
            os.symlink(target, str(link))
    elif link.exists():
        shutil.rmtree(str(link))
        os.symlink(target, str(link))
    else:
        os.symlink(target, str(link))
    return str(link)


def find_rnaseq_pairs(
    rnaseq_dirs: List[str],
    read1_pattern: str,
    read2_pattern: str,
    logger: logging.Logger,
) -> List[Tuple[str, str]]:
    """在多个 RNA-seq 目录中查找 R1/R2 配对|Find R1/R2 pairs across dirs.

    Returns:
        [(r1_path, r2_path), ...]
    """
    pairs = []
    for d in rnaseq_dirs:
        d_path = Path(d)
        if not d_path.is_dir():
            logger.warning(f"RNA-seq 目录不存在,跳过|RNA-seq dir not found, skipping: {d}")
            continue
        r1_files = sorted(d_path.glob(f"*{read1_pattern}"))
        for r1 in r1_files:
            # R2 = R1 名去掉 read1_pattern 再接 read2_pattern|R2 from R1 stem swap
            stem = r1.name[: -len(read1_pattern)] if read1_pattern else r1.stem
            r2 = d_path / f"{stem}{read2_pattern}"
            if r2.exists():
                pairs.append((str(r1), str(r2)))
            else:
                logger.warning(
                    f"R2 不匹配,跳过该对|R2 missing, skipping pair: {r1.name}"
                )
    return pairs


def miniprot_gff_to_protein_hints(miniprot_gff: str, out_hints: str) -> int:
    """miniprot GFF3 -> Augustus protein hints(CDSpart, source=P)
    |Convert miniprot GFF3 to Augustus protein hints (CDSpart, source=P).

    miniprot --gff 输出 mRNA/CDS/stop_codon 三类特征; 取 CDS 行(蛋白编码外显子的
    基因组坐标)作为 CDSpart hints(蛋白编码区正向证据)。
    |miniprot --gff emits mRNA/CDS/stop_codon; use CDS lines (coding-exon genome
    coords) as CDSpart hints (positive coding evidence).

    Returns:
        写入的 hints 行数|number of hint lines written
    """
    n = 0
    with open(miniprot_gff) as fin, open(out_hints, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            if cols[2] != "CDS":
                continue
            seq, _src, _ftype, start, end, _score, strand = cols[:7]
            fout.write(
                f"{seq}\tP\tCDSpart\t{start}\t{end}\t0\t{strand}\t.\tsource=P\n"
            )
            n += 1
    return n
