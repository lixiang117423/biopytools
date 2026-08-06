"""genome2sv 工具函数|genome2sv utilities

日志管理器、conda 环境包装、SV 统计纯函数。
|Logger, conda wrapping, and pure SV-stat helpers.
"""
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# SURVIVOR 统计的 SV 类型|SV types tracked for SURVIVOR stats
_VALID_SVTYPES = {"INS", "DEL", "INV", "DUP", "BND"}
_SUMMARY_COLS = ["sample", "total", "INS", "DEL", "INV", "DUP", "BND", "OTHER"]


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("genome2sv")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)   # WARNING+ → 超算 .err
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:                              # 全级别 → 文件
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在 conda 环境|Detect conda env of a command.

    先从 which 路径的 /envs/<name>/ 解析,否则搜索所有 conda 环境兜底。
    |Parse /envs/<name>/ from which path, else search all envs.
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        m = re.search(r"/envs/([^/]+)/", cmd_path)
        if m:
            return m.group(1)
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        base = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(base, "envs")
        if os.path.isdir(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, "bin", command)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建 conda run 命令(必带 --no-capture-output)|Build conda run command."""
    env = get_conda_env(command)
    if env:
        return ["conda", "run", "-n", env, "--no-capture-output", command] + args
    return [command] + args


def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查关键工具可用性(仅版本探测)|Check key tools via version probe."""
    logger.info("检查依赖工具|Checking dependencies")
    tools = {
        "minimap2": (config.minimap2_path, ["--version"]),
        "samtools": (config.samtools_path, ["--version"]),
        "svim-asm": (config.svim_asm_path, ["--version"]),
        "survivor": (config.survivor_path, []),
    }
    missing = []
    for name, (path, args) in tools.items():
        try:
            cmd = build_conda_command(path, args)
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if r.returncode == 0:
                logger.info(f"{name} 可用|available: {(r.stdout or r.stderr).strip()[:60]}")
            else:
                missing.append(name)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            missing.append(name)
    if missing:
        logger.error(f"缺失工具|Missing tools: {', '.join(missing)} (环境|env: sv_calling)")
        return False
    return True


def parse_svtype_stats(vcf_path: str) -> dict:
    """解析 VCF 统计 SVTYPE 分布|Parse VCF for SVTYPE distribution.

    Returns:
        dict: keys total/INS/DEL/INV/DUP/BND/OTHER;文件不存在返回全零。
        |dict with keys total/INS/DEL/INV/DUP/BND/OTHER; zeros if file missing.
    """
    counts = {k: 0 for k in _SUMMARY_COLS[1:]}
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                svtype = None
                for entry in fields[7].split(";"):
                    if entry.startswith("SVTYPE="):
                        svtype = entry.split("=", 1)[1]
                        break
                counts["total"] += 1
                if svtype in _VALID_SVTYPES:
                    counts[svtype] += 1
                else:
                    counts["OTHER"] += 1
    except FileNotFoundError:
        pass
    return counts


def build_survivor_input(sample_vcf_map: dict, output_file: str) -> List[str]:
    """写 SURVIVOR 输入文件(每行一个 VCF 绝对路径)|Write SURVIVOR input file.

    Args:
        sample_vcf_map: {sample_name: vcf_path} 仅含成功样本|successful samples only
        output_file: 输出文件路径|output path
    Returns:
        写入的绝对路径列表(按样本名排序)|list of absolute paths written (sorted)
    """
    paths = []
    with open(output_file, "w") as fh:
        for sample in sorted(sample_vcf_map):
            vcf = os.path.abspath(sample_vcf_map[sample])
            fh.write(vcf + "\n")
            paths.append(vcf)
    return paths


def format_sv_summary_tsv(rows: List[Tuple[str, dict]]) -> str:
    """格式化 SV 统计为 TSV|Format SV stats as TSV.

    Args:
        rows: [(name, counts_dict), ...]
    Returns:
        TSV 文本(表头 + 数据行)|TSV text with header and rows
    """
    lines = ["\t".join(_SUMMARY_COLS)]
    for name, counts in rows:
        cells = [name if c == "sample" else str(counts.get(c, 0))
                 for c in _SUMMARY_COLS]
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


def write_software_versions(config, logger: logging.Logger, output_path: str,
                            start_time=None) -> None:
    """生成 software_versions.yml|Generate software_versions.yml.

    探测 6 个工具版本 + 记录参数与运行时间,写入 00_pipeline_info/software_versions.yml。
    |Probe 6 tool versions, record parameters and runtime; write yml.
    """
    from datetime import datetime
    import yaml
    tools = {
        "minimap2": (config.minimap2_path, ["--version"]),
        "samtools": (config.samtools_path, ["--version"]),
        "svim-asm": (config.svim_asm_path, ["--version"]),
        "bcftools": (config.bcftools_path, ["--version"]),
        "bedtools": (config.bedtools_path, ["--version"]),
        "survivor": (config.survivor_path, []),
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
    # 关键参数(getattr 容错,适配测试桩)|key params (getattr-tolerant for stubs)
    param_keys = ["ref_sample", "preset", "svim_mode", "threads", "max_dist",
                  "min_support", "survivor_type", "survivor_strand", "est_dist",
                  "min_sv_length", "svim_min_sv_size"]
    info = {
        "pipeline": {"name": "biopytools genome2sv", "version": "1.0.0"},
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
