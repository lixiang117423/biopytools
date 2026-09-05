"""easyhap 核心编排|easyhap orchestration"""
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

from ..common.conda_runner import build_conda_command
from .config import EasyHapConfig
from .utils import (ModuleLogger, parse_region_string, read_region_file,
                    region_label)


@dataclass
class RegionResult:
    """单区域运行结果|Single region result"""
    label: str
    ok: bool
    message: str = ""


class EasyHapCalculator:
    """逐区域调用上游 easyhap analyze|Per-region invocation of upstream easyhap analyze"""

    def __init__(self, config: EasyHapConfig):
        self.config = config
        self.logger = ModuleLogger(config.log_file, config.log_level).get_logger()

    # ---------- 区域解析|region loading ----------

    def _load_regions(self) -> List[Tuple[str, int, int]]:
        if self.config.region:
            return [parse_region_string(self.config.region)]
        return read_region_file(self.config.region_file)

    # ---------- 命令构建|command building ----------

    def _build_command(self, chrom: str, start: int, end: int) -> List[str]:
        args = ["analyze",
                "--vcf", self.config.input_vcf,
                "--group", self.config.group_file,
                "--region", f"{chrom}:{start}-{end}",
                "--outdir", str(self.config.haplotypes_dir),
                "--mode", self.config.mode,
                "--hetero-policy", self.config.hetero_policy,
                "--cluster-threshold", str(self.config.cluster_threshold),
                "--vcf-backend", self.config.vcf_backend,
                "--fisher-adjust", self.config.fisher_adjust,
                "--plot-hap-level", self.config.plot_hap_level,
                "--plot-min-count", str(self.config.plot_min_count),
                "--plot-format", self.config.plot_format]
        if self.config.traits:
            args += ["--traits", self.config.traits]
        if self.config.trait_cols:
            args += ["--trait-cols", self.config.trait_cols]
        if self.config.fisher_groups:
            args += ["--fisher-groups", self.config.fisher_groups]
        if self.config.fisher_alpha is not None:
            args += ["--fisher-alpha", str(self.config.fisher_alpha)]
        if self.config.no_processed:
            args.append("--no-processed")
        if self.config.plot:
            args.append("--plot")
        if self.config.gff:
            args += ["--gff", self.config.gff]
        return build_conda_command(self.config.easyhap_path, args)

    # ---------- 断点续传|checkpointing ----------

    def _done_marker(self, label: str) -> Path:
        # 包装层自管 done-marker, 不用 HapSummary.tsv:
        # --plot 出图阶段失败时 HapSummary 可能已写出, 以它判断会永久跳过缺图的区域
        # |Wrapper-owned marker: HapSummary may exist even when plotting failed
        return self.config.haplotypes_dir / f"{label}_easyhap.done"

    def _is_done(self, label: str) -> bool:
        return self._done_marker(label).exists()

    # ---------- 单区域执行|single region ----------

    def _run_region(self, chrom: str, start: int, end: int) -> RegionResult:
        label = region_label(chrom, start, end)
        if self._is_done(label) and not self.config.force:
            self.logger.info(f"跳过已完成区域|Skipping completed region: {label}")
            return RegionResult(label, True, "skipped")
        cmd = self._build_command(chrom, start, end)
        self.logger.info(f"分析区域|Analyzing region: {chrom}:{start}-{end}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, check=False,
                                    capture_output=True, text=True)
        except (FileNotFoundError, OSError) as e:
            self.logger.error(f"区域执行异常|Region execution error: {label}: {e}")
            return RegionResult(label, False, str(e))
        if result.returncode != 0:
            tail = "\n".join((result.stderr or result.stdout or "")
                             .strip().splitlines()[-5:])
            self.logger.error(
                f"区域分析失败|Region failed: {label} (exit {result.returncode})")
            if tail:
                self.logger.error(f"错误尾部|Error tail:\n{tail}")
            return RegionResult(label, False, f"exit {result.returncode}")
        self._done_marker(label).write_text(
            f"region={chrom}:{start}-{end}\n"
            f"time={datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.logger.info(f"区域完成|Region completed: {label}")
        return RegionResult(label, True, "")

    # ---------- 工具检查|tool check ----------

    def _check_tool(self) -> bool:
        if Path(self.config.easyhap_path).exists():
            return True
        self.logger.error(
            f"easyhap 不存在|easyhap not found: {self.config.easyhap_path}")
        self.logger.error(
            "安装命令|Install: ~/miniforge3/bin/conda run -n pop python -m pip "
            "install --no-deps ~/software/EasyHap")
        return False

    # ---------- 版本记录|versions ----------

    def _probe_tool_version(self) -> str:
        """easyhap --version(经 conda 包装)|Probe easyhap version via conda wrapper"""
        cmd = build_conda_command(self.config.easyhap_path, ["--version"])
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            r = subprocess.run(cmd, shell=False, check=False,
                               capture_output=True, text=True, timeout=60)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return "unknown"
        m = re.search(r"EasyHap\s+([\d.]+)", (r.stdout or "") + (r.stderr or ""))
        return m.group(1) if m else "unknown"

    def _probe_python_version(self) -> str:
        """pop 环境 python 版本|pop env python version"""
        cmd = build_conda_command(
            str(Path(self.config.easyhap_path).parent / "python"), ["--version"])
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            r = subprocess.run(cmd, shell=False, check=False,
                               capture_output=True, text=True, timeout=60)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return "unknown"
        return (r.stdout or r.stderr or "").strip() or "unknown"

    def write_software_versions(self, start_time=None) -> None:
        """生成 software_versions.yml|Generate software_versions.yml"""
        import yaml
        info = {
            "pipeline": {"name": "biopytools easyhap", "version": "1.0.0"},
            "tools": {
                "easyhap": {"version": self._probe_tool_version(),
                            "path": self.config.easyhap_path},
                "python": {"version": self._probe_python_version()},
            },
            "parameters": {k: getattr(self.config, k, None) for k in
                           ("input_vcf", "group_file", "region", "region_file",
                            "output_dir", "mode", "hetero_policy",
                            "cluster_threshold", "vcf_backend", "fisher_groups",
                            "fisher_alpha", "fisher_adjust", "plot",
                            "plot_hap_level", "plot_min_count", "plot_format")},
        }
        end_time = datetime.now()
        if start_time is not None:
            info["execution"] = {
                "start_time": start_time.strftime("%Y-%m-%d %H:%M:%S"),
                "end_time": end_time.strftime("%Y-%m-%d %H:%M:%S"),
                "runtime_seconds": int((end_time - start_time).total_seconds()),
            }
        out = str(self.config.info_dir / "software_versions.yml")
        with open(out, "w") as fh:
            yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本信息已生成|software versions written: {out}")

    # ---------- 警告与主流程|warnings & main run ----------

    def _emit_warnings(self) -> None:
        if self.config.gff and not self.config.plot:
            self.logger.warning(
                "--gff 仅在 --plot 开启时使用|--gff only takes effect with --plot")
        if self.config.mode == "hybrid" and self.config.hetero_policy != "slash":
            self.logger.warning(
                "--hetero-policy 仅 inbred 模式生效|"
                "--hetero-policy only applies in inbred mode")
        if self.config.fisher_alpha is not None and not self.config.fisher_groups:
            self.logger.warning(
                "--fisher-alpha 需 --fisher-groups 才生效|"
                "--fisher-alpha requires --fisher-groups")
        if self.config.fisher_adjust != "none" and not self.config.fisher_groups:
            self.logger.warning(
                "--fisher-adjust 需 --fisher-groups 才生效|"
                "--fisher-adjust requires --fisher-groups")

    def run(self) -> int:
        """端到端运行|Run end-to-end. Returns exit code 0/1."""
        start_time = datetime.now()
        self._emit_warnings()
        if not self._check_tool():
            return 1
        try:
            regions = self._load_regions()
        except (ValueError, OSError) as e:
            self.logger.error(f"区域解析失败|Region parsing failed: {e}")
            return 1
        self.logger.info(f"待分析区域数|Regions to analyze: {len(regions)}")
        results = [self._run_region(c, s, e) for c, s, e in regions]
        ok = sum(1 for r in results if r.ok)
        failed = [r for r in results if not r.ok]
        self.logger.info(f"区域统计|Region summary: 成功|success {ok}/{len(regions)}")
        for r in failed:
            self.logger.error(f"失败区域|Failed region: {r.label}: {r.message}")
        self.write_software_versions(start_time=start_time)
        if failed:
            self.logger.error(
                "存在失败区域(成功产出已保留)|Some regions failed (outputs kept)")
            return 1
        self.logger.info("easyhap 流程完成|easyhap completed")
        return 0
