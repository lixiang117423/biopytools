"""genome2tree 流程编排器|genome2tree orchestrator"""
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import List

from .utils import (decompress_gz_samples, parse_samples_map, write_input_tsv,
                    write_samples_map_tsv)


class Genome2TreeCalculator:
    """目录扫描→waster 建树→(可选)枝长 编排器|Scan → waster → (optional) branch lengths"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    # ---------- 命令构造|command building ----------

    def _root_args(self) -> List[str]:
        """--root 外群参数|Return --root args if set"""
        return ["--root", self.config.root] if self.config.root else []

    def _mapping_args(self) -> List[str]:
        """-a 映射参数(有 --samples-map 时)|Return -a args if samples_map set"""
        return ["-a", self.config.samples_map_tsv] if self.config.samples_map else []

    def waster_command(self) -> List[str]:
        """02 步 waster 命令|Build waster command"""
        cmd = [self.config.waster_path,
               "-i", self.config.input_tsv,
               "-o", self.config.species_tree,
               "-t", str(self.config.threads)]
        return cmd + self._root_args() + self._mapping_args()

    def branchlength_command(self) -> List[str]:
        """03 步 waster_branchlength 命令(-C -c 固定拓扑打分)|Build branch-length command"""
        cmd = [self.config.waster_branchlength_path,
               "-i", self.config.input_tsv,
               "-C", "-c", self.config.species_tree,
               "-o", self.config.bl_tree,
               "-t", str(self.config.threads)]
        return cmd + self._root_args() + self._mapping_args()

    # ---------- 命令执行|command execution ----------

    def _run_tool(self, cmd: List[str], log_path: str) -> bool:
        """执行外部命令,stderr 落盘|Run tool; stderr to log_path. False on failure."""
        self.logger.info(f"命令|Command: {' '.join(str(c) for c in cmd)}")
        try:
            with open(log_path, "w") as fh:
                result = subprocess.run(cmd, shell=False, check=False,
                                        stdout=subprocess.PIPE, stderr=fh, text=True)
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {e}")
            return False
        if result.returncode != 0:
            self.logger.error(
                f"命令失败|Command failed (rc={result.returncode}): "
                f"{self._log_tail(log_path)}")
            return False
        return True

    @staticmethod
    def _log_tail(log_path: str, limit: int = 500) -> str:
        """读日志尾部用于报错|Read log tail for error messages"""
        try:
            with open(log_path) as fh:
                return fh.read()[-limit:].strip()
        except OSError:
            return "(无日志|no log)"

    # ---------- 步骤 01:输入准备|input prep ----------

    def prepare_inputs(self) -> bool:
        """解压 gz、写 input.tsv 与 samples_map.tsv|Decompress & write waster inputs"""
        cfg = self.config
        samples = decompress_gz_samples(cfg.raw_samples,
                                        str(cfg.uncompressed_dir), self.logger)
        write_input_tsv(samples, cfg.input_tsv)
        mapping = parse_samples_map(cfg.samples_map) if cfg.samples_map else {}
        if mapping:
            species_of = {s.stem: mapping.get(s.stem, s.stem) for s in samples}
            write_samples_map_tsv(species_of, cfg.samples_map_tsv)
            self.logger.info(f"样本映射已生成|samples map written: {cfg.samples_map_tsv}")
        species = [mapping.get(s.stem, s.stem) for s in samples]
        self.logger.info(
            f"输入清单已生成|input list written: {cfg.input_tsv} "
            f"(样本|samples: {len(samples)}, 物种|species: {len(set(species))})")
        if cfg.ignored_files:
            self.logger.warning(
                f"忽略非序列文件|Ignored non-sequence files: {cfg.ignored_files}")
        if len(set(species)) < 4:
            self.logger.warning(
                "物种数<4,四聚体方法支持度可能不稳|<4 species; quartet-based "
                "support may be unstable")
        return True

    # ---------- 步骤 02/03:建树/枝长|tree & branch lengths ----------

    def run_waster(self) -> bool:
        """waster 建树(断点续传)|Run waster (checkpointed)"""
        if Path(self.config.species_tree).exists():
            self.logger.info("跳过已完成步骤|Skipping completed step: waster")
            return True
        self.logger.info("开始步骤|Starting step: waster 物种树|species tree")
        ok = self._run_tool(self.waster_command(), self.config.waster_log)
        if ok:
            self.logger.info(
                f"物种树已生成|species tree written: {self.config.species_tree}")
        return ok

    def run_branch_length(self) -> bool:
        """枝长计算(未启用为 no-op,断点续传)|Branch lengths (no-op unless enabled)"""
        if not self.config.branch_length:
            return True
        if Path(self.config.bl_tree).exists():
            self.logger.info("跳过已完成步骤|Skipping completed step: branch length")
            return True
        self.logger.info(
            "开始步骤|Starting step: waster_branchlength 枝长|branch lengths")
        ok = self._run_tool(self.branchlength_command(), self.config.bl_log)
        if ok:
            self.logger.info(
                f"枝长树已生成|branch-length tree written: {self.config.bl_tree}")
        return ok

    # ---------- 版本记录|versions ----------

    def _probe_waster_version(self, tool_path: str) -> str:
        """从 --help 输出抓 Version 行(实测走 stderr)|Probe Version: line from --help."""
        self.logger.info(f"命令|Command: {tool_path} --help")
        try:
            r = subprocess.run([tool_path, "--help"], shell=False, check=False,
                               capture_output=True, text=True, timeout=60)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return "unknown"
        m = re.search(r"Version:\s*(\S+)", r.stderr or "")
        return m.group(1) if m else "unknown"

    def _aster_git_commit(self, tool_path: str) -> str:
        """取 ASTER 仓库 commit(失败 unknown)|ASTER repo commit or unknown"""
        repo = str(Path(tool_path).resolve().parent.parent)
        self.logger.info(f"命令|Command: git -C {repo} rev-parse --short HEAD")
        try:
            r = subprocess.run(["git", "-C", repo, "rev-parse", "--short", "HEAD"],
                               shell=False, check=False,
                               capture_output=True, text=True, timeout=30)
            return r.stdout.strip() if r.returncode == 0 else "unknown"
        except Exception:
            return "unknown"

    def write_software_versions(self, start_time=None) -> None:
        """生成 software_versions.yml|Generate software_versions.yml"""
        import yaml
        versions = {
            "waster": {"version": self._probe_waster_version(self.config.waster_path),
                       "path": self.config.waster_path},
        }
        if self.config.branch_length:
            versions["waster_branchlength"] = {
                "version": self._probe_waster_version(
                    self.config.waster_branchlength_path),
                "path": self.config.waster_branchlength_path}
        info = {
            "pipeline": {"name": "biopytools genome2tree", "version": "1.0.0"},
            "aster_repository": {
                "git_commit": self._aster_git_commit(self.config.waster_path)},
            "tools": versions,
            "parameters": {k: getattr(self.config, k, None) for k in
                           ("input_dir", "output_dir", "threads", "root",
                            "branch_length", "samples_map")},
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

    # ---------- 主流程|main run ----------

    def run(self) -> int:
        """端到端运行|Run end-to-end. Returns exit code 0/1."""
        start_time = datetime.now()
        self.logger.warning(
            "waster 硬编码 K=9 需约 64GB 内存,务必在计算节点运行|"
            "waster hardcodes K=9 and needs ~64GB memory; run on a compute node")
        if not self.prepare_inputs():
            return 1
        if not self.run_waster():
            self.logger.error("waster 建树失败,终止|waster failed; abort")
            return 1
        if not self.run_branch_length():
            self.logger.error("枝长计算失败|branch-length step failed")
            return 1
        self.write_software_versions(start_time=start_time)
        self.logger.info("genome2tree 流程完成|genome2tree completed")
        return 0
