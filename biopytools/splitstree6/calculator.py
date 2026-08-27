"""splitstree6 流程编排器|SplitsTree6 orchestrator

VCF→距离矩阵(可选) → 构建工作流命令 → Xvfb 下运行 workflow-run → 导出多格式
|VCF→distance matrix (optional) → build command → run under Xvfb → export formats
"""

import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

from ..common.paths import expand_path
from .utils import ModuleLogger, XvfbDisplay, read_vcf_gt_matrix, write_distance_csv


class Splitstree6Calculator:
    """SplitsTree6 建网编排器|SplitsTree6 network/tree orchestrator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        out = Path(config.output_dir)
        self.wf_dir = out / "01_workflow"
        self.info_dir = out / "00_pipeline_info"
        self.logs_dir = out / "99_logs"

    # ---------- 工具路径|tool paths ----------

    @property
    def tools_dir(self) -> str:
        return self.config.splitstree_tools_dir

    def _runtime_env(self) -> dict:
        """用 jpackage 自带运行时调用工具包 jars(SplitsTree6 主程序未用)
        |Use jpackage-bundled runtime for the tools jars (GUI main unused)"""
        runtime_lib = Path(expand_path("~/software/splitstree6/app/SplitsTree/lib/runtime"))
        java_bin = runtime_lib / "bin" / "java"
        if not java_bin.exists():
            # 退化到系统 java|fall back to system java
            for cand in ("~/.local/bin/java", "/usr/bin/java"):
                if os.path.exists(os.path.expanduser(cand)):
                    return {"java": os.path.expanduser(cand), "prefix": []}
            raise FileNotFoundError(
                f"java 未找到|java not found: {java_bin}, ~/.local/bin/java, /usr/bin/java")
        return {"java": str(java_bin), "prefix": []}

    def modpath(self) -> str:
        """tools 包全部 jar 拼 module-path|Concatenate all tool jars"""
        jars_dir = Path(self.tools_dir) / "lib" / "jars"
        jars = sorted(jars_dir.glob("*.jar"))
        return ":".join(str(p.resolve()) for p in jars)

    def prepare_workflow(self) -> str:
        """把官方工作流模板复制进输出目录并定制|Copy template workflow into output dir

        - 追加 RunWorkflow 需要的 node 标题映射保持 Splits;
        - 模板本身即 Loader→Hamming→NeighborNet→ShowSplits 管线。
        |Template is Loader→Hamming→NeighborNet→ShowSplits pipeline.
        """
        src = Path(self.config.workflow)
        dst = self.wf_dir / src.name
        if not dst.exists():
            dst.write_text(src.read_text())
            self.logger.info(f"工作流已就位|Workflow staged: {dst}")
        return str(dst)

    # ---------- VCF 输入|VCF input ----------

    def convert_vcf_if_needed(self) -> Tuple[str, str]:
        """如输入是 VCF,转为 CSV 距离矩阵;否则原样返回
        |Convert VCF to distance CSV; otherwise pass through

        Returns:
            (转换后输入路径, 实际 input format|(converted path, effective format))
        """
        cfg = self.config
        suffix = os.path.basename(cfg.input).lower()
        if suffix.endswith(".vcf") or suffix.endswith(".vcf.gz"):
            import time
            t0 = time.time()
            labels, matrix = read_vcf_gt_matrix(cfg.input)
            csv_out = self.wf_dir / (os.path.basename(suffix.replace(".gz", ""))
                                     .rsplit(".", 1)[0] + ".distances.csv")
            write_distance_csv(labels, matrix, str(csv_out))
            dt = time.time() - t0
            self.logger.info(
                f"VCF 已转距离矩阵|VCF converted to distance matrix: {csv_out} "
                f"(样本|samples={len(labels)}, 变异位点|sites used={len(matrix)}, "
                f"耗时|took {dt:.1f}s)")
            return str(csv_out), "CSV"
        return cfg.input, cfg.input_format or ""

    # ---------- 命令构建|command building ----------

    def build_command(self, workflow: str, data_input: str,
                      fmt: str, export_format: str, out_file: str) -> List[str]:
        """构建 RunWorkflow 命令|Build RunWorkflow command

        注意:SplitsTree 合法导出格式不含 VCF——本模块默认输入 VCF、
        输出 Newick/Nexus/GML 等;若需要"网络→VCF"请另行转换。
        """
        rt = self._runtime_env()
        cmd = [*rt["prefix"], rt["java"], "-Xmx8G", "-Dfile.encoding=UTF8",
               "--module-path", self.modpath(),
               "-m", "splitstreesix/splitstree6.tools.RunWorkflow"]
        cmd += ["-w", workflow]
        cmd += ["-i", data_input]
        if fmt:
            cmd += ["-f", fmt]
        cmd += ["-n", self.config.node_name]
        cmd += ["-e", export_format]
        cmd += ["-o", out_file]
        return cmd

    def run_command(self, cmd: List[str], log_path: str, env: dict) -> bool:
        """执行命令(stderr 落盘)|Run command with stderr captured"""
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        with open(log_path, "w") as fh:
            result = subprocess.run(cmd, shell=False, check=False,
                                    stdout=subprocess.PIPE, stderr=fh,
                                    text=True, env=env)
        if result.returncode != 0:
            tail = ""
            try:
                tail = open(log_path).read()[-500:].strip()
            except OSError:
                pass
            self.logger.error(f"命令失败|Command failed (rc={result.returncode}): {tail}")
            return False
        return True

    # ---------- 主流程|main run ----------

    def run(self) -> int:
        """端到端建网/建树|Run end-to-end. Returns exit code."""
        start_time = datetime.now()
        config = self.config

        try:
            config.validate()
        except ValueError as e:
            self.logger.error(str(e))
            return 1

        # 0. Xvfb(JFX 必需 DISPLAY)|Xvfb required by JavaFX stack
        # CentOS 6 兼容库(libcrypto.so.10)路径提示,Xvfb 启动需要
        # |CentOS-6 compat lib hints required by Xvfb startup
        lib_hints = [os.path.expanduser("~/tmp/x11libs")]
        xvfb = XvfbDisplay(config.xvfb_path, logger=self.logger)
        started = xvfb.start()
        if not started:
            self.logger.error("Xvfb 启动失败,无法运行 SplitsTree6|Xvfb failed to start")
            return 1

        try:
            workflow = self.prepare_workflow()
            data_input, fmt = self.convert_vcf_if_needed()

            ok_all = True
            for fmt_name in config.export_formats_list:
                ext = {"PlainText": "txt",
                       "NexusWithTaxa": "nexus"}.get(fmt_name, fmt_name.lower())
                out_file = self.wf_dir / f"{config.node_name}.{ext}"
                log_path = self.logs_dir / f"run_{fmt_name}.log"
                cmd = self.build_command(workflow, data_input, fmt, fmt_name,
                                         str(out_file))
                self.logger.info(f"开始步骤|Starting step: 导出|Export {fmt_name}")
                if not self.run_command(cmd, str(log_path), xvfb.env):
                    ok_all = False
                    continue
                size = os.path.getsize(out_file) if os.path.exists(out_file) else 0
                if size == 0:
                    self.logger.warning(
                        f"导出为空|Export produced empty file: {fmt_name}")
                    ok_all = False
                else:
                    self.logger.info(f"完成|Done: {out_file} ({size} bytes)")

            self.write_versions(start_time)
            report = f"{config.input} -> " + \
                     ", ".join(config.export_formats_list)
            self.logger.info(f"输入与格式|input/formats: {report}")
            self.logger.info("splitstree6 完成|splitstree6 completed")
            return 0 if ok_all else 1
        finally:
            xvfb.stop()

    # ---------- 版本记录|versions ----------

    def write_versions(self, start_time=None) -> None:
        """写 software_versions.yml|Write software_versions.yml"""
        import yaml
        versions = {
            "pipeline": {"name": "biopytools splitstree6"},
            "splits_tree_tools_dir": self.tools_dir,
            "workflow": self.config.workflow,
        }
        if start_time is not None:
            end_time = datetime.now()
            versions["execution"] = {
                "start_time": start_time.strftime("%Y-%m-%d %H:%M:%S"),
                "end_time": end_time.strftime("%Y-%m-%d %H:%M:%S"),
                "runtime_seconds": int((end_time - start_time).total_seconds()),
            }
        out = self.info_dir / "software_versions.yml"
        with open(out, "w") as fh:
            yaml.dump(versions, fh, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本信息已生成|software versions written: {out}")
