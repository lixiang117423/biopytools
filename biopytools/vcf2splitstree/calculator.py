"""vcf2splitstree 流程编排器|vcf2splitstree orchestrator

VCF → p-distance 距离矩阵 CSV(SplitsTree6 GUI 可直接打开)
|VCF → p-distance distance matrix CSV (directly openable in SplitsTree6 GUI)
"""

import os
import time
from datetime import datetime
from pathlib import Path

from .utils import read_vcf_gt_matrix, write_distance_csv


class Vcf2SplitstreeCalculator:
    """VCF 距离矩阵转换器|VCF distance-matrix converter"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run(self) -> int:
        """端到端转换|Run end-to-end. Returns exit code 0/1."""
        start = datetime.now()
        cfg = self.config

        try:
            cfg.validate()
        except ValueError as e:
            self.logger.error(str(e))
            return 1

        self.logger.info("开始 VCF 距离矩阵转换|Starting VCF distance conversion")
        self.logger.info(f"输入|Input: {cfg.input}")
        self.logger.info(f"输出|Output: {cfg.output_csv}")

        # 断点续传:输出已存在且非空则跳过|resume: skip when output present
        if cfg.output_csv.exists() and cfg.output_csv.stat().st_size > 0:
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: 输出已存在|output "
                f"exists: {cfg.output_csv}")
            return 0

        t0 = time.time()
        try:
            labels, matrix = read_vcf_gt_matrix(cfg.input)
        except ValueError as e:
            self.logger.error(str(e))
            return 1
        dt = time.time() - t0

        n = len(labels)
        m = matrix.shape[0]
        write_distance_csv(labels, matrix, str(cfg.output_csv))
        size = cfg.output_csv.stat().st_size

        self.logger.info(
            f"转换完成|Converted: 样本|samples={n}, 位点|sites={m}, "
            f"耗时|took {dt:.1f}s, 输出|output={cfg.output_csv} ({size} bytes)")
        self.logger.info(
            f"下一步|Next: 在 Mac 上用 SplitsTree6 GUI 打开该 CSV "
            f"(File→Open), 会自动识别距离矩阵并运行 NeighborNet")
        self._write_versions(start, n, m, dt)
        return 0

    def _write_versions(self, start, n, m, dt) -> None:
        """写 software_versions.yml|Write software_versions.yml"""
        import yaml
        end = datetime.now()
        info = {
            "pipeline": {"name": "biopytools vcf2splitstree", "version": "1.0.0"},
            "input": self.config.input,
            "output_csv": str(self.config.output_csv),
            "samples": n,
            "variant_sites": m,
            "conversion_seconds": round(dt, 1),
            "execution": {
                "start_time": start.strftime("%Y-%m-%d %H:%M:%S"),
                "end_time": end.strftime("%Y-%m-%d %H:%M:%S"),
                "runtime_seconds": int((end - start).total_seconds()),
            },
        }
        info_dir = Path(self.config.output_dir) / '00_pipeline_info'
        info_dir.mkdir(parents=True, exist_ok=True)
        out = info_dir / 'software_versions.yml'
        with open(out, "w") as fh:
            yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本信息已生成|software versions written: {out}")
