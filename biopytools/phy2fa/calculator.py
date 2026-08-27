"""phy2fa 流程编排器|phy2fa orchestrator

Phylip → FASTA 转换(断点续传:输出存在即跳过)
|Phylip → FASTA conversion (resume: skip when output exists)
"""

import time
from datetime import datetime
from pathlib import Path

import yaml

from .utils import open_text, parse_header, parse_phy, write_fasta


class Phy2FaCalculator:
    """Phylip→FASTA 转换器|Phylip to FASTA converter"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.info_dir = Path(config.output_dir) / '00_pipeline_info'

    def run(self) -> int:
        """端到端转换|Run end-to-end. Returns exit code 0/1."""
        start = datetime.now()
        cfg = self.config

        try:
            cfg.validate()
        except ValueError as e:
            self.logger.error(str(e))
            return 1

        self.logger.info("开始 Phylip→FASTA 转换|Starting Phylip→FASTA conversion")
        self.logger.info(f"输入|Input: {cfg.input}")
        self.logger.info(f"输出|Output: {cfg.output_fasta}")

        # 断点续传|resume when output present and non-empty
        if cfg.output_fasta.exists() and cfg.output_fasta.stat().st_size > 0:
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: 输出已存在|output "
                f"exists: {cfg.output_fasta}")
            return 0

        t0 = time.time()
        try:
            with open_text(cfg.input) as fh:
                header = fh.readline()
                n_tax, n_char = parse_header(header)
                records = parse_phy(fh, self.logger)
        except ValueError as e:
            self.logger.error(str(e))
            return 1
        dt = time.time() - t0

        # 一致性校验|consistency check against the header
        problems = []
        if len(records) != n_tax:
            problems.append(
                f"样本数不符|taxa count mismatch: 声明|declared={n_tax}, "
                f"实际|found={len(records)}")
        for name, seq in records.items():
            if len(seq) != n_char:
                problems.append(
                    f"样本 '{name}' 序列长度不符|length mismatch for '{name}': "
                    f"声明|declared={n_char}, 实际|actual={len(seq)}")
        if problems:
            self.logger.error("\n".join(problems))
            return 1

        write_fasta(records, str(cfg.output_fasta), cfg.line_width, self.logger)
        size = cfg.output_fasta.stat().st_size
        total_chars = sum(len(s) for s in records.values())
        self.logger.info(
            f"转换完成|Converted: 样本|samples={len(records)}, "
            f"总字符|total chars={total_chars}, "
            f"耗时|took {dt:.1f}s, 大小|size={size} bytes")

        self._write_versions(start, len(records), n_char, dt)
        self.logger.info("phy2fa 完成|phy2fa completed")
        return 0

    def _write_versions(self, start, n_records, seq_len, dt) -> None:
        """写 software_versions.yml|Write software_versions.yml"""
        end = datetime.now()
        info = {
            "pipeline": {"name": "biopytools phy2fa", "version": "1.0.0"},
            "input": self.config.input,
            "output_fasta": str(self.config.output_fasta),
            "n_tax": n_records,
            "seq_length": seq_len,
            "conversion_seconds": round(dt, 1),
            "execution": {
                "start_time": start.strftime("%Y-%m-%d %H:%M:%S"),
                "end_time": end.strftime("%Y-%m-%d %H:%M:%S"),
                "runtime_seconds": int((end - start).total_seconds()),
            },
        }
        self.info_dir.mkdir(parents=True, exist_ok=True)
        out = self.info_dir / 'software_versions.yml'
        with open(out, "w") as fh:
            yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本信息已生成|software versions written: {out}")
