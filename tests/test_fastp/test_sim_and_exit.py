"""fastp模拟数据逐样本检测 + 失败退出码测试|Per-sample simulated detection + failure exit code tests

根因|Root causes:
1. detect_simulated_data 只读第一个样品, 命中即把整批 quality_threshold 改 0,
   混入的 wgsim 模拟材料污染真实样品的 QC;
2. 42/42 全部失败仍正常返回、退出码 0, 上游(fastq2vcf_gtx)误判成功并写检查点。
"""

import gzip
import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.fastp.data_processing import SampleFinder
from biopytools.fastp.main import FastpProcessor
from biopytools.fastp.processing import FastpCore

_LOGGER = logging.getLogger("test_fastp")


def _write_fastq(path: Path, quals):
    """写一个 FASTQ(gzip), quals 为每条read的质量字符串列表|Write a gzipped FASTQ"""
    with gzip.open(str(path), "wt") as f:
        for i, qual in enumerate(quals):
            f.write(f"@read{i} desc\nACGTACGTAC\n+\n{qual}\n")


def _sim_fastq(path: Path):
    """模拟数据: 质量行全为'!'(wgsim默认)|Simulated: uniform '!' quality"""
    _write_fastq(path, ["!" * 10] * 10)


def _real_fastq(path: Path):
    """真实数据: 质量行字符多样|Real data: varied quality characters"""
    _write_fastq(path, [chr(35 + (i % 20)) * 10 for i in range(10)])


class TestDetectSimulatedSamples:
    """逐样本模拟数据检测|Per-sample simulated data detection"""

    def test_mixed_batch_detects_only_simulated_samples(self, tmp_path):
        sim = tmp_path / "sim_1.fq.gz"
        real = tmp_path / "real_1.fq.gz"
        _sim_fastq(sim)
        _real_fastq(real)
        finder = SampleFinder(SimpleNamespace(), _LOGGER)
        pairs = [("sim", sim, None), ("real", real, None)]

        detected = finder.detect_simulated_samples(pairs)

        assert set(detected) == {"sim"}, "只应标记模拟样品, 不得整批判定"

    def test_all_real_batch_detects_nothing(self, tmp_path):
        real1 = tmp_path / "a_1.fq.gz"
        real2 = tmp_path / "b_1.fq.gz"
        _real_fastq(real1)
        _real_fastq(real2)
        finder = SampleFinder(SimpleNamespace(), _LOGGER)
        pairs = [("a", real1, None), ("b", real2, None)]

        assert finder.detect_simulated_samples(pairs) == []


class TestPerSampleQualityThreshold:
    """逐样品质量阈值: 模拟样品 -q 0, 真实样品保持配置值|Per-sample -q override"""

    @staticmethod
    def _make_core(tmp_path):
        class _RecordingRunner:
            def __init__(self):
                self.commands = []

            def run(self, cmd, description=""):
                self.commands.append([str(c) for c in cmd])
                return True

        cfg = SimpleNamespace(
            fastp_path="fastp", threads=12, quality_threshold=30, min_length=50,
            unqualified_percent=40, n_base_limit=10, single_end=True,
            output_path=tmp_path, report_path=tmp_path / "fastp_reports",
            enable_pair=False, force=False, dry_run=False, seqkit_path="seqkit",
        )
        (tmp_path / "fastp_reports").mkdir()
        runner = _RecordingRunner()
        return FastpCore(cfg, _LOGGER, runner), runner, cfg

    def test_simulated_sample_gets_q0_real_keeps_config(self, tmp_path):
        core, runner, cfg = self._make_core(tmp_path)
        sim_r1 = tmp_path / "sim_raw.fq.gz"
        real_r1 = tmp_path / "real_raw.fq.gz"
        sim_r1.write_bytes(b"x")
        real_r1.write_bytes(b"x")

        core.process_sample("sim", sim_r1, None, quality_threshold=0)
        core.process_sample("real", real_r1, None)

        q_values = [c[c.index("-q") + 1] for c in runner.commands]
        assert q_values == ["0", "30"]
        assert cfg.quality_threshold == 30, "全局配置不得被模拟样品改写"


class TestFailureExitCode:
    """失败必须以非零退出码向上传递|Failures must propagate via non-zero exit"""

    @staticmethod
    def _make_processor(tmp_path):
        raw = tmp_path / "raw"
        raw.mkdir()
        (raw / "s1_1.fq.gz").write_bytes(b"x")
        (raw / "s1_2.fq.gz").write_bytes(b"x")
        return FastpProcessor(input_dir=str(raw), output_dir=str(tmp_path / "clean"))

    def test_all_failed_exits_nonzero(self, tmp_path, monkeypatch):
        proc = self._make_processor(tmp_path)
        monkeypatch.setattr(proc.fastp_core, "validate_fastp", lambda: True)
        monkeypatch.setattr(proc.fastp_core, "process_sample", lambda *a, **k: False)

        with pytest.raises(SystemExit) as excinfo:
            proc.run_batch_processing()
        assert excinfo.value.code == 1

    def test_all_succeeded_completes_normally(self, tmp_path, monkeypatch):
        proc = self._make_processor(tmp_path)
        monkeypatch.setattr(proc.fastp_core, "validate_fastp", lambda: True)
        monkeypatch.setattr(proc.fastp_core, "process_sample", lambda *a, **k: True)

        # 全部成功时不得抛 SystemExit|No SystemExit when everything succeeds
        proc.run_batch_processing()
