"""mixrace reads统计守卫测试|host_filter stats guard tests.

契约:BAM 未解析(空串)时不得下发注定失败的 samtools 命令——mapped/覆盖广度记 0
+WARNING,统计表照常写出(优雅降级);BAM 在场时 coverage 正常执行。
|Contract: unresolved (empty) BAM must not fire doomed samtools calls — mapped /
breadth recorded as 0 + WARNING, stats table still written (graceful degradation);
coverage runs normally when BAM is present.
"""
import logging
from pathlib import Path

from biopytools.mixrace.config import MixraceConfig
from biopytools.mixrace.host_filter import pathogen_alignment_stats, _read_field_tsv


class _FakeRunner:
    """记录命令的假 runner|command-recording fake runner."""

    def __init__(self):
        self.logger = logging.getLogger("test_mixrace_host_filter")
        self.commands = []

    def run_conda(self, tool_path, args, description="", timeout=None):
        self.commands.append(list(args))
        cov_out = ("#rname\tstartpos\tendpos\tnumreads\tcovbases\n"
                   "chr1\t1\t1000\t500\t900\n")
        return True, cov_out, ""


def _setup(tmp_path: Path, sample: str = "S1", with_mapq: bool = True) -> MixraceConfig:
    """搭 output 目录与阶段缓存表|build output dir + stage cache tables."""
    out = tmp_path / "out"
    (out / "02_host_filter").mkdir(parents=True)
    (out / "04_het_eval").mkdir(parents=True)
    (out / "02_host_filter" / f"{sample}_host_filter.tsv").write_text(
        "field\tvalue\ntotal_reads\t1000\nhost_mapped_reads\t100\nnonhost_reads\t900\n",
        encoding="utf-8")
    if with_mapq:
        (out / "04_het_eval" / f"{sample}_mapq_stats.tsv").write_text(
            "field\tvalue\ntotal_primary_reads\t1000\nmapped_q_reads\t800\n",
            encoding="utf-8")
    return MixraceConfig(output_dir=str(out))


class TestPathogenAlignmentStatsNoBam:
    def test_empty_bam_fires_no_command_and_breadth_zero(self, tmp_path):
        # 空BAM:不下任何samtools命令,breadth记0,表照常写|empty BAM: no samtools
        # call, breadth=0, table still written
        cfg = _setup(tmp_path)
        runner = _FakeRunner()
        stats = pathogen_alignment_stats(cfg, runner, "S1", "", 1000, mean_depth=5.0)
        assert runner.commands == []
        assert stats["breadth_1x"] == 0.0
        saved = _read_field_tsv(tmp_path / "out" / "02_host_filter" / "S1_host_stats.tsv")
        assert saved["breadth_1x"] == 0.0
        assert saved["sample"] == "S1"

    def test_empty_bam_without_mapq_stats_mapped_zero_no_command(self, tmp_path):
        # 空BAM且无mapq缓存:现场计数回退也须跳过(不下doomed view -c)|empty BAM
        # with no cache: live-count fallback also skipped (no doomed view -c)
        cfg = _setup(tmp_path, with_mapq=False)
        runner = _FakeRunner()
        stats = pathogen_alignment_stats(cfg, runner, "S1", "", 1000)
        assert runner.commands == []
        assert stats["pathogen_mapped_reads"] == 0

    def test_real_bam_runs_coverage(self, tmp_path):
        # 对照:非空BAM正常执行coverage并计算breadth|control: real BAM runs
        # coverage and computes breadth
        cfg = _setup(tmp_path)
        runner = _FakeRunner()
        stats = pathogen_alignment_stats(cfg, runner, "S1", "x.bam", 1000)
        assert len(runner.commands) == 1
        assert runner.commands[0][0] == "coverage"
        assert stats["breadth_1x"] == 90.0   # 900/1000
