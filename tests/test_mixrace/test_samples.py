"""mixrace 样本发现新旧命名双兼容测试|sample discovery dual-naming tests.

契约:命名迁移(95d775d)后,discover_samples 必须同时识别新下划线风格
(_1_clean.fq.gz)与旧点号风格(_1.clean.fq.gz,含 nohost 变体)——旧产物目录
(如迁移前产出的 01_qc/、02_host_filter/)重跑/续跑时不得发现 0 样本。
|Contract: after the naming migration, discover_samples must recognize BOTH the
new underscore style and the legacy dot style (incl. nohost variants), so pre-
migration output dirs still resume instead of discovering 0 samples.
"""
from pathlib import Path

import pytest

from biopytools.mixrace.samples import discover_samples


def _make_pair(d: Path, sample: str, r1_suffix: str, r2_suffix: str):
    (d / f"{sample}{r1_suffix}").write_text("@r\nACGT\n+\nI\n", encoding="utf-8")
    (d / f"{sample}{r2_suffix}").write_text("@r\nACGT\n+\nI\n", encoding="utf-8")


class TestDiscoverSamplesLegacyNaming:
    def test_old_dot_clean_style(self, tmp_path):
        # 旧点号风格 clean(迁移前 01_qc 产物)|legacy dot-style clean
        _make_pair(tmp_path, "15", "_1.clean.fq.gz", "_2.clean.fq.gz")
        _make_pair(tmp_path, "Pb16", "_1.clean.fq.gz", "_2.clean.fq.gz")
        samples = discover_samples(str(tmp_path))
        assert sorted(s["sample"] for s in samples) == ["15", "Pb16"]

    def test_old_dot_r1r2_style(self, tmp_path):
        _make_pair(tmp_path, "S1", "_R1.clean.fq.gz", "_R2.clean.fq.gz")
        assert discover_samples(str(tmp_path))[0]["sample"] == "S1"

    def test_old_dot_nohost_style(self, tmp_path):
        # 旧点号风格 nohost(迁移前寄主剔除产物)|legacy dot-style nohost
        _make_pair(tmp_path, "15", "_1.nohost.fq.gz", "_2.nohost.fq.gz")
        assert discover_samples(str(tmp_path))[0]["sample"] == "15"

    def test_new_underscore_style_still_works(self, tmp_path):
        # 新下划线风格不受影响|new underscore style unaffected
        _make_pair(tmp_path, "15", "_1_clean.fq.gz", "_2_clean.fq.gz")
        assert discover_samples(str(tmp_path))[0]["sample"] == "15"

    def test_mixed_styles_same_sample_no_duplicate(self, tmp_path):
        # 新旧并存同一样本不重复计数|same sample in both styles counts once
        _make_pair(tmp_path, "15", "_1.clean.fq.gz", "_2.clean.fq.gz")
        _make_pair(tmp_path, "15", "_1_clean.fq.gz", "_2_clean.fq.gz")
        assert len(discover_samples(str(tmp_path))) == 1

    def test_empty_dir_returns_empty(self, tmp_path):
        assert discover_samples(str(tmp_path)) == []
