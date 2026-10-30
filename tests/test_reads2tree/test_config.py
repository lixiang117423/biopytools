"""reads2tree config 测试|reads2tree config tests"""

import pytest

from biopytools.reads2tree.config import Reads2TreeConfig


def make_fastq_dir(tmp_path, pairs=("A", "B"), single=None):
    """构造假 fastq 输入目录|Build a fake fastq input dir"""
    d = tmp_path / "reads"
    d.mkdir()
    for sample in pairs:
        (d / f"{sample}_R1.fq.gz").write_text("r1")
        (d / f"{sample}_R2.fq.gz").write_text("r2")
    if single:
        (d / f"{single}.fq").write_text("s")
    return d


def make_config(tmp_path, **kwargs):
    """构造配置|Build a config"""
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    return Reads2TreeConfig(input_dir=str(tmp_path / "reads"), **kwargs)


class TestScan:
    """输入扫描|input scanning"""

    def test_paired_detected(self, tmp_path):
        make_fastq_dir(tmp_path)
        cfg = make_config(tmp_path)
        assert {p[0] for p in cfg.paired} == {"A", "B"}
        assert cfg.singles == []
        assert cfg.ignored_files == []

    def test_single_end_kept(self, tmp_path):
        make_fastq_dir(tmp_path, single="C")
        cfg = make_config(tmp_path)
        assert len(cfg.singles) == 1

    def test_no_fastq_validation_error(self, tmp_path):
        d = tmp_path / "reads"
        d.mkdir()
        (d / "x.txt").write_text("x")
        cfg = make_config(tmp_path)
        with pytest.raises(ValueError, match="fastq"):
            cfg.validate()

    def test_missing_input_dir_error(self, tmp_path):
        cfg = make_config(tmp_path)
        cfg.input_dir = str(tmp_path / "nope")
        with pytest.raises(ValueError, match="输入目录"):
            cfg.validate()

    def test_tilde_expansion(self, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path))
        d = tmp_path / "reads"
        d.mkdir()
        (d / "A_R1.fq").write_text("x")
        (d / "A_R2.fq").write_text("x")
        cfg = Reads2TreeConfig(input_dir="~/reads", output_dir=str(tmp_path / "out"))
        assert cfg.input_dir.startswith(str(tmp_path))


class TestValidate:
    """参数校验|parameter validation"""

    def test_valid_config_passes(self, tmp_path):
        make_fastq_dir(tmp_path)
        cfg = make_config(tmp_path)
        assert cfg.validate() is True

    def test_merge_needs_bbmerge(self, tmp_path):
        make_fastq_dir(tmp_path)
        cfg = make_config(tmp_path, merge=True,
                          bbmerge_path=str(tmp_path / "no" / "bbmerge.sh"))
        with pytest.raises(ValueError, match="bbmerge"):
            cfg.validate()

    def test_merge_with_existing_bbmerge_passes(self, tmp_path):
        make_fastq_dir(tmp_path)
        bb = tmp_path / "bbmerge.sh"
        bb.write_text("#!/bin/sh\n")
        cfg = make_config(tmp_path, merge=True, bbmerge_path=str(bb))
        assert cfg.validate() is True

    def test_missing_samples_map_rejected(self, tmp_path):
        make_fastq_dir(tmp_path)
        cfg = make_config(tmp_path, samples_map=str(tmp_path / "nope.tsv"))
        with pytest.raises(ValueError, match="映射"):
            cfg.validate()


class TestPaths:
    """产物路径|artifact paths"""

    def test_step_layout(self, tmp_path):
        make_fastq_dir(tmp_path)
        cfg = make_config(tmp_path)
        assert cfg.species_tree.endswith("02_waster/waster_species_tree.nw")
        assert cfg.input_tsv.endswith("01_input/input.tsv")
        assert str(cfg.uncompressed_dir).endswith("01_input/uncompressed")
