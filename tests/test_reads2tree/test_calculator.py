"""reads2tree calculator 测试|reads2tree calculator tests"""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from biopytools.reads2tree.calculator import Reads2TreeCalculator
from biopytools.reads2tree.config import Reads2TreeConfig
from biopytools.reads2tree.utils import ModuleLogger


def make_calculator(tmp_path, **kwargs):
    """构造计算器|Build a calculator"""
    d = tmp_path / "reads"
    d.mkdir()
    for sample in ("A", "B", "C"):
        (d / f"{sample}_R1.fq").write_text(f"@{sample}_r1\nACGT\n+\n!!!!\n")
        (d / f"{sample}_R2.fq").write_text(f"@{sample}_r2\nTGCA\n+\n!!!!\n")
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    config = Reads2TreeConfig(input_dir=str(d), **kwargs)
    config.validate()
    logger = ModuleLogger(log_file=str(config.logs_dir / "reads2tree.log"),
                          log_level="ERROR").get_logger()
    return Reads2TreeCalculator(config, logger)


class TestPrepareInputs:
    """输入准备|input preparation"""

    def test_paired_reads_merged_into_single_fq(self, tmp_path):
        calc = make_calculator(tmp_path)
        assert calc.prepare_inputs() is True
        # 每个样本一个合并 fq|one merged fq per sample
        for sample in ("A", "B", "C"):
            merged = calc.config.uncompressed_dir / f"{sample}.fq"
            assert merged.exists()
            content = merged.read_text()
            assert f"@{sample}_r1" in content and f"@{sample}_r2" in content

    def test_input_tsv_written(self, tmp_path):
        calc = make_calculator(tmp_path)
        calc.prepare_inputs()
        lines = Path(calc.config.input_tsv).read_text().strip().splitlines()
        assert len(lines) == 3
        for line in lines:
            stem, path = line.split("\t")
            assert path.endswith(f"{stem}.fq")
            assert Path(path).exists()

    def test_samples_map_written(self, tmp_path):
        map_file = tmp_path / "map.tsv"
        map_file.write_text("A\tSpeciesX\nB\tSpeciesX\nC\tSpeciesY\n")
        calc = make_calculator(tmp_path, samples_map=str(map_file))
        calc.prepare_inputs()
        lines = Path(calc.config.samples_map_tsv).read_text().strip().splitlines()
        assert len(lines) == 3
        assert "A\tSpeciesX" in lines

    def test_name_collision_fails(self, tmp_path):
        # 双端样本 A + 单端文件 A.fq → 冲突|paired A + single A.fq → collision
        d = tmp_path / "reads"
        d.mkdir()
        (d / "A_R1.fq").write_text("x")
        (d / "A_R2.fq").write_text("x")
        (d / "A.fq").write_text("x")
        config = Reads2TreeConfig(input_dir=str(d), output_dir=str(tmp_path / "out"))
        logger = ModuleLogger(log_level="ERROR").get_logger()
        calc = Reads2TreeCalculator(config, logger)
        assert calc.prepare_inputs() is False

    def test_merge_mode_uses_bbmerge(self, tmp_path, monkeypatch):
        bb = tmp_path / "bbmerge.sh"
        bb.write_text("#!/bin/sh\n")
        calc = make_calculator(tmp_path, merge=True, bbmerge_path=str(bb))

        def fake_merge(sample, r1s, r2s, dest_dir, bbmerge_path, logger=None):
            dest = Path(dest_dir) / f"{sample}.fq"
            dest.write_text(f"merged_{sample}")
            return str(dest)

        import biopytools.reads2tree.calculator as calc_mod
        monkeypatch.setattr(calc_mod, "merge_paired_reads", fake_merge)
        assert calc.prepare_inputs() is True
        lines = Path(calc.config.input_tsv).read_text().strip().splitlines()
        assert all(Path(l.split("\t")[1]).read_text().startswith("merged_") for l in lines)


class TestWaster:
    """waster 步骤|waster steps"""

    def test_command_construction(self, tmp_path):
        calc = make_calculator(tmp_path)
        cmd = calc.waster_command()
        assert cmd[0] == calc.config.waster_path
        assert cmd[1:3] == ["-i", calc.config.input_tsv]
        assert cmd[3:5] == ["-o", calc.config.species_tree]
        assert "-t" in cmd

    def test_root_args_appended(self, tmp_path):
        calc = make_calculator(tmp_path, root="SpeciesX")
        cmd = calc.waster_command()
        assert cmd[-2:] == ["--root", "SpeciesX"]

    def test_branchlength_command_uses_topology(self, tmp_path):
        calc = make_calculator(tmp_path, branch_length=True)
        cmd = calc.branchlength_command()
        assert "-C" in cmd and "-c" in cmd
        assert calc.config.species_tree in cmd

    @patch.object(Reads2TreeCalculator, "_run_tool", return_value=True)
    def test_run_waster_checkpoint(self, mock_run, tmp_path):
        calc = make_calculator(tmp_path)
        assert calc.run_waster() is True
        mock_run.assert_called_once()
        # 断点:树已存在则跳过|checkpoint: skip when tree exists
        Path(calc.config.species_tree).write_text("(A,B);")
        assert calc.run_waster() is True
        assert mock_run.call_count == 1

    @patch.object(Reads2TreeCalculator, "_run_tool", return_value=True)
    def test_branch_length_noop_unless_enabled(self, mock_run, tmp_path):
        calc = make_calculator(tmp_path)
        assert calc.run_branch_length() is True
        mock_run.assert_not_called()

    def test_version_probe(self, tmp_path):
        calc = make_calculator(tmp_path)
        with patch("biopytools.reads2tree.calculator.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                [], 0, stdout="", stderr="waster\nVersion: v9.9.9\n")
            assert calc._probe_waster_version("waster") == "v9.9.9"


class TestRun:
    """主流程|main run"""

    @patch.object(Reads2TreeCalculator, "run_branch_length", return_value=True)
    @patch.object(Reads2TreeCalculator, "run_waster", return_value=True)
    @patch.object(Reads2TreeCalculator, "write_software_versions")
    def test_run_end_to_end(self, mock_versions, mock_waster, mock_bl, tmp_path):
        calc = make_calculator(tmp_path)
        assert calc.run() == 0
        mock_waster.assert_called_once()

    @patch.object(Reads2TreeCalculator, "run_waster", return_value=False)
    def test_waster_failure_aborts(self, mock_waster, tmp_path):
        calc = make_calculator(tmp_path)
        assert calc.run() == 1
